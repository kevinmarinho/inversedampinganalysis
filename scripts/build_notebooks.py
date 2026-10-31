"""Build and *execute* Jupyter notebooks from percent-format sources.

Sources live in ``notebooks/<phase>/_src/NN_name.py`` using the jupytext "percent"
format (``# %%`` code cells, ``# %% [markdown]`` text cells), so they diff cleanly in
git and can also be opened directly with jupytext/VS Code. This script executes each
source top to bottom in a fresh namespace, captures printed output and every
matplotlib figure, and writes a standard ``.ipynb`` (nbformat 4.5) with the outputs
embedded, plus the figures as PNG under ``docs/reports/<phase>/figs/``.

It exists because Jupyter itself could not be installed in the audit sandbox
(PyPI blocked). On a normal machine you can simply open the .ipynb and re-run it.

Usage:  python scripts/build_notebooks.py fase0 [NN ...]
"""

from __future__ import annotations

import base64
import contextlib
import io
import json
import os
import sys
import time
import traceback
import uuid
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


def parse_percent(text: str):
    cells, kind, buf = [], None, []

    def flush():
        if kind is None:
            return
        src = "\n".join(buf).strip("\n")
        if kind == "markdown":
            src = "\n".join(
                line[2:] if line.startswith("# ") else line.lstrip("#") for line in src.splitlines()
            )
        if src.strip():
            cells.append((kind, src))

    for line in text.splitlines():
        if line.startswith("# %%"):
            flush()
            kind, buf = ("markdown" if "[markdown]" in line else "code"), []
        else:
            if kind is None:
                kind = "code"
            buf.append(line)
    flush()
    return cells


def _lines(s: str):
    parts = s.splitlines(keepends=True)
    return parts if parts else [""]


def execute(src_path: Path, phase: str):
    stem = src_path.stem
    fig_dir = ROOT / "docs" / "reports" / phase / "figs"
    fig_dir.mkdir(parents=True, exist_ok=True)
    cells = parse_percent(src_path.read_text(encoding="utf-8"))
    ns: dict = {"__name__": "__main__"}
    nb_cells, n_exec, n_fig = [], 0, 0
    t0 = time.time()
    old_cwd = os.getcwd()
    os.chdir(ROOT)
    sys.path[:0] = [str(ROOT / "src")]
    try:
        for kind, src in cells:
            cid = uuid.uuid4().hex[:8]
            if kind == "markdown":
                nb_cells.append(
                    {"cell_type": "markdown", "id": cid, "metadata": {}, "source": _lines(src)}
                )
                continue
            n_exec += 1
            out = io.StringIO()
            outputs, error = [], None
            with contextlib.redirect_stdout(out), contextlib.redirect_stderr(out):
                try:
                    exec(compile(src, f"<{stem} cell {n_exec}>", "exec"), ns)
                except Exception:  # keep going? no: a failing cell is a bug in the notebook
                    error = traceback.format_exc()
            text = out.getvalue()
            if text:
                outputs.append({"output_type": "stream", "name": "stdout", "text": _lines(text)})
            for num in plt.get_fignums():
                fig = plt.figure(num)
                n_fig += 1
                png = io.BytesIO()
                fig.savefig(png, format="png", dpi=110, bbox_inches="tight")
                fig.savefig(fig_dir / f"{stem}_fig{n_fig:02d}.png", dpi=130, bbox_inches="tight")
                outputs.append(
                    {
                        "output_type": "display_data",
                        "data": {
                            "image/png": base64.b64encode(png.getvalue()).decode(),
                            "text/plain": [f"<Figure {stem}_fig{n_fig:02d}>"],
                        },
                        "metadata": {},
                    }
                )
                plt.close(fig)
            if error:
                outputs.append(
                    {
                        "output_type": "error",
                        "ename": "Exception",
                        "evalue": "see traceback",
                        "traceback": error.splitlines(),
                    }
                )
            nb_cells.append(
                {
                    "cell_type": "code",
                    "id": cid,
                    "metadata": {},
                    "execution_count": n_exec,
                    "outputs": outputs,
                    "source": _lines(src),
                }
            )
            if error:
                print(error, file=sys.stderr)
                raise SystemExit(f"[{stem}] cell {n_exec} failed")
    finally:
        os.chdir(old_cwd)
    nb = {
        "cells": nb_cells,
        "metadata": {
            "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
            "language_info": {"name": "python", "version": sys.version.split()[0]},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }
    out_path = src_path.parent.parent / f"{stem}.ipynb"
    out_path.write_text(json.dumps(nb, indent=1, ensure_ascii=False), encoding="utf-8")
    elapsed = time.time() - t0
    rel = out_path.relative_to(ROOT)
    print(f"[{stem}] {n_exec} code cells, {n_fig} figures, {elapsed:.1f} s -> {rel}")


if __name__ == "__main__":
    phase = sys.argv[1] if len(sys.argv) > 1 else "fase0"
    only = set(sys.argv[2:])
    srcs = sorted((ROOT / "notebooks" / phase / "_src").glob("*.py"))
    for s in srcs:
        if only and s.stem[:2] not in only:
            continue
        execute(s, phase)
