# %% [markdown]
# # Fase 0 · Etapa 0.1 — Setup do repositório e integridade
#
# > **Base de dados usada neste notebook:** apenas os **dados do pesquisador** (CSVs da
# > 1ª campanha `exp_5/` e da 2ª campanha livre-livre `AÇO/`, `Al/`, `ACRILICO/`).
# > Nenhuma base externa. Nenhum dado é modificado: só leitura de cabeçalhos e hashes.
#
# **Objetivo.** Registrar o ambiente em que a auditoria rodou, provar que o legado e os
# dados brutos estão íntegros (hashes) e fazer um primeiro inventário dos cabeçalhos.

# %%
import platform
import subprocess
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import pydantic
import scipy

from dampid.audit import legacy_harness as lh
from dampid.audit.reference import read_ni_csv

ROOT = lh.REPO_ROOT
print("Python   :", sys.version.split()[0], "|", platform.platform())
for mod in (np, scipy, pd, matplotlib, pydantic):
    print(f"{mod.__name__:<10}: {mod.__version__}")
print("SymPy disponível?", lh.sympy_available())

# %% [markdown]
# ## 1. Integridade do legado
# O hash SHA-256 de cada arquivo em `legacy/` é recalculado e comparado ao registrado em
# `legacy/README.md`. Se alguém editar o legado por engano, esta célula acusa.

# %%
import hashlib
import re

readme = (ROOT / "legacy" / "README.md").read_text(encoding="utf-8")
registered = dict(re.findall(r"`(\w+\.py)` \| [^|]+\| `([0-9a-f]{64})`", readme))
for name, h in registered.items():
    now = hashlib.sha256((ROOT / "legacy" / name).read_bytes()).hexdigest()
    print(f"{name:<20} {'OK ' if now == h else 'ALTERADO!'} {now[:16]}…")

# %% [markdown]
# ## 2. Integridade dos dados brutos
# `data/raw/MANIFEST.sha256` lista o hash de cada CSV. Os arquivos estão sem permissão de
# escrita (`chmod a-w`), conforme a regra "dados brutos são somente leitura".

# %%
out = subprocess.run(["sha256sum", "-c", "--quiet", "MANIFEST.sha256"], cwd=ROOT / "data/raw",
                     capture_output=True, text=True)
n_files = len((ROOT / "data/raw/MANIFEST.sha256").read_text().splitlines())
print(f"{n_files} arquivos no manifesto ->", "todos íntegros" if out.returncode == 0 else out.stdout)

# %% [markdown]
# ## 3. Inventário dos cabeçalhos (primeiro olhar)
# Cada CSV tem 4 linhas de cabeçalho (`Timestamp`, `Interval`, `Channel name`, `Unit`) e
# três pares (tempo, valor). Aqui só lemos o cabeçalho e o número de amostras. O
# catálogo completo é tarefa da Fase 1.

# %%
rows = []
for p in sorted((ROOT / "data/raw").rglob("*.csv")):
    with open(p, encoding="utf-8-sig") as fh:
        head = [next(fh).rstrip("\r\n").split(";") for _ in range(4)]
        n = sum(1 for _ in fh)
    dt = float(head[1][1].replace(",", "."))
    rows.append({
        "pasta": p.parent.name, "arquivo": p.name, "data": head[0][1][:10],
        "fs [Hz]": round(1 / dt), "amostras": n, "duração [s]": round(n * dt, 3),
        "canais": " | ".join(s.strip('"') for s in head[2][1::2]),
        "unidades": " | ".join(s.strip('"') for s in head[3][1::2]),
    })
inv = pd.DataFrame(rows)
summary = (inv.groupby(["pasta", "data", "fs [Hz]", "duração [s]", "canais", "unidades"])
              .size().rename("n arquivos").reset_index())
with pd.option_context("display.width", 200, "display.max_colwidth", 60):
    print(summary.to_string(index=False))

# %% [markdown]
# **Leitura.** Na 2ª campanha (nov/2024) os canais têm nome: `Acelerômetro 1`,
# `Acelerometro 2`, `Martelo`, com unidades `g, g, N` — o martelo é o **3º** canal. Na 1ª
# campanha (`exp_5/`, jan/2023) os canais se chamam `Input 0/1/2` com unidades
# `g, g, N`: o canal em **N** (força) também é o **3º**. O legado, porém, lê o **2º**
# canal (`r2`, unidade g) como força — isso é o achado **C9**, investigado no notebook 06.

# %%
for f in ["exp_5/Al_engastado_livre_mponta_aponta.csv", "exp_5/Aco_engastado_livre_mponta_aponta.csv"]:
    d = read_ni_csv(ROOT / "data/raw" / f)
    print(f, "| canais:", d["names"], "| unidades:", d["units"])
inv.to_csv(ROOT / "docs/reports/fase0/00_inventario_cabecalhos.csv", index=False)
print("Inventário completo salvo em docs/reports/fase0/00_inventario_cabecalhos.csv")
