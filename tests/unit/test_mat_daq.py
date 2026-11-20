"""Reader for the MATLAB DAQ log files of the external UFRGS data set (D1.5).

Skipped when the external data are not present (they are not in git)."""

from pathlib import Path

import numpy as np
import pytest

from dampid.io.mat_daq import read_mat_daq

ROOT = Path(__file__).resolve().parents[2] / "data/external/ufrgs_portico/ExperimentalData"
SETP = np.array([2.090, 2.211, 2.150])  # acquisition_3ch.m
SENS = np.array([3.205, 3.205, 3.205]) / 9.81

pytestmark = pytest.mark.skipif(not (ROOT / "Leituras.dat").exists(), reason="external data absent")


def _leituras():
    lines = (ROOT / "Leituras.dat").read_text(encoding="latin-1").splitlines()
    starts = [i for i, s in enumerate(lines) if s.startswith("Dados gravados")] + [len(lines)]
    out = []
    for a, b in zip(starts[:-1], starts[1:], strict=True):
        rows = [s.split() for s in lines[a + 4 : b] if s.strip()]
        out.append(np.array(rows, dtype=float))
    return out


def test_header_and_size():
    r = read_mat_daq(ROOT / "datastream100.daq")
    assert r.device == "USB-1208FS" and r.bits == 12
    assert r.volts.shape == (3, 10001)
    assert r.fs_hz == pytest.approx(1000.1, rel=1e-6)
    assert r.trigger_time[:5] == (2022, 4, 8, 16, 39)
    assert r.native.min() >= 0 and r.native.max() <= 4095  # 12-bit counts


def test_all_runs_match_text_export():
    """Every run 51-100 equals one Leituras.dat block to the printed precision."""
    blocks = [b for b in _leituras() if b.shape[0] == 10001]
    used = set()
    for k in range(51, 101):
        r = read_mat_daq(ROOT / f"datastream{k}.daq")
        a = (r.volts.T - SETP) / SENS
        hits = [j for j, b in enumerate(blocks) if np.max(np.abs(b[:, 1:4] - a)) < 1e-4]
        assert len(hits) == 1, k
        j = hits[0]
        assert j not in used
        used.add(j)
        assert np.max(np.abs(blocks[j][:, 0] - r.t_s)) <= 5.1e-5  # half a unit of %12.4e
