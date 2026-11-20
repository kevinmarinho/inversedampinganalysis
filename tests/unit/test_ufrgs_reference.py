"""Our port of the authors' processing reproduces their published table (external data)."""

from pathlib import Path

import numpy as np
import pytest

from dampid.external.ufrgs import authors_method, published_table

ROOT = Path(__file__).resolve().parents[2] / "data/external/ufrgs_portico/ExperimentalData"
pytestmark = pytest.mark.skipif(
    not (ROOT / "ResultadosExperimentais.xlsx").exists(), reason="external data absent"
)


def test_reproduces_published_frequencies_and_damping():
    pytest.importorskip("openpyxl")
    tab = published_table(ROOT / "ResultadosExperimentais.xlsx")
    res = [authors_method(ROOT / f"datastream{int(k)}.daq", int(k)) for k in tab[:, 0]]
    f = np.array([r[0][:4] for r in res])  # the table keeps the first four peaks
    z = np.array([r[1][:4] for r in res])
    assert np.max(np.abs(f / tab[:, 1:5] - 1)) < 5e-4
    dz = np.abs(z / tab[:, 5:9] - 1)
    assert np.max(dz[:, 1:]) < 0.05  # modes 2-4
    assert np.median(dz[:, 0]) < 0.05  # mode 1: envelope approximation, see module doc
