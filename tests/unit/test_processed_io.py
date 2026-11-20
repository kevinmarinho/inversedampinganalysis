import numpy as np
import pytest

from dampid.io.processed import has_netcdf, load_frf, save_frf
from dampid.signal.frf import estimate_frf


def _est():
    rng = np.random.default_rng(0)
    F = rng.normal(size=(3, 20)) + 1j * rng.normal(size=(3, 20))
    X = rng.normal(size=(3, 2, 20)) + 1j * rng.normal(size=(3, 2, 20))
    return estimate_frf(np.arange(20.0), F, X, meta={"config": {"exp_tau_s": None}, "id": "x"})


@pytest.mark.parametrize("backend", ["npz", "netcdf"])
def test_roundtrip(tmp_path, backend):
    if backend == "netcdf" and not has_netcdf():
        pytest.skip("xarray/netCDF4 not installed in this environment (D2.4)")
    e = _est()
    p = save_frf(e, tmp_path / "frf", backend)
    r = load_frf(p)
    for k in ("f_hz", "H1", "H2", "coh", "Gff", "Gxx", "Gfx", "F_rec", "X_rec"):
        np.testing.assert_allclose(getattr(r, k), getattr(e, k))
    assert r.meta == e.meta and r.n_avg == 3


def test_single_record_frfs():
    e = _est()
    np.testing.assert_allclose(e.h_single()[:, 1, :], e.X_rec[:, 1, :] / e.F_rec)
