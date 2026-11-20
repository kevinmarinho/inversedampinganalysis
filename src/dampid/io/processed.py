"""Processed FRF storage (D2.4): NetCDF through xarray when available, else .npz + JSON.

Both backends store the same variables (f_hz, H1, H2, coh, Gff, Gxx, Gfx with complex
values split into real/imag) and the same metadata (JSON string attribute ``meta``), so a
file written by one environment can be regenerated identically in the other.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np

from dampid.signal.frf import FRFEstimate

_COMPLEX = ("H1", "H2", "Gfx")
_REAL = ("coh", "Gxx")


def has_netcdf() -> bool:
    return all(importlib.util.find_spec(m) is not None for m in ("xarray", "netCDF4"))


def _arrays(est: FRFEstimate) -> dict:
    out = {"f_hz": est.f_hz, "Gff": est.Gff}
    for k in _COMPLEX:
        v = getattr(est, k)
        out[f"{k}_re"], out[f"{k}_im"] = np.real(v), np.imag(v)
    for k in _REAL:
        out[k] = np.real(getattr(est, k))
    for k in ("F_rec", "X_rec"):  # per-record spectra, kept for per-record identification
        v = getattr(est, k)
        if v is not None:
            out[f"{k}_re"], out[f"{k}_im"] = np.real(v), np.imag(v)
    return out


def _dims(name: str, v) -> tuple:
    if name.startswith("F_rec"):
        return ("record", "f_hz")
    return {1: ("f_hz",), 2: ("channel", "f_hz"), 3: ("record", "channel", "f_hz")}[v.ndim]


def save_frf(est: FRFEstimate, path, backend: str | None = None) -> Path:
    """Write `est`; backend 'netcdf' | 'npz' | None (auto). Returns the path written."""
    backend = backend or ("netcdf" if has_netcdf() else "npz")
    path = Path(path)
    meta = json.dumps(est.meta, default=str, ensure_ascii=False)
    arr = _arrays(est)
    if backend == "netcdf":
        import xarray as xr

        ds = xr.Dataset(
            {k: (_dims(k, v), v) for k, v in arr.items() if k != "f_hz"},
            coords={"f_hz": arr["f_hz"], "channel": np.arange(est.H1.shape[0])},
            attrs={"meta": meta, "n_avg": est.n_avg},
        )
        path = path.with_suffix(".nc")
        ds.to_netcdf(path)
    else:
        path = path.with_suffix(".npz")
        np.savez_compressed(path, **arr, meta=np.array(meta), n_avg=np.array(est.n_avg))
    return path


def load_frf(path) -> FRFEstimate:
    path = Path(path)
    if path.suffix == ".nc":
        import xarray as xr

        ds = xr.open_dataset(path)
        get = {k: ds[k].values for k in ds.data_vars}
        get["f_hz"] = ds["f_hz"].values
        meta, n_avg = json.loads(ds.attrs["meta"]), int(ds.attrs["n_avg"])
    else:
        z = np.load(path, allow_pickle=False)
        get = {k: z[k] for k in z.files}
        meta, n_avg = json.loads(str(get.pop("meta"))), int(get.pop("n_avg"))
    cplx = {k: get[f"{k}_re"] + 1j * get[f"{k}_im"] for k in _COMPLEX}
    rec = {
        k: (get[f"{k}_re"] + 1j * get[f"{k}_im"]) if f"{k}_re" in get else None
        for k in ("F_rec", "X_rec")
    }
    return FRFEstimate(
        f_hz=get["f_hz"],
        H1=cplx["H1"],
        H2=cplx["H2"],
        coh=get["coh"],
        Gff=get["Gff"],
        Gxx=get["Gxx"],
        Gfx=cplx["Gfx"],
        n_avg=n_avg,
        meta=meta,
        F_rec=rec["F_rec"],
        X_rec=rec["X_rec"],
    )
