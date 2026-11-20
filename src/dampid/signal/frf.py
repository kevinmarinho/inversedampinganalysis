"""FRF estimators with averaging over records (never over magnitude or phase).

H1 = G_fx / G_ff,  H2 = G_xx / G_xf,  gamma^2 = |G_fx|^2 / (G_ff G_xx) = H1 / H2.
Bendat & Piersol, *Random Data*, 4th ed., ch. 9; Ewins, *Modal Testing*, §4.3.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class FRFEstimate:
    f_hz: np.ndarray  # (n_f,)
    H1: np.ndarray  # (n_ch, n_f) complex, units of x per unit of force
    H2: np.ndarray
    coh: np.ndarray  # (n_ch, n_f)
    Gff: np.ndarray  # (n_f,)  sum over records (scaling irrelevant for the ratios)
    Gxx: np.ndarray  # (n_ch, n_f)
    Gfx: np.ndarray  # (n_ch, n_f)
    n_avg: int
    meta: dict = field(default_factory=dict)
    F_rec: np.ndarray | None = None  # (n_rec, n_f) per-record force spectra (optional)
    X_rec: np.ndarray | None = None  # (n_rec, n_ch, n_f) per-record response spectra

    def h_single(self) -> np.ndarray:
        """Per-record FRFs X/F, shape (n_rec, n_ch, n_f) (requires the record spectra)."""
        if self.F_rec is None or self.X_rec is None:
            raise ValueError("record spectra were not kept")
        with np.errstate(divide="ignore", invalid="ignore"):
            return self.X_rec / self.F_rec[:, None, :]


def estimate_frf(f_hz, F, X, meta=None) -> FRFEstimate:
    """F: (n_rec, n_f) force spectra; X: (n_rec, n_ch, n_f) response spectra."""
    F = np.asarray(F)
    X = np.asarray(X)
    if X.ndim == 2:
        X = X[:, None, :]
    Gff = np.sum(np.abs(F) ** 2, axis=0)
    Gxx = np.sum(np.abs(X) ** 2, axis=0)
    Gfx = np.sum(np.conj(F)[:, None, :] * X, axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        H1 = Gfx / Gff
        H2 = Gxx / np.conj(Gfx)
        coh = np.abs(Gfx) ** 2 / (Gff * Gxx)
    return FRFEstimate(
        np.asarray(f_hz), H1, H2, coh, Gff, Gxx, Gfx, F.shape[0], dict(meta or {}), F, X
    )


def useful_band(f_hz, Gff, coh, force_drop_db=20.0, coh_min=0.9, f_min_hz=1.0):
    """D2.3: band where the force auto-spectrum is within `force_drop_db` of its maximum
    (over f >= f_min) AND the coherence is >= coh_min.

    Returns (f_hi_hz, mask): f_hi is the first frequency above the force maximum where the
    force drops below the threshold; mask marks the bins in [f_min, f_hi] with good coherence.
    """
    f = np.asarray(f_hz)
    level = 10 * np.log10(np.maximum(Gff, 1e-300))
    sel = f >= f_min_hz
    ref = np.max(level[sel])
    i_ref = int(np.nonzero(sel)[0][np.argmax(level[sel])])
    below = np.nonzero((level < ref - force_drop_db) & (np.arange(f.size) > i_ref))[0]
    f_hi = float(f[below[0]]) if below.size else float(f[-1])
    mask = sel & (f <= f_hi) & (np.asarray(coh) >= coh_min)
    return f_hi, mask


def impulse_response(H, n_time: int):
    """IRF by inverse real FFT of a one-sided FRF evaluated on the rfft grid of length n_time.
    (Units: [x per force] per second when H is an FRF and the time step is 1/fs.)"""
    return np.fft.irfft(H, n=n_time)
