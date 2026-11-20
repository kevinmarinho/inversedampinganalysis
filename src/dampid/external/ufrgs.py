"""UFRGS aluminium portal frame (external control data set, D1.5).

Barazzutti, Gomes & Drehmer, Mendeley Data v3 (2023), DOI 10.17632/npcn4c5htt.3, CC BY 4.0.
Output-only data: 3 accelerometers (z), no force channel.

:func:`authors_method` is a Python port of the authors' ``processamento.m`` (sum of the three
channels, FFT peaks above a per-run hand-tuned threshold, spline-refined peak, rectangular
band-pass of +-15 Hz in the frequency domain, upper envelope over 0.5 s from the maximum,
exponential fit -> zeta = lambda / (2 pi f)). It exists only to check that we read the data
exactly as the authors did and to reproduce their published table; it is NOT our EMA method
(Phase 3). Known deviation: MATLAB's ``envelope(..., 'peaks')`` is approximated by a cubic
spline through local maxima, which changes zeta of mode 1 (~36 Hz, ~18 cycles in the 0.5 s
window) by up to ~25 % in some runs; modes 2-4 agree within a few percent.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize_scalar
from scipy.signal import find_peaks

from dampid.io.mat_daq import read_mat_daq

G = 9.81  # value used by the authors
SETPOINT_V = np.array([2.090, 2.211, 2.150])  # acquisition_3ch.m
SENS_V_PER_MS2 = np.array([3.205, 3.205, 3.205]) / G

# per-run peak threshold (fraction of the largest FFT peak) hand-tuned by the authors
_LEVEL = {k: 0.35 for k in range(1, 101)}
_LEVEL.update(
    {50: .2, 51: .8, 52: .12, 53: .68, 54: .25, 55: .6, 56: .2, 57: .85, 58: .45, 59: .1,
     60: .4, 61: .7, 62: .5, 63: .5, 64: .2, 65: .2, 66: .3, 67: .5, 68: .4, 69: .4, 70: .4,
     71: .5, 72: .5, 74: .5, 76: .5, 78: .5, 84: .3, 88: .3, 91: .25, 93: .38, 94: .38,
     95: .2, 96: .2, 98: .2, 99: .55, 100: .35}
)  # fmt: skip


def acceleration(path) -> tuple[np.ndarray, np.ndarray]:
    """(t [s], a [m/s^2] shape (n, 3)) with the authors' conversion (V - setp) / sens."""
    r = read_mat_daq(path)
    return r.t_s, (r.volts.T - SETPOINT_V) / SENS_V_PER_MS2


def authors_method(path, run: int) -> tuple[np.ndarray, np.ndarray]:
    """Port of processamento.m for one run. Returns (f [Hz], zeta [-]) per detected peak."""
    t, a = acceleration(path)
    dt, df = t[1] - t[0], 1.0 / (t[-1] - t[0])
    a = a - a.mean(axis=0)
    s = a.sum(axis=1)
    s = s - s.mean()
    n2 = round(len(s) / 2)
    spec = np.fft.fft(s) / n2
    mag = np.abs(spec[:n2])
    freq = np.linspace(0, 0.5 / dt, n2)
    locs, _ = find_peaks(
        mag, height=_LEVEL[run] * mag.max(), distance=round((20 if run == 54 else 40) / df)
    )
    if run == 54:  # pks(1)=[]; pks(2)=[] in MATLAB removes the 1st and the original 3rd
        locs = np.delete(locs, [0, 2])
    npt = round(15 / df)
    f_out, z_out = [], []
    for L in locs:
        xs = freq[L - npt : L + npt + 1]
        fi = np.linspace(xs[0], xs[-1], 500)
        f0 = fi[np.argmax(CubicSpline(xs, mag[L - npt : L + npt + 1])(fi))]
        band = spec.copy()
        band[: L - npt] = 0  # MATLAB temp(1:locs-npt)=0 (1-based)
        band[L + npt - 1 :] = 0  # MATLAB temp(locs+npt:end)=0
        sf = np.real(np.fft.ifft(band * n2))
        i0 = int(np.argmax(sf))
        seg = sf[i0 : i0 + round(0.5 / dt) + 1]
        ts = t[i0 : i0 + len(seg)] - t[i0]
        pk, _ = find_peaks(seg, distance=2 * round((1 / f0) / dt))
        knots = np.unique(np.r_[0, pk, len(seg) - 1])
        up = CubicSpline(knots, seg[knots])(np.arange(len(seg)))
        up[0] = seg[0]
        lam = minimize_scalar(
            lambda lm, up=up, ts=ts: np.sum((up - up[0] * np.exp(-lm * ts)) ** 2),
            bounds=(0.0, 50.0),
            method="bounded",
        ).x
        f_out.append(f0)
        z_out.append(lam / (2 * np.pi * f0))
    return np.array(f_out), np.array(z_out)


def published_table(path) -> np.ndarray:
    """ResultadosExperimentais.xlsx -> array (run, f1..f4, z1..z4), runs 51-100."""
    import pandas as pd

    x = pd.read_excel(Path(path), header=None)
    return x.iloc[1:51, 1:10].astype(float).to_numpy()
