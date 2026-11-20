"""Synthetic impact records from a modal model (verification of the pipeline).

Acceleration at a sensor for a force at one point, viscous modal damping:
    q_r'' + 2 zeta_r w_r q_r' + w_r^2 q_r = f(t),   a(t) = sum_r A_r q_r''(t)
with A_r = phi_r(sensor) phi_r(force) (mass-normalised). Each mode is integrated exactly
for a piecewise-linear input (scipy.signal.lsim, first-order hold), so the FRF of the
continuous system is
    A(w) = sum_r -w^2 A_r / (w_r^2 - w^2 + 2 i zeta_r w_r w).
"""

from __future__ import annotations

import numpy as np
from scipy import signal


def half_sine(fs: float, n: int, t0_s: float, T_s: float, amp: float) -> np.ndarray:
    t = np.arange(n) / fs
    x = (t >= t0_s) & (t <= t0_s + T_s)
    return np.where(x, amp * np.sin(np.pi * (t - t0_s) / T_s), 0.0)


def modal_acceleration(force, fs, f_hz, zeta, residues):
    """Acceleration response (sum of modes) to the sampled force."""
    t = np.arange(len(force)) / fs
    out = np.zeros(len(force))
    for fr, z, A in zip(
        np.atleast_1d(f_hz), np.atleast_1d(zeta), np.atleast_1d(residues), strict=True
    ):
        w = 2 * np.pi * fr
        # states [q, q']; output q'' = f - 2 z w q' - w^2 q
        sys = signal.StateSpace(
            [[0, 1], [-(w**2), -2 * z * w]], [[0], [1]], [[-(w**2), -2 * z * w]], [[1]]
        )
        _, y, _ = signal.lsim(sys, force, t, interp=True)
        out += A * y
    return out


def modal_accelerance(f_eval_hz, f_hz, zeta, residues):
    w = 2 * np.pi * np.asarray(f_eval_hz, dtype=float)[:, None]
    wr = 2 * np.pi * np.atleast_1d(f_hz)[None, :]
    z = np.atleast_1d(zeta)[None, :]
    A = np.atleast_1d(residues)[None, :]
    return np.sum(-(w**2) * A / (wr**2 - w**2 + 2j * z * wr * w), axis=1)
