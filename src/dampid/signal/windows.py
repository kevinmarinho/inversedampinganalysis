"""Windows for impact testing and the exact correction of the exponential window.

References
----------
* Ewins, D. J., *Modal Testing*, 2nd ed., §3.9.4 (force and exponential windows).
* Brandt, A., *Noise and Vibration Analysis*, Wiley 2011, §13.3 (exponential window,
  pole shift by -1/tau).
"""

from __future__ import annotations

import numpy as np


def force_window(n: int, i_lo: int, i_hi: int, taper_n: int) -> np.ndarray:
    """1 on [i_lo, i_hi], half-cosine tapers of `taper_n` samples outside, 0 elsewhere.

    Keeps the hammer pulse(s) and removes the force-channel noise of the rest of the record.
    """
    w = np.zeros(n)
    lo, hi = max(0, i_lo), min(n - 1, i_hi)
    w[lo : hi + 1] = 1.0
    if taper_n > 0:
        k = np.arange(1, taper_n + 1)
        ramp = 0.5 * (1 + np.cos(np.pi * k / (taper_n + 1)))  # 1 -> 0
        left = lo - k
        ok = left >= 0
        w[left[ok]] = ramp[ok]
        right = hi + k
        ok = right < n
        w[right[ok]] = ramp[ok]
    return w


def exponential_window(n: int, fs: float, tau_s: float, i0: int) -> np.ndarray:
    """w(t) = 1 for t < t0 and exp(-(t - t0)/tau) for t >= t0, t0 = i0/fs.

    Applied to BOTH force and response (with the force zero before t0), the windowed FRF
    equals the true FRF with every pole shifted by -1/tau in its real part — see
    :func:`correct_pole_exponential`.
    """
    t = (np.arange(n) - i0) / fs
    return np.where(t < 0, 1.0, np.exp(-np.clip(t, 0, None) / tau_s))


def tau_for_target_decay(t_end_s: float, natural_decay_db: float, target_db: float = -60.0):
    """tau [s] so that window x natural decay reaches `target_db` at t_end (both in dB,
    negative). Returns None if the record already decays below the target by itself."""
    extra_db = target_db - natural_decay_db  # still needed (negative)
    if extra_db >= 0:
        return None
    nepers = -extra_db / (20 / np.log(10))  # dB -> Np
    return t_end_s / nepers


def correct_pole_exponential(f_meas_hz, zeta_meas, tau_s):
    """Exact removal of the exponential window from an identified pole.

    Measured pole p_m = -zeta_m w_m + i w_m sqrt(1 - zeta_m^2); the window shifted it by
    -1/tau, so the true pole is p = p_m + 1/tau. Returns (f_hz, zeta) of p. First order:
    zeta ~= zeta_m - 1/(w tau).
    """
    w_m = 2 * np.pi * np.asarray(f_meas_hz, dtype=float)
    z_m = np.asarray(zeta_meas, dtype=float)
    p = -z_m * w_m + 1j * w_m * np.sqrt(1 - z_m**2) + 1.0 / tau_s
    w = np.abs(p)
    return w / (2 * np.pi), -np.real(p) / w
