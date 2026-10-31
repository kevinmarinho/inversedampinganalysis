"""Reference (correct) formulas used as the *yardstick* in the Phase-0 audit.

These are deliberately small and textbook-level. They are NOT the Phase-2 signal
pipeline nor the Phase-4/6 FE and damping modules; those will be written, tested and
decided on in their own phases. Here they only answer "is the legacy result right?".

References
----------
* Ewins, D. J., *Modal Testing: Theory, Practice and Application*, 2nd ed. (FRF, H1,
  half-power, modal damping matrix).
* Bendat & Piersol, *Random Data* (H1 / coherence).
"""

from __future__ import annotations

import numpy as np
import scipy.linalg as sla

G0 = 9.80665  # standard gravity [m/s^2], used only in explicit unit conversions


# ----------------------------------------------------------------------------- damping
def modal_damping_matrix(M: np.ndarray, Phi: np.ndarray, omega_rad: np.ndarray, zeta: np.ndarray):
    """Viscous modal damping, Phi mass-normalised (Phi^T M Phi = I):

        C = M Phi diag(2 zeta_r omega_r) Phi^T M

    so that Phi^T C Phi = diag(2 zeta_r omega_r)  (Ewins, sec. 2.3).
    """
    return M @ Phi @ np.diag(2.0 * zeta * omega_rad) @ Phi.T @ M


def effective_modal_damping(K: np.ndarray, M: np.ndarray, C: np.ndarray):
    """Natural frequencies [rad/s] and damping ratios from the *damped* eigenproblem.

    State-space eigenvalues lambda = -zeta omega +/- i omega sqrt(1-zeta^2), so
    omega_r = |lambda| and zeta_r = -Re(lambda)/|lambda|. Works for any C
    (proportional or not). Returns arrays sorted by omega, one entry per mode.
    """
    n = K.shape[0]
    Minv = np.linalg.inv(M)
    A = np.block([[np.zeros((n, n)), np.eye(n)], [-Minv @ K, -Minv @ C]])
    lam = sla.eigvals(A)
    lam = lam[np.imag(lam) > 0]  # one of each conjugate pair
    w = np.abs(lam)
    z = -np.real(lam) / w
    order = np.argsort(w)
    return w[order], z[order]


# ----------------------------------------------------------------------------- SDOF
def sdof_accelerance(f_hz: np.ndarray, fn_hz: float, zeta: float, m: float = 1.0):
    """A(w) = -w^2 / (m (wn^2 - w^2 + 2 i zeta wn w))  [m/s^2/N], e^{+iwt} convention."""
    w = 2 * np.pi * np.asarray(f_hz)
    wn = 2 * np.pi * fn_hz
    return -(w**2) / (m * (wn**2 - w**2 + 2j * zeta * wn * w))


def half_power_zeta(f_hz: np.ndarray, mag: np.ndarray, f_lo: float, f_hi: float):
    """Half-power (-3 dB) damping estimate on |H| inside [f_lo, f_hi], with linear
    interpolation of the crossing points. zeta = (f2 - f1) / (2 fn)."""
    sel = (f_hz >= f_lo) & (f_hz <= f_hi)
    f, h = f_hz[sel], mag[sel]
    ip = int(np.argmax(h))
    hp = h[ip] / np.sqrt(2)
    # left crossing
    i = ip
    while i > 0 and h[i] > hp:
        i -= 1
    f1 = np.interp(hp, [h[i], h[i + 1]], [f[i], f[i + 1]])
    j = ip
    while j < len(h) - 1 and h[j] > hp:
        j += 1
    f2 = np.interp(hp, [h[j], h[j - 1]], [f[j], f[j - 1]])
    return (f2 - f1) / (2 * f[ip]), f[ip]


# ----------------------------------------------------------------------------- signals
def read_ni_csv(path):
    """Read an NI SignalExpress-style CSV: 4 header lines, 3 (time, value) pairs.

    Returns dict with 'fs', 'names', 'units', 't' (s) and 'ch' (3 x n array, engineering
    units as written in the header). No unit conversion is applied here.
    """
    import pandas as pd

    with open(path, encoding="utf-8-sig") as fh:
        head = [next(fh).rstrip("\r\n").split(";") for _ in range(4)]
    dt = float(head[1][1].replace(",", "."))
    names = [s.strip('"') for s in head[2][1::2]]
    units = [s.strip('"') for s in head[3][1::2]]
    df = pd.read_csv(path, sep=";", decimal=",", skiprows=4, header=None, encoding="utf-8-sig")
    arr = df.to_numpy(dtype=float)
    return {
        "fs": 1.0 / dt,
        "dt": dt,
        "names": names,
        "units": units,
        "t": arr[:, 0],
        "ch": arr[:, 1::2].T,
    }


def force_onset(force: np.ndarray, fs: float, frac: float = 0.05):
    """Index of the first sample where |F| exceeds `frac` of its peak."""
    a = np.abs(force - np.median(force[: int(0.05 * fs)]))
    return int(np.argmax(a > frac * a.max()))


def count_impacts(force: np.ndarray, fs: float, frac: float = 0.1, min_gap_s: float = 0.01):
    """Number of distinct force pulses above `frac` x peak, separated by >= min_gap_s."""
    from scipy.signal import find_peaks

    base = np.median(force[: int(0.05 * fs)])
    a = force - base
    peaks, _ = find_peaks(a, height=frac * a.max(), distance=int(min_gap_s * fs))
    return len(peaks), peaks


def simple_h1(force: np.ndarray, resp: np.ndarray, fs: float):
    """H1 = G_fa / G_ff for a *set* of records (rows), rectangular window, one block each.

    force, resp : arrays (n_records, n_samples), same length.
    Returns f [Hz], H1, coherence gamma^2. Audit yardstick only (no windows/trigger).
    """
    F = np.fft.rfft(np.atleast_2d(force), axis=1)
    X = np.fft.rfft(np.atleast_2d(resp), axis=1)
    Gff = np.mean(np.abs(F) ** 2, axis=0)
    Gxx = np.mean(np.abs(X) ** 2, axis=0)
    Gfx = np.mean(np.conj(F) * X, axis=0)
    H1 = Gfx / Gff
    coh = np.abs(Gfx) ** 2 / (Gff * Gxx)
    f = np.fft.rfftfreq(np.atleast_2d(force).shape[1], 1 / fs)
    return f, H1, coh


# ----------------------------------------------------------------------------- DOFs
LEGACY_NODE_DOFS = ("u", "v", "v'", "w", "w'", "theta")


def describe_legacy_dof(global_index: int, n_points: int, length: float):
    """Translate a legacy global DOF index into (node, x [m], dof name)."""
    node, local = divmod(global_index, 6)
    x = length * node / (n_points - 1)
    return node, x, LEGACY_NODE_DOFS[local]


def projected_modal_damping(
    K: np.ndarray, M: np.ndarray, C: np.ndarray, n_modes: int | None = None
):
    """Modal damping ratio of each *undamped* mode r by projection:

        zeta_r = phi_r^T C phi_r / (2 omega_r),   Phi mass-normalised.

    Exact when Phi^T C Phi is diagonal (proportional damping); for non-proportional C
    it is the usual first-order estimate. zeta_r > 1 means the mode is overdamped (it
    no longer oscillates). Unlike the state-space eigenvalues, the result stays attached
    to the undamped mode r, so modes cannot be mislabelled when some become overdamped.
    Returns (omega_r [rad/s], zeta_r, n_overdamped_state_eigs).
    """
    w2, Phi = sla.eigh(K, M)
    w = np.sqrt(np.abs(w2))
    z = np.einsum("ir,ij,jr->r", Phi, C, Phi) / (2 * w)
    n = K.shape[0]
    Minv = np.linalg.inv(M)
    A = np.block([[np.zeros((n, n)), np.eye(n)], [-Minv @ K, -Minv @ C]])
    lam = sla.eigvals(A)
    n_real = int(np.sum(np.abs(lam.imag) < 1e-9 * np.abs(lam)))
    if n_modes is not None:
        w, z = w[:n_modes], z[:n_modes]
    return w, z, n_real


def fit_sdof_zeta_logmag(f_hz, H_meas, fn_model_hz, band, m_free=True):
    """Least-squares fit of zeta for a 1-DOF accelerance model with a FIXED natural
    frequency, on log10|H| inside `band` (Hz). If m_free, the modal mass (amplitude scale)
    is also fitted (it enters log|H| additively, so it is solved in closed form).

    Used only in the audit to show how a frequency mismatch biases the identified zeta.
    Returns (zeta_hat, m_hat).
    """
    from scipy.optimize import minimize_scalar

    sel = (f_hz >= band[0]) & (f_hz <= band[1])
    f, y = f_hz[sel], np.log10(np.abs(H_meas[sel]))

    def resid(logz):
        h = np.log10(np.abs(sdof_accelerance(f, fn_model_hz, 10**logz)))
        off = np.mean(y - h) if m_free else 0.0
        return y - h - off, off

    res = minimize_scalar(
        lambda lz: np.sum(resid(lz)[0] ** 2),
        bounds=(-6, 0),
        method="bounded",
        options={"xatol": 1e-9},
    )
    _, off = resid(res.x)
    return 10**res.x, 10 ** (-off)
