"""Closed-form references for uniform Euler–Bernoulli beams and thin sections.

All inputs/outputs in SI. Frequencies are returned in Hz (suffix ``_hz``).

References
----------
* Blevins, R. D., *Formulas for Natural Frequency and Mode Shape* (beam roots beta_n L).
* Timoshenko & Goodier, *Theory of Elasticity* (torsion of rectangular sections).
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

# Roots of the frequency equations (first five), computed numerically by
# ``bending_roots`` and kept here only as a readable cross-check.
_ROOTS_REF = {
    "clamped-free": [1.87510407, 4.69409113, 7.85475744, 10.99554073, 14.13716839],
    "free-free": [4.73004074, 7.85320462, 10.99560784, 14.13716549, 17.27875966],
    "clamped-clamped": [4.73004074, 7.85320462, 10.99560784, 14.13716549, 17.27875966],
}


def bending_roots(bc: str, n: int = 5) -> np.ndarray:
    """Roots beta_n L of the characteristic equation for the given boundary condition.

    clamped-free:      cos(x) cosh(x) + 1 = 0
    free-free / clamped-clamped: cos(x) cosh(x) - 1 = 0  (x > 0, elastic modes only)
    """
    if bc == "clamped-free":
        f = lambda x: np.cos(x) * np.cosh(x) + 1.0  # noqa: E731
    elif bc in ("free-free", "clamped-clamped"):
        f = lambda x: np.cos(x) * np.cosh(x) - 1.0  # noqa: E731
    else:
        raise ValueError(f"unknown boundary condition {bc!r}")
    # Scan for sign changes on a fine grid, then refine each bracket with Brent.
    grid = np.linspace(0.5, (n + 2) * np.pi, 20000 * (n + 2))
    vals = f(grid)
    idx = np.nonzero(np.sign(vals[:-1]) != np.sign(vals[1:]))[0]
    roots = [brentq(f, grid[i], grid[i + 1], xtol=1e-14) for i in idx]
    return np.array(roots[:n])


def bending_frequencies_hz(E: float, I: float, rho: float, A: float, L: float, bc: str, n: int = 5):
    """f_n = (beta_n L)^2 / (2 pi L^2) * sqrt(E I / (rho A))."""
    bl = bending_roots(bc, n)
    return bl**2 / (2 * np.pi * L**2) * np.sqrt(E * I / (rho * A))


def torsion_constant_rectangle(a: float, b: float) -> float:
    """Saint-Venant torsion constant of a solid rectangle a x b (a >= b), series form.

    J_t = a b^3 / 3 * [1 - (192/pi^5)(b/a) sum_{n odd} tanh(n pi a / (2 b)) / n^5].
    For thin sections this tends to a b^3/3 (1 - 0.63 b/a).
    """
    a, b = max(a, b), min(a, b)
    s = sum(np.tanh(n * np.pi * a / (2 * b)) / n**5 for n in range(1, 40, 2))
    return a * b**3 / 3 * (1 - 192 / np.pi**5 * (b / a) * s)


def torsion_frequencies_hz(
    G: float, Jt: float, rho: float, Ip: float, L: float, bc: str, n: int = 5
):
    """Uniform shaft in torsion: f_n = lambda_n / (2 pi L) sqrt(G J_t / (rho I_p)).

    clamped-free: lambda_n = (2n-1) pi/2 ; free-free / clamped-clamped: lambda_n = n pi.
    """
    k = np.arange(1, n + 1)
    lam = (2 * k - 1) * np.pi / 2 if bc == "clamped-free" else k * np.pi
    return lam / (2 * np.pi * L) * np.sqrt(G * Jt / (rho * Ip))


def axial_frequencies_hz(E: float, rho: float, L: float, bc: str, n: int = 5):
    k = np.arange(1, n + 1)
    lam = (2 * k - 1) * np.pi / 2 if bc == "clamped-free" else k * np.pi
    return lam / (2 * np.pi * L) * np.sqrt(E / rho)


def free_free_end_mass_factor(m_end_total: float, m_beam: float) -> float:
    """Rayleigh estimate of f1 / f1(no mass) for a free-free beam with point masses
    *at the two ends*.

    With the 1st free-free mode normalised so that (1/L) int phi^2 dx = 1, the end
    amplitude is |phi(0)| = |phi(L)| = 2, so the modal mass grows by 4 * m_end_total.
    f/f0 ~= 1 / sqrt(1 + 4 m_end_total / m_beam). First-order (shape unchanged).
    """
    return 1.0 / np.sqrt(1.0 + 4.0 * m_end_total / m_beam)


def mode_shape(x_over_L, bc: str, n: int = 1):
    """Euler–Bernoulli bending mode shape n, normalised so that (1/L) int phi^2 dx = 1.

    free-free:    phi = cosh bx + cos bx - s (sinh bx + sin bx),
                  s = (cosh bL - cos bL) / (sinh bL - sin bL)
    clamped-free: phi = cosh bx - cos bx - s (sinh bx - sin bx),
                  s = (cosh bL + cos bL) / (sinh bL + sin bL)
    """
    bl = bending_roots(bc, n)[n - 1]
    xi = np.asarray(x_over_L, dtype=float)

    def raw(z):
        bx = bl * z
        if bc == "free-free":
            s = (np.cosh(bl) - np.cos(bl)) / (np.sinh(bl) - np.sin(bl))
            return np.cosh(bx) + np.cos(bx) - s * (np.sinh(bx) + np.sin(bx))
        if bc == "clamped-free":
            s = (np.cosh(bl) + np.cos(bl)) / (np.sinh(bl) + np.sin(bl))
            return np.cosh(bx) - np.cos(bx) - s * (np.sinh(bx) - np.sin(bx))
        raise ValueError(bc)

    zz = np.linspace(0, 1, 4001)
    norm = np.sqrt(np.trapezoid(raw(zz) ** 2, zz))
    return raw(xi) / norm


def point_mass_frequency_factor(positions_over_L, masses, m_beam: float, bc: str, n: int = 1):
    """First-order (Rayleigh) estimate of f_n / f_n(bare) for point masses on a uniform beam:

        f/f0 ~= 1 / sqrt(1 + sum_i m_i phi_n(x_i)^2 / m_beam),  phi normalised as in mode_shape.

    Neglects rotary inertia of the masses and the change of mode shape (valid for small
    mass ratios; it is an audit-level sanity check, not the Phase-4 model).
    """
    phi = mode_shape(np.asarray(positions_over_L), bc, n)
    return 1.0 / np.sqrt(1.0 + np.sum(np.asarray(masses) * phi**2) / m_beam)


def free_free_with_point_masses(
    L: float, m_beam: float, EI: float, positions_over_L, masses, n_elastic: int = 6
):
    """Free-free Euler–Bernoulli beam carrying point masses, by Rayleigh–Ritz on the bare
    beam's own modes (rigid translation, rigid rotation and `n_elastic` bending modes).

    Basis: psi_0 = 1, psi_1 = sqrt(3)(2 xi - 1) (rotation about the centre), psi_r = mode_shape
    (all normalised to (1/L) int psi^2 dx = 1, so the bare generalised mass is m_beam * I).
    Point masses add sum_i m_i psi(x_i) psi(x_i)^T to the mass matrix (rotary inertia of the
    masses neglected). Returns (f_hz, Phi_at(xi) callable) with Phi mass-normalised in kg^-1/2:
    the accelerance between xi_f and xi_r is
        A(w) = sum_r -w^2 phi_r(xi_f) phi_r(xi_r) / (w_r^2 - w^2 + 2 i zeta w_r w)
    (rigid modes have w_r = 0). Audit/EDA-level model, not the Phase-4 FE.
    """
    pos = np.atleast_1d(np.asarray(positions_over_L, dtype=float))
    ms = np.atleast_1d(np.asarray(masses, dtype=float))

    def basis(xi):
        xi = np.atleast_1d(np.asarray(xi, dtype=float))
        rows = [np.ones_like(xi), np.sqrt(3.0) * (2 * xi - 1)]
        rows += [mode_shape(xi, "free-free", r) for r in range(1, n_elastic + 1)]
        return np.array(rows)  # (n_basis, n_points)

    bl = bending_roots("free-free", n_elastic)
    w_bare = bl**2 / L**2 * np.sqrt(EI / (m_beam / L))
    nb = n_elastic + 2
    M = m_beam * np.eye(nb)
    K = np.zeros((nb, nb))
    K[2:, 2:] = np.diag(m_beam * w_bare**2)
    B = basis(pos)
    M = M + (B * ms) @ B.T
    import scipy.linalg as sla

    w2, Q = sla.eigh(K, M)  # Q^T M Q = I
    w = np.sqrt(np.clip(w2, 0.0, None))

    def phi_at(xi):
        return basis(xi).T @ Q  # (n_points, n_modes), units kg^-1/2

    return w / (2 * np.pi), phi_at


def accelerance(f_hz, f_modes_hz, phi_f, phi_r, zeta=1e-3):
    """A(f) = sum_r -w^2 phi_r(f) phi_r(r) / (w_r^2 - w^2 + 2 i zeta w_r w)  [(m/s^2)/N]."""
    w = 2 * np.pi * np.asarray(f_hz, dtype=float)[:, None]
    wr = 2 * np.pi * np.asarray(f_modes_hz)[None, :]
    with np.errstate(invalid="ignore", divide="ignore"):
        terms = -(w**2) * (phi_f * phi_r)[None, :] / (wr**2 - w**2 + 2j * zeta * wr * w)
    # 0/0 at f = 0 for rigid modes: the limit is the mass-line value +phi_f phi_r
    terms = np.where(np.isnan(terms), (phi_f * phi_r)[None, :] + 0j, terms)
    return terms.sum(axis=1)
