"""Harness to load and run the *legacy* code (``legacy/Biblioteca.py``) unchanged.

Why a harness is needed
-----------------------
1. The legacy code opens data with relative paths (``exp_5/...``, ``AÇO/...``).
   :func:`legacy_workdir` creates a scratch working directory whose sub-folders are
   symlinks to ``data/raw/<folder>`` so the legacy readers work without editing them.
2. ``Biblioteca.matrixes()`` builds the beam element with SymPy. SymPy could not be
   installed in the audit sandbox (PyPI blocked by the network policy). To still run
   every *other* legacy function unmodified, :func:`load_legacy` can replace
   ``matrixes`` with :func:`legacy_matrixes_numeric`, a **line-by-line numerical
   reproduction** of the same integrals (same shape functions, same terms, same
   mistakes). The reproduction is verified against closed-form element entries in
   ``tests/legacy/test_legacy_element.py``. When SymPy *is* available the original
   symbolic function is used, and the two are compared by the same test.

Nothing in ``legacy/`` is ever modified.
"""

from __future__ import annotations

import contextlib
import importlib.util
import os
import sys
import tempfile
import types
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[3]
LEGACY_DIR = REPO_ROOT / "legacy"
RAW_DIR = REPO_ROOT / "data" / "raw"

# Gauss-Legendre with 6 points integrates polynomials up to degree 11 exactly.
# The highest-degree legacy integrand is a product of two cubics (degree 6).
_GAUSS_X, _GAUSS_W = np.polynomial.legendre.leggauss(6)


def sympy_available() -> bool:
    mod = sys.modules.get("sympy")
    if mod is not None and getattr(mod, "__dampid_stub__", False):
        return False
    try:
        return importlib.util.find_spec("sympy") is not None
    except ValueError:  # a module registered without __spec__
        return False


def _hermite(xi: np.ndarray, L: float):
    """Cubic Hermite functions exactly as ``matrixes()`` builds them.

    ``matrixes()`` solves N a = e_i with rows [w(0), w'(0), w(L), w'(L)], which gives
    the classical functions below (xi = x/L).
    Returns (psi, dpsi/dx, d2psi/dx2), each of shape (4, n).
    """
    psi = np.array(
        [
            1 - 3 * xi**2 + 2 * xi**3,
            L * (xi - 2 * xi**2 + xi**3),
            3 * xi**2 - 2 * xi**3,
            L * (-(xi**2) + xi**3),
        ]
    )
    dpsi = np.array(
        [
            (-6 * xi + 6 * xi**2) / L,
            1 - 4 * xi + 3 * xi**2,
            (6 * xi - 6 * xi**2) / L,
            -2 * xi + 3 * xi**2,
        ]
    )
    ddpsi = np.array(
        [
            (-6 + 12 * xi) / L**2,
            (-4 + 6 * xi) / L,
            (6 - 12 * xi) / L**2,
            (-2 + 6 * xi) / L,
        ]
    )
    return psi, dpsi, ddpsi


def _legacy_element(E, A, Iy, Iz, G, L, r, J, rotary: str = "legacy"):
    """Numerical reproduction of the ``k`` and ``m`` expressions of ``matrixes()``.

    Local DOF order (per element, 12): [u1, v1, v1', w1, w1', t1, u2, v2, v2', w2, w2', t2].
    Kept faithfully, including:
      * ``phi = [x/L, 1 - x/L]`` (node order reversed; harmless, see audit notebook 05);
      * torsional stiffness ``G*J`` with ``J = Iy + Iz`` passed by the caller (finding I3);
      * "rotary inertia" ``r*Iz*int(Nv Nv^T)`` / ``r*Iy*int(Nw Nw^T)`` (finding I4).

    ``rotary`` only exists for the audit comparison: "legacy" (default, as the code),
    "correct" (r*I*int(N' N'^T), the Rayleigh-beam term) or "none".
    """
    x = 0.5 * L * (_GAUSS_X + 1.0)
    w = 0.5 * L * _GAUSS_W
    xi = x / L
    phi = np.array([xi, 1 - xi])
    dphi = np.array([np.full_like(xi, 1 / L), np.full_like(xi, -1 / L)])
    psi, dpsi, ddpsi = _hermite(xi, L)
    n = len(x)

    def vec(entries):
        out = np.zeros((12, n))
        for idx, val in entries:
            out[idx] = val
        return out

    Bn = vec([(0, dphi[0]), (6, dphi[1])])
    Bv = vec([(1, ddpsi[0]), (2, ddpsi[1]), (7, ddpsi[2]), (8, ddpsi[3])])
    Bw = vec([(3, ddpsi[0]), (4, ddpsi[1]), (9, ddpsi[2]), (10, ddpsi[3])])
    Bt = vec([(5, dphi[0]), (11, dphi[1])])
    Nu = vec([(0, phi[0]), (6, phi[1])])
    Nv = vec([(1, psi[0]), (2, psi[1]), (7, psi[2]), (8, psi[3])])
    Nw = vec([(3, psi[0]), (4, psi[1]), (9, psi[2]), (10, psi[3])])
    Nt = vec([(5, phi[0]), (11, phi[1])])

    def integ(a, b):
        return (a * w) @ b.T

    k = (
        E * A * integ(Bn, Bn)
        + E * Iy * integ(Bw, Bw)
        + E * Iz * integ(Bv, Bv)
        + G * J * integ(Bt, Bt)
    )
    m = (
        r * A * integ(Nu, Nu)
        + r * A * integ(Nv, Nv)
        + r * A * integ(Nw, Nw)
        + r * J * integ(Nt, Nt)
    )
    if rotary == "legacy":
        m = m + r * Iz * integ(Nv, Nv) + r * Iy * integ(Nw, Nw)
    elif rotary == "correct":
        dNv = vec([(1, dpsi[0]), (2, dpsi[1]), (7, dpsi[2]), (8, dpsi[3])])
        dNw = vec([(3, dpsi[0]), (4, dpsi[1]), (9, dpsi[2]), (10, dpsi[3])])
        m = m + r * Iz * integ(dNv, dNv) + r * Iy * integ(dNw, dNw)
    elif rotary != "none":
        raise ValueError(rotary)
    return k, m


def legacy_matrixes_numeric():
    """Drop-in replacement for ``Biblioteca.matrixes()`` (same call signature)."""

    def local_stiffness(E, A, Iy, Iz, G, L, r, J):
        return _legacy_element(E, A, Iy, Iz, G, L, r, J)[0]

    def local_mass(E, A, Iy, Iz, G, L, r, J):
        return _legacy_element(E, A, Iy, Iz, G, L, r, J)[1]

    return local_stiffness, local_mass


def load_legacy(force_numeric_matrixes: bool | None = None) -> types.ModuleType:
    """Import ``legacy/Biblioteca.py`` as module ``Biblioteca`` without editing it.

    If SymPy is missing, a stub ``sympy`` module is registered so the top-level
    ``import sympy as sp`` succeeds, and ``matrixes`` is replaced by the numerical
    reproduction. The attribute ``Biblioteca.__matrixes_source__`` records which one
    is in use ("sympy" or "numeric-reproduction").
    """
    use_numeric = (
        (not sympy_available()) if force_numeric_matrixes is None else force_numeric_matrixes
    )
    if not sympy_available() and "sympy" not in sys.modules:
        stub = types.ModuleType("sympy")
        stub.__dampid_stub__ = True
        sys.modules["sympy"] = stub
    if str(LEGACY_DIR) not in sys.path:
        sys.path.insert(0, str(LEGACY_DIR))
    spec = importlib.util.spec_from_file_location("Biblioteca", LEGACY_DIR / "Biblioteca.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["Biblioteca"] = mod
    spec.loader.exec_module(mod)
    if use_numeric:
        mod.matrixes = legacy_matrixes_numeric
        mod.__matrixes_source__ = "numeric-reproduction"
    else:
        mod.__matrixes_source__ = "sympy"
    return mod


@contextlib.contextmanager
def legacy_workdir():
    """chdir into a scratch dir where ``exp_5/``, ``AÇO/``... point to ``data/raw``."""
    old = Path.cwd()
    with tempfile.TemporaryDirectory(prefix="legacy_run_") as tmp:
        for folder in RAW_DIR.iterdir():
            if folder.is_dir():
                os.symlink(folder, Path(tmp) / folder.name)
        os.chdir(tmp)
        try:
            yield Path(tmp)
        finally:
            os.chdir(old)


def assemble_legacy(mod, material: str, n_points: int, rotary: str | None = None, props=None):
    """Global K, M exactly as the legacy assembly loop (no boundary conditions).

    ``rotary`` (audit only) bypasses ``mod.matrixes`` and uses the numerical element with
    the chosen rotary-inertia option; ``props`` overrides ``material_prop`` output.
    """
    E, _, r, L, _, _, Iy, Iz, A, G, J = props if props is not None else mod.material_prop(material)
    le = L / (n_points - 1)
    if rotary is None:
        k_fun, m_fun = mod.matrixes()
        k = k_fun(E, A, Iy, Iz, G, le, r, J)
        m = m_fun(E, A, Iy, Iz, G, le, r, J)
    else:
        k, m = _legacy_element(E, A, Iy, Iz, G, le, r, J, rotary=rotary)
    nd = 6 * n_points
    K = np.zeros((nd, nd))
    M = np.zeros((nd, nd))
    for i in range(n_points - 1):
        idx = np.arange(6 * i, 6 * i + 12)
        K[np.ix_(idx, idx)] += k
        M[np.ix_(idx, idx)] += m
    return K, M
