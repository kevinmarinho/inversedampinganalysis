"""Characterisation of the legacy 12-DOF beam element.

Checks that the numerical reproduction used by the harness equals the textbook
closed-form entries that ``matrixes()`` must produce symbolically (Przemieniecki /
Ewins). If SymPy is installed, also checks reproduction == original symbolic code.
"""
import numpy as np
import pytest

from dampid.audit import legacy_harness as lh

E, A, Iy, Iz, G, L, r, J = 71e9, 1.6e-4, 1.3e-10, 3.4e-8, 26.7e9, 0.0542, 2710.0, 3.42e-8


def element():
    return lh._legacy_element(E, A, Iy, Iz, G, L, r, J)


def test_bending_stiffness_entries():
    k, _ = element()
    # w-bending block (DOFs 3,4,9,10) uses E*Iy
    np.testing.assert_allclose(k[3, 3], 12 * E * Iy / L**3, rtol=1e-12)
    np.testing.assert_allclose(k[3, 4], 6 * E * Iy / L**2, rtol=1e-12)
    np.testing.assert_allclose(k[4, 4], 4 * E * Iy / L, rtol=1e-12)
    np.testing.assert_allclose(k[4, 10], 2 * E * Iy / L, rtol=1e-12)
    # v-bending block uses E*Iz
    np.testing.assert_allclose(k[1, 1], 12 * E * Iz / L**3, rtol=1e-12)
    # axial and torsion
    np.testing.assert_allclose(k[0, 0], E * A / L, rtol=1e-12)
    np.testing.assert_allclose(k[0, 6], -E * A / L, rtol=1e-12)
    np.testing.assert_allclose(k[5, 5], G * J / L, rtol=1e-12)


def test_mass_entries_including_legacy_rotary_term():
    _, m = element()
    # consistent translational mass 156/420 * rho A L, plus the legacy extra term
    # r*Iy*int(Nw Nw^T) which has the SAME shape-function integral (finding I4)
    np.testing.assert_allclose(m[3, 3], (r * A + r * Iy) * 156 * L / 420, rtol=1e-12)
    np.testing.assert_allclose(m[1, 1], (r * A + r * Iz) * 156 * L / 420, rtol=1e-12)
    np.testing.assert_allclose(m[0, 0], r * A * L / 3, rtol=1e-12)
    np.testing.assert_allclose(m[5, 5], r * J * L / 3, rtol=1e-12)


def test_symmetry():
    k, m = element()
    np.testing.assert_allclose(k, k.T, atol=0)
    np.testing.assert_allclose(m, m.T, atol=0)


@pytest.mark.skipif(not lh.sympy_available(), reason="SymPy not installed")
def test_numeric_reproduction_equals_sympy_original():
    mod = lh.load_legacy(force_numeric_matrixes=False)
    ks, ms = mod.matrixes()
    k, m = element()
    np.testing.assert_allclose(np.array(ks(E, A, Iy, Iz, G, L, r, J), float), k, rtol=1e-10, atol=1e-6)
    np.testing.assert_allclose(np.array(ms(E, A, Iy, Iz, G, L, r, J), float), m, rtol=1e-10, atol=1e-12)


def test_correct_rotary_term_is_rayleigh_beam_matrix():
    _, m_none = lh._legacy_element(E, A, Iy, Iz, G, L, r, J, rotary="none")
    _, m_ok = lh._legacy_element(E, A, Iy, Iz, G, L, r, J, rotary="correct")
    # rho I / (30 L) * 36 on the translational diagonal (Przemieniecki)
    np.testing.assert_allclose(m_ok[3, 3] - m_none[3, 3], r * Iy * 36 / (30 * L), rtol=1e-12)
    np.testing.assert_allclose(m_ok[4, 4] - m_none[4, 4], r * Iy * 4 * L / 30, rtol=1e-12)
