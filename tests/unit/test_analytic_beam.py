"""Verification of the closed-form beam references (Blevins tables)."""

import numpy as np
import pytest

from dampid.analytic import beam


@pytest.mark.parametrize("bc", ["clamped-free", "free-free", "clamped-clamped"])
def test_bending_roots_match_tables(bc):
    np.testing.assert_allclose(beam.bending_roots(bc, 5), beam._ROOTS_REF[bc], rtol=1e-7)


def test_torsion_constant_limits():
    # square: J_t = 0.1406 a^4 (Timoshenko & Goodier)
    assert beam.torsion_constant_rectangle(1.0, 1.0) == pytest.approx(0.1406, rel=2e-3)
    # thin strip tends to a b^3 / 3 (1 - 0.63 b/a)
    a, b = 50e-3, 3e-3
    approx = a * b**3 / 3 * (1 - 0.63 * b / a)
    assert beam.torsion_constant_rectangle(a, b) == pytest.approx(approx, rel=5e-3)


def test_end_mass_factor_limits():
    assert beam.free_free_end_mass_factor(0.0, 1.0) == 1.0
    assert beam.free_free_end_mass_factor(0.25, 1.0) == pytest.approx(1 / np.sqrt(2))


def test_free_free_mode_shape_end_value():
    # normalised 1st free-free mode has |phi| = 2 at the free ends
    assert abs(beam.mode_shape(0.0, "free-free", 1)) == pytest.approx(2.0, rel=1e-3)
    assert abs(beam.mode_shape(1.0, "free-free", 1)) == pytest.approx(2.0, rel=1e-3)
    # 2nd mode is antisymmetric: node at the centre
    assert beam.mode_shape(0.5, "free-free", 2) == pytest.approx(0.0, abs=1e-6)


def test_point_mass_factor_matches_end_formula():
    f1 = beam.point_mass_frequency_factor([0.0, 1.0], [0.01, 0.01], 0.13, "free-free", 1)
    assert f1 == pytest.approx(beam.free_free_end_mass_factor(0.02, 0.13), rel=1e-3)


def test_ritz_bare_beam_reproduces_analytic_and_rigid_body():
    L, m, EI = 0.3, 0.36, 200e9 * 0.05 * 0.003**3 / 12
    f, phi = beam.free_free_with_point_masses(L, m, EI, [0.0], [0.0])
    ref = beam.bending_frequencies_hz(200e9, 0.05 * 0.003**3 / 12, m / (L * 0.05 * 0.003),
                                      0.05 * 0.003, L, "free-free", 3)
    np.testing.assert_allclose(f[:2], 0.0, atol=1e-6)
    np.testing.assert_allclose(f[2:5], ref, rtol=1e-9)
    # rigid-body accelerance between centre and end: (1/m)(1 + 12 x_f x_r / L^2) = 1/m
    ff = np.array([1e-3])
    A = beam.accelerance(ff, f[:2], phi(0.5)[0, :2], phi(1.0)[0, :2])
    # a rigid mode contributes -w^2 phi phi / (0 - w^2) = +phi phi (mass line, in phase with F)
    assert A.real[0] == pytest.approx(1.0 / m, rel=1e-6)
    A_end = beam.accelerance(ff, f[:2], phi(1.0)[0, :2], phi(1.0)[0, :2])
    assert abs(A_end[0]) == pytest.approx(4.0 / m, rel=1e-6)


def test_ritz_end_masses_lower_f1_like_rayleigh():
    L, m, EI = 0.3, 0.13, 70e9 * 0.05 * 0.0032**3 / 12
    f0, _ = beam.free_free_with_point_masses(L, m, EI, [0.0, 1.0], [0.0, 0.0])
    f1, _ = beam.free_free_with_point_masses(L, m, EI, [0.0, 1.0], [0.0075, 0.0075])
    ratio = f1[2] / f0[2]
    assert ratio < 1.0
    # the first-order estimate keeps the bare shape (ignores the rigid-body coupling the
    # masses create), so it only agrees to a few per cent
    assert ratio == pytest.approx(beam.free_free_end_mass_factor(0.015, m), rel=0.03)
    # Ritz convergence: doubling the basis changes f1 by < 0.1 %
    f2, _ = beam.free_free_with_point_masses(L, m, EI, [0.0, 1.0], [0.0075, 0.0075], n_elastic=12)
    assert f2[2] == pytest.approx(f1[2], rel=1e-3)
