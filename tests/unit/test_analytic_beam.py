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
