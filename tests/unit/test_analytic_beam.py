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
