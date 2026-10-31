"""Verification of the audit yardsticks against closed-form results."""
import numpy as np
import pytest

from dampid.audit import reference as ref


def test_modal_damping_matrix_recovers_prescribed_zeta():
    rng = np.random.default_rng(0)
    n = 6
    Q = rng.standard_normal((n, n))
    K = Q @ Q.T + n * np.eye(n)
    M = np.diag(rng.uniform(0.5, 2.0, n))
    import scipy.linalg as sla
    w2, Phi = sla.eigh(K, M)
    zeta = np.linspace(0.001, 0.05, n)
    C = ref.modal_damping_matrix(M, Phi, np.sqrt(w2), zeta)
    w, z = ref.effective_modal_damping(K, M, C)
    np.testing.assert_allclose(z, zeta, rtol=1e-8)
    np.testing.assert_allclose(w, np.sqrt(w2), rtol=1e-6)


def test_half_power_on_sdof():
    f = np.arange(0, 400, 0.01)
    H = ref.sdof_accelerance(f, 150.0, 0.01)
    z, fn = ref.half_power_zeta(f, np.abs(H), 100, 200)
    # accelerance peaks at fn / sqrt(1 - 2 zeta^2) ~ 150.015 Hz
    assert fn == pytest.approx(150.0 / np.sqrt(1 - 2 * 0.01**2), abs=0.01)
    assert z == pytest.approx(0.01, rel=0.02)


def test_simple_h1_identity_system():
    rng = np.random.default_rng(1)
    x = rng.standard_normal((4, 1024))
    f, H, coh = ref.simple_h1(x, 2.0 * x, 1000.0)
    np.testing.assert_allclose(H, 2.0, rtol=1e-12)
    np.testing.assert_allclose(coh, 1.0, rtol=1e-12)


def test_describe_legacy_dof():
    assert ref.describe_legacy_dof(34, 6, 0.271)[2] == "w'"
    assert ref.describe_legacy_dof(33, 6, 0.271)[2] == "w"


def test_projected_damping_matches_prescribed():
    rng = np.random.default_rng(2)
    n = 5
    Q = rng.standard_normal((n, n))
    K = Q @ Q.T + n * np.eye(n)
    M = np.diag(rng.uniform(0.5, 2.0, n))
    import scipy.linalg as sla
    w2, Phi = sla.eigh(K, M)
    zeta = np.linspace(0.002, 0.03, n)
    C = ref.modal_damping_matrix(M, Phi, np.sqrt(w2), zeta)
    w, z, n_real = ref.projected_modal_damping(K, M, C)
    np.testing.assert_allclose(z, zeta, rtol=1e-8)
    assert n_real == 0


def test_fit_sdof_zeta_recovers_truth_without_mismatch():
    f = np.arange(100, 200, 0.05)
    H = ref.sdof_accelerance(f, 152.0, 0.002, m=0.15)
    z, m = ref.fit_sdof_zeta_logmag(f, H, 152.0, (140, 165))
    assert z == pytest.approx(0.002, rel=1e-4)
    assert m == pytest.approx(0.15, rel=1e-4)
