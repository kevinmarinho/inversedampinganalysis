"""Verification of the Phase-2 pipeline on synthetic records with known FRF and damping."""

import numpy as np
import pytest

from dampid.signal import frf as F
from dampid.signal import synthetic as syn
from dampid.signal import windows as W
from dampid.signal.pipeline import ProcessingConfig, process

FS = 25600.0
MODES = dict(f_hz=[150.0, 420.0, 820.0], zeta=[0.002, 0.003, 0.004], residues=[8.0, -5.0, 3.0])


def _record(n_s, t_imp, T=5e-3, amp=30.0, second=None, modes=MODES, noise=0.0, rng=None):
    n = int(n_s * FS)
    f = syn.half_sine(FS, n, t_imp, T, amp)
    if second is not None:
        f += syn.half_sine(FS, n, t_imp + second[0], T, amp * second[1])
    a = syn.modal_acceleration(f, FS, **modes)
    if noise:
        a = a + noise * rng.standard_normal(n)
    return f, a


def test_force_window_and_taper():
    w = W.force_window(100, 40, 50, 5)
    assert np.all(w[40:51] == 1) and w[0] == 0 and w[99] == 0
    assert np.all(np.diff(w[34:41]) >= 0) and 0 < w[37] < 1


def test_exponential_pole_correction_roundtrip():
    f, z, tau = 158.0, 7e-4, 0.5
    w = 2 * np.pi * f
    p = -z * w + 1j * w * np.sqrt(1 - z**2) - 1 / tau  # what the window produces
    fm, zm = abs(p) / (2 * np.pi), -p.real / abs(p)
    fc, zc = W.correct_pole_exponential(fm, zm, tau)
    assert fc == pytest.approx(f, rel=1e-12) and zc == pytest.approx(z, rel=1e-9)
    assert zm - z == pytest.approx(1 / (w * tau), rel=1e-3)  # first-order added damping


def test_tau_for_target_decay():
    assert W.tau_for_target_decay(4.0, -70.0, -60.0) is None
    tau = W.tau_for_target_decay(4.0, -25.0, -60.0)
    assert 20 * np.log10(np.exp(-4.0 / tau)) == pytest.approx(-35.0)


def _band(est, drop=20.0):
    f_hi, _ = F.useful_band(est.f_hz, est.Gff, np.ones_like(est.f_hz), drop, 0.0, 5.0)
    return (est.f_hz > 5) & (est.f_hz < f_hi)


def test_h1_recovers_modal_frf_noise_free():
    recs = [_record(5.0, 0.5 + 0.1 * k, T=3e-3) for k in range(3)]
    est = process(recs, FS, ProcessingConfig())
    sel = _band(est)
    ref = syn.modal_accelerance(est.f_hz[sel], **MODES)
    # relative error, with antiresonances (|ref| -> 0) measured against the band's scale
    scale = np.maximum(np.abs(ref), 0.1 * np.median(np.abs(ref)))
    err = np.abs(est.H1[0, sel] - ref) / scale
    assert np.median(err) < 1e-3 and np.max(err) < 0.01
    np.testing.assert_allclose(est.coh[0, sel], 1.0, atol=1e-9)


def test_truncated_record_exponential_window_is_exact_pole_shift():
    steel = dict(f_hz=[158.0], zeta=[7e-4], residues=[5.0])
    recs = [_record(5.0, 0.6, modes=steel)]
    t_end = 5.0 - 0.6
    tau = W.tau_for_target_decay(t_end, natural_decay_db=-25.0, target_db=-60.0)
    est_rect = process(recs, FS, ProcessingConfig())
    est_exp = process(recs, FS, ProcessingConfig(exp_tau_s=tau))
    f = est_exp.f_hz
    sel = (f > 140) & (f < 176)
    # windowed FRF = true FRF evaluated at the complex frequency w - i/tau
    wt = 2 * np.pi * f[sel] - 1j / tau
    wr = 2 * np.pi * 158.0
    ref_shift = -(wt**2) * 5.0 / (wr**2 - wt**2 + 2j * 7e-4 * wr * wt)
    err_exp = np.max(np.abs(est_exp.H1[0, sel] - ref_shift) / np.abs(ref_shift))
    ref_true = syn.modal_accelerance(f[sel], **steel)
    err_rect = np.max(np.abs(est_rect.H1[0, sel] - ref_true) / np.abs(ref_true))
    assert err_exp < 5e-3
    assert err_rect > 10 * err_exp  # truncation (leakage) is what the window removes


def test_double_impact_keeps_h1_valid():
    recs = [_record(4.0, 0.4, second=(0.012, 0.3))]
    est = process(recs, FS, ProcessingConfig())
    assert est.meta["records"][0]["n_pulses_in_segment"] == 2
    sel = _band(est, drop=15.0)
    ref = syn.modal_accelerance(est.f_hz[sel], **MODES)
    err = np.abs(est.H1[0, sel] - ref) / np.abs(ref)
    assert np.median(err) < 2e-3


def test_output_noise_h1_unbiased_h2_above_coherence_below_one():
    rng = np.random.default_rng(3)
    # 5-s records: the 150-Hz mode (zeta 0.2 %) decays > 60 dB, so no truncation bias
    recs = [_record(5.0, 0.3 + 0.01 * k, noise=0.05, rng=rng) for k in range(16)]
    est = process(recs, FS, ProcessingConfig())
    f = est.f_hz
    near = np.abs(f - 150.0) < 1.0
    off = (f > 250) & (f < 350)
    ref = syn.modal_accelerance(f[near], **MODES)
    assert np.median(np.abs(est.H1[0, near] / ref)) == pytest.approx(1.0, abs=0.01)
    assert np.all(np.abs(est.H2[0]) >= np.abs(est.H1[0]) * (1 - 1e-12))
    # noise lowers the coherence most where the response is small (between modes)
    assert np.min(est.coh[0, near]) > 0.99
    assert np.median(est.coh[0, off]) < np.min(est.coh[0, near])


def test_useful_band_cuts_at_force_drop():
    f = np.linspace(0, 1000, 1001)
    Gff = 1.0 / (1 + (f / 200.0) ** 4)  # -20 dB at ~ 200 * 99^(1/4) ~ 631 Hz
    coh = np.where((f > 400) & (f < 410), 0.5, 0.99)
    f_hi, mask = F.useful_band(f, Gff, coh, 20.0, 0.9, 5.0)
    assert f_hi == pytest.approx(200 * 99**0.25, abs=1.5)
    assert not mask[405] and mask[300] and not mask[700]
