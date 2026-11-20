"""Tests for dampid.io (reader, units) and dampid.eda.quality (Phase 1)."""

import numpy as np
import pytest
from scipy.optimize import brentq

from dampid.eda import quality as q
from dampid.io.ni_csv import read_ni_csv
from dampid.io.units import G0, full_scale, to_si

FS = 25600.0


def _write_ni_csv(
    path, data, dt=1 / FS, names=("A1", "A2", "F"), units=("g", "g", "N"), extra=None
):
    n = data.shape[1]
    t = np.arange(n) * dt
    with open(path, "w", encoding="utf-8-sig") as fh:
        fh.write(";".join("Timestamp;01/01/2024 10:00:00" for _ in names) + "\n")
        fh.write(";".join(f"Interval;{dt:.6E}".replace(".", ",") for _ in names) + "\n")
        fh.write(";".join(f'Channel name;"{s}"' for s in names) + "\n")
        fh.write(";".join(f'Unit;"{u}"' for u in units) + "\n")
        rows = list(zip(t, *data, strict=False))
        for row in rows:
            cells = []
            for k in range(len(names)):
                cells += [repr(float(row[0])), repr(float(row[1 + k]))]
            fh.write(";".join(c.replace(".", ",") for c in cells) + "\n")
        if extra is not None:  # appended, overlapping block (seen in one exp_5 file)
            for tt in extra:
                fh.write(";".join([repr(float(tt)).replace(".", ","), "0,5"] * len(names)) + "\n")


def test_reader_roundtrip_and_units(tmp_path):
    rng = np.random.default_rng(0)
    data = rng.normal(size=(3, 500))
    p = tmp_path / "r.csv"
    _write_ni_csv(p, data)
    r = read_ni_csv(p)
    assert r.names == ("A1", "A2", "F") and r.units == ("g", "g", "N")
    assert r.fs_hz == pytest.approx(FS, rel=1e-6)
    np.testing.assert_allclose(r.data, data, rtol=1e-12)
    assert r.time_max_rel_dev < 1e-6 and not r.anomalies
    a_si, unit = to_si(r.channel("A1"), "g")
    assert unit == "m/s^2"
    np.testing.assert_allclose(a_si, data[0] * G0)


def test_reader_discards_appended_block(tmp_path):
    data = np.ones((3, 400))
    p = tmp_path / "bad.csv"
    _write_ni_csv(p, data, extra=np.arange(5) / FS + 0.001)
    r = read_ni_csv(p)
    assert r.n_samples == 400
    assert len(r.anomalies) == 3 and "trailing rows discarded" in r.anomalies[0]


def test_full_scale_ni9234():
    assert full_scale(100.0) == pytest.approx(50.0)  # 100 mV/g accelerometer -> 50 g
    assert full_scale(2.25) == pytest.approx(2222.2, rel=1e-4)  # 2.25 mV/N hammer


def _half_sine(T, amp, fs, n, t0):
    t = np.arange(n) / fs
    f = np.where((t >= t0) & (t <= t0 + T), amp * np.sin(np.pi * (t - t0) / T), 0.0)
    return f


def test_find_pulses_single_double_and_ringing():
    n = int(0.5 * FS)
    f1 = _half_sine(5e-3, 100.0, FS, n, 0.1)
    p = q.find_pulses(f1, FS)
    assert len(p) == 1
    assert p[0].peak == pytest.approx(100.0, rel=1e-3)
    # width above 10 % of a half sine: T (1 - 2 asin(0.1)/pi)
    assert p[0].width_s == pytest.approx(5e-3 * (1 - 2 * np.arcsin(0.1) / np.pi), abs=2 / FS)
    assert p[0].impulse == pytest.approx(100.0 * 2 * 5e-3 / np.pi, rel=0.02)
    f2 = f1 + _half_sine(5e-3, 30.0, FS, n, 0.13)
    p2 = q.find_pulses(f2, FS)
    assert len(p2) == 2 and q.main_pulse(p2) == 0
    # a dip that stays above threshold for < merge gap is not a second impact
    f3 = f1.copy()
    i = int(0.1025 * FS)
    f3[i : i + 3] = 0.0
    assert len(q.find_pulses(f3, FS)) == 1


def test_force_bandwidth_half_sine_matches_analytic():
    T = 5e-3
    n = int(1.0 * FS)
    f = _half_sine(T, 100.0, FS, n, 0.2)
    p = q.find_pulses(f, FS)[0]
    bw = q.force_bandwidth(f, FS, p)

    # |F(f)| / F(0) of a half sine: |cos(pi f T)| / |1 - (2 f T)^2|
    def rel_db(x):
        return 20 * np.log10(abs(np.cos(np.pi * x * T)) / abs(1 - (2 * x * T) ** 2))

    f10 = brentq(lambda x: rel_db(x) + 10, 10, 1.4 / T)
    f20 = brentq(lambda x: rel_db(x) + 20, f10, 1.45 / T)
    assert bw["f_10db"] == pytest.approx(f10, abs=1.5)
    assert bw["f_20db"] == pytest.approx(f20, abs=1.5)


def test_clipping_and_response_quality():
    t = np.arange(int(2 * FS)) / FS
    i0 = int(0.5 * FS)
    x = np.where(t >= 0.5, 60 * np.exp(-3 * (t - 0.5)) * np.sin(2 * np.pi * 150 * (t - 0.5)), 0)
    c = q.clipping(np.clip(x, -50.0, 50.0), 50.0)
    assert c["n_clipped"] > 0 and c["peak_over_fs"] == pytest.approx(1.0)
    rq = q.response_quality(x, FS, i0)
    # envelope decays 3 Np/s: end window centred ~1.45 s after impact vs. first 25 ms
    expected = 20 * np.log10(np.exp(-3 * 1.45) / np.exp(-3 * 0.025))
    assert rq.decay_db == pytest.approx(expected, abs=0.6)
    assert rq.pre_rms == 0.0 and rq.remaining_s == pytest.approx(1.5, abs=1e-3)


def test_mass_line_fit_recovers_rigid_body_level():
    f = np.linspace(1, 400, 4000)
    m, f1, fs_ = 0.4, 150.0, 3.0
    w = 2 * np.pi * f
    rigid = 2.0 / m * (1 + 0 * w)  # e.g. transfer between the two ends of a free beam
    susp = -(w**2) / ((2 * np.pi * fs_) ** 2 - w**2) * 0  # included through am2 below
    elastic = -(w**2) * 4.0 / m / ((2 * np.pi * f1) ** 2 - w**2) * (-0.5)
    H = rigid * (1 + (fs_ / f) ** 2) + elastic + susp + 1e-6j
    r = q.mass_line_fit(f, H, (20, 50))
    assert r["a0"] == pytest.approx(2.0 / m, rel=0.01)


def test_ensemble_coherence_is_one_for_exact_repeats():
    rng = np.random.default_rng(1)
    F = rng.normal(size=(4, 64)) + 1j * rng.normal(size=(4, 64))
    H = rng.normal(size=64) + 1j * rng.normal(size=64)
    np.testing.assert_allclose(q.ensemble_coherence(F, F * H), 1.0, rtol=1e-12)
    X = F * H + 0.5 * (rng.normal(size=(4, 64)) + 1j * rng.normal(size=(4, 64)))
    assert np.mean(q.ensemble_coherence(F, X)) < 0.99
