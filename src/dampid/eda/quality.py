"""Per-record quality metrics for impact tests (Phase 1, EDA).

These metrics *describe* a record; they do not process it. The Phase-2 signal pipeline
(trigger, windows, H1/H2) will reuse some of these ideas but is decided and tested
separately (prompt master §9, Fase 2).

Conventions: inputs in the file's engineering units (g, N) or SI — every function is
unit-agnostic unless stated. Times in seconds, frequencies in Hz.

References
----------
* Ewins, D. J., *Modal Testing*, 2nd ed., §3.9 (impact excitation, force spectrum,
  double hits, windows).
* Brandt, A., *Noise and Vibration Analysis*, Wiley 2011, §13 (impact testing).
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np

# ----------------------------------------------------------------------------- noise


def robust_sigma(x: np.ndarray) -> float:
    """Noise standard deviation from the median absolute deviation (1.4826 * MAD)."""
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    return float(1.4826 * np.median(np.abs(x - np.median(x))))


def short_time_rms(x: np.ndarray, fs: float, win_s: float = 0.01) -> np.ndarray:
    """RMS over consecutive non-overlapping windows of `win_s` seconds."""
    n = max(1, int(round(win_s * fs)))
    m = len(x) // n
    seg = np.asarray(x[: m * n], dtype=float).reshape(m, n)
    return np.sqrt(np.mean(seg**2, axis=1))


# ----------------------------------------------------------------------------- force


@dataclass(frozen=True)
class Pulse:
    i_start: int  # first sample above threshold
    i_end: int  # last sample above threshold
    i_peak: int
    peak: float  # signed peak value (sign of the main pulse)
    width_s: float  # time above 10 % of this pulse's peak
    impulse: float  # integral of F over [i_start, i_end] (unit * s)


def find_pulses(
    force: np.ndarray,
    fs: float,
    rel_threshold: float = 0.05,
    noise_k: float = 8.0,
    merge_gap_s: float = 0.5e-3,
) -> list[Pulse]:
    """Detect force pulses (hammer contacts) in a record.

    A pulse is a contiguous run of samples with s*F > thr, where s is the sign of the
    largest |F| sample (hammers normally give a positive pulse) and
    thr = max(rel_threshold * |F|max, noise_k * sigma_noise). Runs separated by less
    than `merge_gap_s` are merged (ringing inside one contact is not a second impact).
    Pulses are returned in time order. The main pulse is the one with the largest peak.
    """
    f = np.nan_to_num(np.asarray(force, dtype=float))
    f = f - np.median(f)
    k = int(np.argmax(np.abs(f)))
    s = 1.0 if f[k] >= 0 else -1.0
    g = s * f
    thr = max(rel_threshold * g[k], noise_k * robust_sigma(f))
    above = g > thr
    if not above.any():
        return []
    edges = np.diff(above.astype(int))
    starts = list(np.nonzero(edges == 1)[0] + 1)
    ends = list(np.nonzero(edges == -1)[0])
    if above[0]:
        starts.insert(0, 0)
    if above[-1]:
        ends.append(len(g) - 1)
    runs = [[a, b] for a, b in zip(starts, ends, strict=True)]
    merged = [runs[0]]
    gap = int(round(merge_gap_s * fs))
    for a, b in runs[1:]:
        if a - merged[-1][1] <= gap:
            merged[-1][1] = b
        else:
            merged.append([a, b])
    out = []
    for a, b in merged:
        seg = g[a : b + 1]
        ip = a + int(np.argmax(seg))
        # width above 10 % of this pulse's own peak, extended outside the run if needed
        lvl = 0.1 * g[ip]
        lo = ip
        while lo > 0 and g[lo - 1] > lvl:
            lo -= 1
        hi = ip
        while hi < len(g) - 1 and g[hi + 1] > lvl:
            hi += 1
        out.append(
            Pulse(
                i_start=int(a),
                i_end=int(b),
                i_peak=int(ip),
                peak=float(s * g[ip]),
                width_s=(hi - lo + 1) / fs,
                impulse=float(s * np.sum(g[lo : hi + 1]) / fs),
            )
        )
    return out


def main_pulse(pulses: list[Pulse]) -> int:
    return int(np.argmax([abs(p.peak) for p in pulses]))


def force_bandwidth(
    force: np.ndarray, fs: float, pulse: Pulse, pad_s: float = 2e-3, nfft_s: float = 1.0
) -> dict:
    """Spectrum of the *main* pulse alone (rectangular force window around it) and the
    frequencies where it first drops 10 dB and 20 dB below its low-frequency level.

    Reference level = mean |F(f)| over 0 < f <= 5 Hz (for a short pulse this equals the
    impulse). Returns f_10db, f_20db [Hz] (NaN if never reached below fs/2).
    """
    f = np.nan_to_num(np.asarray(force, dtype=float))
    f = f - np.median(f)
    lo = max(0, pulse.i_start - int(pad_s * fs))
    hi = min(len(f), pulse.i_end + int(pad_s * fs) + 1)
    w = np.zeros(int(nfft_s * fs))
    seg = f[lo:hi]
    w[: len(seg)] = seg
    F = np.abs(np.fft.rfft(w)) / fs
    fr = np.fft.rfftfreq(len(w), 1 / fs)
    ref = np.mean(F[(fr > 0) & (fr <= 5.0)])
    db = 20 * np.log10(np.maximum(F, 1e-300) / ref)

    def first_below(level):
        idx = np.nonzero((fr > 5.0) & (db <= level))[0]
        return float(fr[idx[0]]) if idx.size else float("nan")

    return {"f_10db": first_below(-10.0), "f_20db": first_below(-20.0), "f": fr, "db": db}


# ----------------------------------------------------------------------------- response


def clipping(x: np.ndarray, full_scale: float, frac: float = 0.98) -> dict:
    """Samples at or beyond `frac` x full scale (ADC/sensor saturation)."""
    a = np.abs(np.nan_to_num(np.asarray(x, dtype=float)))
    idx = np.nonzero(a >= frac * full_scale)[0]
    return {
        "n_clipped": int(idx.size),
        "peak_over_fs": float(a.max() / full_scale),
        "i_first": int(idx[0]) if idx.size else -1,
        "i_last": int(idx[-1]) if idx.size else -1,
    }


@dataclass(frozen=True)
class ResponseQuality:
    pre_rms: float  # RMS before the impact (window ending 2 ms before the main pulse)
    pre_len_s: float
    post_rms: float  # RMS of the first 50 ms after the main pulse
    end_rms: float  # RMS of the last 0.1 s of the record
    min_strms: float  # smallest 10-ms RMS in the record (quietest stretch)
    dc_pre: float  # median before the impact
    decay_db: float  # 20 log10(end_rms / post_rms)
    pre_db: float  # 20 log10(pre_rms / post_rms): vibration already present before impact
    remaining_s: float  # record length after the impact

    def asdict(self):
        return asdict(self)


def response_quality(x: np.ndarray, fs: float, i_impact: int) -> ResponseQuality:
    x = np.nan_to_num(np.asarray(x, dtype=float))
    gap = int(2e-3 * fs)
    pre = x[: max(0, i_impact - gap)]
    post = x[i_impact : i_impact + int(0.05 * fs)]
    end = x[-int(0.1 * fs) :]

    def rms(v):
        return float(np.sqrt(np.mean(v**2))) if v.size else float("nan")

    pre_rms, post_rms, end_rms = rms(pre), rms(post), rms(end)

    def db(a, b):
        with np.errstate(divide="ignore"):
            return float(20 * np.log10(a / b))

    return ResponseQuality(
        pre_rms=pre_rms,
        pre_len_s=pre.size / fs,
        post_rms=post_rms,
        end_rms=end_rms,
        min_strms=float(np.min(short_time_rms(x, fs))),
        dc_pre=float(np.median(pre)) if pre.size else float("nan"),
        decay_db=db(end_rms, post_rms),
        pre_db=db(pre_rms, post_rms) if pre.size else float("nan"),
        remaining_s=(len(x) - i_impact) / fs,
    )


# ----------------------------------------------------------------------------- FRF (EDA only)


def quick_h1(force, resp, fs, i_impact: int, n: int, pre_s: float = 0.01):
    """Single-record FRF X/F for EDA (rectangular window, no averaging, no correction).

    Segment = [i_impact - pre_s*fs, ... + n). NOT the Phase-2 estimator: it only serves to
    compare repetitions and sensors and to read the mass line. Returns (f [Hz], H).
    """
    i0 = max(0, i_impact - int(pre_s * fs))
    f_seg = np.nan_to_num(np.asarray(force, dtype=float)[i0 : i0 + n])
    x_seg = np.nan_to_num(np.asarray(resp, dtype=float)[i0 : i0 + n])
    if f_seg.size < n:
        raise ValueError("record too short for the requested segment")
    F = np.fft.rfft(f_seg - np.median(f_seg[: max(1, int(pre_s * fs) // 2)]))
    X = np.fft.rfft(x_seg)
    return np.fft.rfftfreq(n, 1 / fs), X / F


def ensemble_coherence(F_list, X_list):
    """Ordinary coherence across an ensemble of records (one spectrum per record):
    gamma^2 = |sum conj(F) X|^2 / (sum |F|^2 sum |X|^2). Equals 1 if every record obeys
    X = H F exactly; drops with noise, nonlinearity, or records that are not repeats."""
    F = np.asarray(F_list)
    X = np.asarray(X_list)
    Gfx = np.sum(np.conj(F) * X, axis=0)
    return np.abs(Gfx) ** 2 / (np.sum(np.abs(F) ** 2, axis=0) * np.sum(np.abs(X) ** 2, axis=0))


def mass_line_fit(f_hz, H, band):
    """Fit Re H(f) = a0 + a2 f^2 + am2 f^-2 inside `band` (Hz).

    Between the suspension modes and the 1st elastic mode a free-free structure behaves
    as a rigid body: a0 is the rigid-body accelerance (1/m)(1 + 12 x_f x_r / L^2) for a
    uniform beam (x measured from the centre). The f^2 term absorbs the stiffness line
    of the elastic modes above, the f^-2 term the tail of the suspension modes.
    Returns dict(a0, a2, am2, rel_resid, imag_over_real).
    """
    f = np.asarray(f_hz)
    sel = (f >= band[0]) & (f <= band[1])
    ff, y = f[sel], np.real(H[sel])
    A = np.column_stack([np.ones_like(ff), ff**2, ff**-2.0])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    resid = y - A @ coef
    return {
        "a0": float(coef[0]),
        "a2": float(coef[1]),
        "am2": float(coef[2]),
        "rel_resid": float(np.sqrt(np.mean(resid**2)) / abs(coef[0])),
        "imag_over_real": float(np.median(np.abs(np.imag(H[sel]))) / abs(coef[0])),
    }
