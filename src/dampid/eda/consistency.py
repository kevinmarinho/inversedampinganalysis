"""Phase-1 consistency checks across repetitions of one test configuration.

Every repetition gives one single-record FRF (``quality.quick_h1``: rectangular window,
segment from 10 ms before the impact, common length per configuration). From the set:

* ensemble H1 and ordinary coherence across repetitions (how repeatable a/F is);
* dominant peaks of the ensemble FRF where the coherence is high;
* per repetition: peak frequency (parabolic interpolation of |H|) near each ensemble peak
  and a log-magnitude deviation from the ensemble H1 inside the force band.

EDA only: the Phase-2 estimator (trigger, windows, H1/H2 with proper averaging) is
decided separately (D2.1–D2.3).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.signal import find_peaks

from dampid.eda import quality as q
from dampid.io.catalog import Entry
from dampid.io.ni_csv import read_ni_csv
from dampid.io.units import to_si


@dataclass
class ConfigFRF:
    config: str
    f_hz: np.ndarray
    ids: list[str]
    F: np.ndarray  # (n_rec, n_f) force spectra
    X: dict[int, np.ndarray]  # channel -> (n_rec, n_f) response spectra [m/s^2]
    clipped: dict[int, np.ndarray]  # channel -> bool per record
    n_seg: int
    excluded: dict[int, list[str]] = field(default_factory=dict)

    def h_single(self, ch: int) -> np.ndarray:
        return self.X[ch] / self.F

    def h1(self, ch: int, use: np.ndarray | None = None):
        """Ensemble H1 and coherence over the records in `use` (default: not clipped)."""
        use = ~self.clipped[ch] if use is None else use
        if not np.any(use):  # every record saturated: describe them all, caller flags it
            use = np.ones_like(use, dtype=bool)
        F, X = self.F[use], self.X[ch][use]
        Gff = np.sum(np.abs(F) ** 2, axis=0)
        Gfx = np.sum(np.conj(F) * X, axis=0)
        return Gfx / Gff, q.ensemble_coherence(F, X)


def load_config(entries: list[Entry], raw_dir="data/raw", max_s: float = 3.7, pre_s=0.01):
    recs, starts = [], []
    for e in entries:
        r = read_ni_csv(f"{raw_dir}/{e.path}")
        P = q.find_pulses(r.data[2], r.fs_hz)
        recs.append((e, r))
        starts.append(P[q.main_pulse(P)].i_start)
    fs = recs[0][1].fs_hz
    avail = min(r.n_samples - s + int(pre_s * fs) for (_, r), s in zip(recs, starts, strict=True))
    n = min(int(max_s * fs), avail - int(0.01 * fs))
    F, X, clip = [], {0: [], 1: []}, {0: [], 1: []}
    f = np.fft.rfftfreq(n, 1 / fs)
    for (e, r), s in zip(recs, starts, strict=True):
        # same segment and baseline as quality.quick_h1, kept as spectra for averaging
        i0 = max(0, s - int(pre_s * fs))
        seg = np.nan_to_num(r.data[2][i0 : i0 + n])
        F.append(np.fft.rfft(seg - np.median(seg[: max(1, int(pre_s * fs) // 2)])))
        for k in (0, 1):
            xs = to_si(np.nan_to_num(r.data[k][i0 : i0 + n]), "g")[0]
            X[k].append(np.fft.rfft(xs))
            fsc = next(c.full_scale for c in e.channels if c.index == k) or 50.0
            clip[k].append(q.clipping(r.data[k], fsc)["n_clipped"] > 0)
    return ConfigFRF(
        config=entries[0].config_id,
        f_hz=f,
        ids=[e.id for e in entries],
        F=np.array(F),
        X={k: np.array(v) for k, v in X.items()},
        clipped={k: np.array(v) for k, v in clip.items()},
        n_seg=n,
    )


def ensemble_peaks(f, H, coh, band, coh_min=0.8, rel_prom=0.1, max_peaks=6):
    """Peaks of |H| inside `band` with coherence >= coh_min (prominence relative to the
    band's maximum, in log10 units)."""
    sel = (f >= band[0]) & (f <= band[1])
    mag = np.log10(np.abs(H[sel]) + 1e-30)
    idx, props = find_peaks(mag, prominence=rel_prom)
    idx = [i for i in idx if coh[sel][i] >= coh_min]
    idx = sorted(idx, key=lambda i: -np.abs(H[sel])[i])[:max_peaks]
    return np.sort(f[sel][idx])


def peak_frequency(f, H, f0, rel_win=0.03):
    """Peak of |H| within f0*(1 ± rel_win), refined by a parabola through log|H|."""
    sel = np.nonzero((f >= f0 * (1 - rel_win)) & (f <= f0 * (1 + rel_win)))[0]
    i = sel[int(np.argmax(np.abs(H[sel])))]
    if 0 < i < len(f) - 1:
        y0, y1, y2 = np.log(np.abs(H[i - 1 : i + 2]))
        den = y0 - 2 * y1 + y2
        d = 0.5 * (y0 - y2) / den if den != 0 else 0.0
        return float(f[i] + d * (f[1] - f[0]))
    return float(f[i])


def logmag_deviation(f, H_i, H_ref, band) -> float:
    """Median over `band` of |20 log10 |H_i / H_ref|| [dB]."""
    sel = (f >= band[0]) & (f <= band[1])
    return float(np.median(np.abs(20 * np.log10(np.abs(H_i[sel]) / np.abs(H_ref[sel])))))
