"""Phase-1 survey: every quality metric of every catalogued record in one table.

``survey_record`` reads one file and returns a flat dict; ``survey_catalog`` runs it on
the whole catalogue (≈1 min for 156 files). The table is cached in
``data/processed/fase1_survey.csv`` and is fully regenerable from ``data/raw``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from dampid.eda import quality as q
from dampid.io.catalog import Catalog, Entry
from dampid.io.ni_csv import read_ni_csv


def survey_record(e: Entry, raw_dir="data/raw") -> dict:
    rec = read_ni_csv(Path(raw_dir) / e.path)
    fs = rec.fs_hz
    row: dict = {"id": e.id, "config": e.config_id, "campaign": e.campaign, "rep": e.repetition}
    # Force channel: the hammer channel. In the 1st campaign the hammer was wired to two
    # inputs (D1.4); take the copy whose pulse stands highest above its own noise (one copy
    # is dead in two files). Descriptive only: the scale to N is in the catalogue.
    hammers = [c.index for c in e.channels if c.role == "hammer"]
    snr = [
        np.nanmax(np.abs(rec.data[k])) / max(q.robust_sigma(rec.data[k]), 1e-12) for k in hammers
    ]
    kf = hammers[int(np.argmax(snr))]
    force = rec.data[kf]
    pulses = q.find_pulses(force, fs)
    im = q.main_pulse(pulses)
    p = pulses[im]
    others = [x for j, x in enumerate(pulses) if j != im]
    row.update(
        {
            "force_channel": kf,
            "n_pulses": len(pulses),
            "t_impact_s": p.i_start / fs,
            "F_peak_N": p.peak,
            "F_width_ms": p.width_s * 1e3,
            "F_impulse_Ns": p.impulse,
            "second_pulse_ratio": max((abs(x.peak) / abs(p.peak) for x in others), default=0.0),
            "second_pulse_dt_ms": (
                (others[int(np.argmax([abs(x.peak) for x in others]))].i_peak - p.i_peak) / fs * 1e3
                if others
                else np.nan
            ),
            "pulses_before_main": im,
            "F_noise_sigma": q.robust_sigma(force),
        }
    )
    bw = q.force_bandwidth(force, fs, p)
    row["F_f10db_hz"], row["F_f20db_hz"] = bw["f_10db"], bw["f_20db"]
    fsc = next((c.full_scale for c in e.channels if c.index == kf), None)
    row["F_peak_over_fs"] = abs(p.peak) / fsc if fsc else np.nan
    for c in e.channels:
        if c.role != "accelerometer":
            continue
        x = rec.data[c.index]
        tag = f"a{c.index}"
        cl = q.clipping(x, c.full_scale or 50.0)
        rq = q.response_quality(x, fs, p.i_start)
        row[f"{tag}_n_clipped"] = cl["n_clipped"]
        row[f"{tag}_peak_over_fs"] = cl["peak_over_fs"]
        row[f"{tag}_clip_span_ms"] = (
            (cl["i_last"] - cl["i_first"] + 1) / fs * 1e3 if cl["n_clipped"] else 0.0
        )
        row[f"{tag}_clip_end_after_impact_ms"] = (
            (cl["i_last"] - p.i_start) / fs * 1e3 if cl["n_clipped"] else np.nan
        )
        for k, v in rq.asdict().items():
            row[f"{tag}_{k}"] = v
    row["duration_s"] = rec.duration_s
    row["n_anomalies"] = len(rec.anomalies)
    return row


def survey_catalog(cat: Catalog, raw_dir="data/raw") -> pd.DataFrame:
    return pd.DataFrame([survey_record(e, raw_dir) for e in cat.entries])
