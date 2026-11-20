"""1st-campaign force channel (finding C9, decision D1.4).

The hammer was wired to two DAQ inputs (researcher, 2026-11-20). Columns ``Input 1 [g]`` and
``Input 2 [N]`` carry the same force with different conversion factors. The transfer
|In1/In2| is flat over 5–300 Hz, so its median is the scale factor between the two copies.
Force = Input 2 (labelled N); where Input 2 is dead, Input 1 divided by the median ratio of
the files where both are alive.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from dampid.eda import quality as q
from dampid.io.catalog import Catalog
from dampid.io.ni_csv import read_ni_csv

BAND_HZ = (5.0, 300.0)


def duplicate_hammer_table(cat: Catalog, raw_dir="data/raw", alive_snr: float = 20.0):
    rows = []
    for e in [x for x in cat.entries if x.campaign == 1]:
        r = read_ni_csv(f"{raw_dir}/{e.path}")
        fs = r.fs_hz
        x = np.nan_to_num(r.data)
        P = q.find_pulses(x[1], fs)
        m = P[q.main_pulse(P)]
        n = min(int(1.5 * fs), x.shape[1] - m.i_start)
        alive = bool(np.max(np.abs(x[2])) > alive_snr * q.robust_sigma(x[2]))
        row = {"id": e.id, "input2_alive": alive, "ratio_in1_in2": np.nan, "span_db": np.nan}
        if alive:
            f, T = q.quick_h1(x[2], x[1], fs, m.i_start, n)
            b = (f > BAND_HZ[0]) & (f < BAND_HZ[1])
            db = 20 * np.log10(np.abs(T[b]))
            row["ratio_in1_in2"] = float(np.median(np.abs(T[b])))
            row["span_db"] = float(np.percentile(db, 95) - np.percentile(db, 5))
        rows.append(row)
    return pd.DataFrame(rows)


def assign_force_channels(cat: Catalog, table: pd.DataFrame) -> float:
    """Set force_channel / force_scale_to_N / force_note on 1st-campaign entries.
    Returns the median In1/In2 ratio used for the files whose Input 2 is dead."""
    t = table.set_index("id")
    r_med = float(t.ratio_in1_in2.dropna().median())
    r_min, r_max = t.ratio_in1_in2.min(), t.ratio_in1_in2.max()
    for e in cat.entries:
        if e.campaign != 1:
            continue
        if bool(t.loc[e.id, "input2_alive"]):
            e.force_channel, e.force_scale_to_N = 2, 1.0
            e.force_note = "Input 2 [N]; Input 1 é cópia do mesmo martelo (D1.4)"
        else:
            e.force_channel, e.force_scale_to_N = 1, 1.0 / r_med
            e.force_note = (
                f"Input 2 sem sinal: força = Input 1 / {r_med:.3f} (razão In1/In2 mediana dos "
                f"outros arquivos, {r_min:.3f}–{r_max:.3f}) (D1.4)"
            )
    return r_med
