"""Reader for the NI DAQExpress / SignalExpress CSV exports used in both campaigns.

File layout (confirmed on all 156 files in Phase 1):

    Timestamp;<date time>;Timestamp;<date time>;...        (one pair per channel)
    Interval;<dt [s]>;Interval;<dt>;...
    Channel name;"<name>";Channel name;"<name>";...
    Unit;"<unit>";Unit;"<unit>";...
    t;value;t;value;...                                      (';' separator, ',' decimal)

Each channel has its own time column. In the 1st campaign the channels can have
different lengths (shorter ones are padded with empty cells -> NaN here).

The values are already in engineering units (the DAQ applied the sensor sensitivity
configured in the task). No unit conversion happens here: see :mod:`dampid.io.units`.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class NIRecord:
    path: Path
    timestamp: str  # as written in the file (dd/mm/yyyy hh:mm:ss)
    dt_s: float
    names: tuple[str, ...]
    units: tuple[str, ...]
    data: np.ndarray  # (n_channels, n_samples), NaN-padded
    n_valid: tuple[int, ...]  # valid (non-NaN) samples per channel
    time_max_rel_dev: float  # max |diff(t) - dt| / dt over all channels
    time_offsets_s: tuple[float, ...] = field(default=())  # first time stamp per channel
    anomalies: tuple[str, ...] = field(default=())  # structural problems found while reading

    @property
    def fs_hz(self) -> float:
        return 1.0 / self.dt_s

    @property
    def n_samples(self) -> int:
        return self.data.shape[1]

    @property
    def duration_s(self) -> float:
        return self.n_samples * self.dt_s

    def channel(self, name_or_index) -> np.ndarray:
        if isinstance(name_or_index, int):
            return self.data[name_or_index]
        return self.data[self.names.index(name_or_index)]


def read_header(path) -> dict:
    with open(path, encoding="utf-8-sig") as fh:
        rows = [next(fh).rstrip("\r\n").split(";") for _ in range(4)]
    stamps = rows[0][1::2]
    intervals = [float(s.replace(",", ".")) for s in rows[1][1::2]]
    names = [s.strip('"') for s in rows[2][1::2]]
    units = [s.strip('"') for s in rows[3][1::2]]
    return {"timestamps": stamps, "intervals": intervals, "names": names, "units": units}


def read_ni_csv(path, check_time: bool = True) -> NIRecord:
    """Read one file. Raises if the channels declare different sampling intervals."""
    path = Path(path)
    head = read_header(path)
    if len(set(head["intervals"])) != 1:
        raise ValueError(f"{path.name}: channels with different intervals {head['intervals']}")
    dt = head["intervals"][0]
    df = pd.read_csv(path, sep=";", decimal=",", skiprows=4, header=None, encoding="utf-8-sig")
    arr = df.to_numpy(dtype=float)
    n_ch = len(head["names"])
    t = arr[:, 0 : 2 * n_ch : 2].T
    x = arr[:, 1 : 2 * n_ch : 2].T
    anomalies: list[str] = []
    dev = 0.0
    offsets = []
    if check_time:
        # A record is the leading block of rows whose time stamps advance by exactly dt.
        # Anything after the first break (one 1st-campaign file has ~2200 rows of a second,
        # overlapping block appended) is reported and discarded, never silently used.
        cut = arr.shape[0]
        for k, tc in enumerate(t):
            ok = np.isfinite(tc)
            last = int(np.nonzero(ok)[0][-1]) + 1 if ok.any() else 0
            offsets.append(float(tc[0]) if ok[0] else float("nan"))
            d = np.abs(np.diff(tc[:last]) - dt) / dt
            bad = np.nonzero(~(d < 1e-3))[0]
            c_k = int(bad[0]) + 1 if bad.size else last
            if bad.size:
                anomalies.append(
                    f"channel {k}: time break after sample {c_k} "
                    f"(t={tc[c_k - 1]:.6g} -> {tc[c_k]:.6g}); {last - c_k} trailing rows discarded"
                )
            cut = min(cut, c_k)
            if c_k > 1:
                dev = max(dev, float(np.max(d[: c_k - 1])))
        x = x[:, :cut]
    n_valid = tuple(int(np.sum(np.isfinite(c))) for c in x)
    return NIRecord(
        path=path,
        timestamp=head["timestamps"][0],
        dt_s=dt,
        names=tuple(head["names"]),
        units=tuple(head["units"]),
        data=x,
        n_valid=n_valid,
        time_max_rel_dev=dev,
        time_offsets_s=tuple(offsets),
        anomalies=tuple(anomalies),
    )
