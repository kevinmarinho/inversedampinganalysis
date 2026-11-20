"""Reader for legacy MATLAB Data Acquisition Toolbox log files (``*.daq``).

Used for the external control data set (UFRGS portal frame, D1.5): MCC USB-1208FS,
``analoginput`` object, ``LogToDiskMode='Index'``. The format is undocumented; this reader
was reverse-engineered on those files and is validated sample by sample against the
``Leituras.dat`` text export written by the authors' own script (test
``tests/unit/test_mat_daq.py`` + notebook ``fase2/03``). It is NOT a general ``daqread``.

Layout found:

* ASCII signature ``"MATLAB Data Acquisition File."`` followed by binary blocks;
* object descriptions as MATLAB assignment strings (``x.SampleRate=1000.1;...``), one per
  object (``AnalogInput``, ``DaqHwInfo``, ``Channel``...). Before the acquisition they
  describe the configured object; after it, the final state (``InitialTriggerTime``,
  ``SamplesAcquired``);
* data blocks: little-endian header of 48 bytes

  ======  =======  ===========================================
  offset  type     meaning
  ======  =======  ===========================================
  0       int32    block length in bytes (header + payload)
  4       int32    block type (= 1 for data)
  8       int32    header length (= 48)
  12      int32    (unknown, 9 in all files)
  16      float64  time stamp of the block start (POSIX s)
  24      float64  time stamp of the block end (POSIX s)
  32      int64    index of the first sample in the block
  40      int32    (unknown, 0)
  44      int32    (unknown)
  ======  =======  ===========================================

  followed by ``n_samples x n_channels`` native samples (``uint16``), channel-interleaved.

Volts = native * NativeScaling + NativeOffset (per channel, from the ``Channel`` strings).
The sample instants are ``k / SampleRate`` (the actual rate, 1000.1 Hz here).
"""

from __future__ import annotations

import re
import struct
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

SIGNATURE = b"MATLAB Data Acquisition File."
_HDR = 48


@dataclass(frozen=True)
class DaqRecord:
    path: Path
    fs_hz: float
    volts: np.ndarray  # (n_channels, n_samples) [V]
    native: np.ndarray  # (n_channels, n_samples) raw counts
    scaling: tuple[float, ...]
    offset: tuple[float, ...]
    trigger_time: tuple[float, ...]  # datevec [Y, M, D, h, m, s] (empty if absent)
    device: str
    bits: int
    info: dict = field(default_factory=dict)

    @property
    def n_samples(self) -> int:
        return self.volts.shape[1]

    @property
    def t_s(self) -> np.ndarray:
        return np.arange(self.n_samples) / self.fs_hz


def _num(s: str) -> float:
    return float(s.strip().rstrip(";").strip("[]"))


def _assignments(d: bytes) -> list[str]:
    """All printable runs that look like MATLAB assignment strings ``x.Name=...;``."""
    return [m.decode("latin-1") for m in re.findall(rb"x\.[ -~]{3,}", d)]


def _field(text: str, name: str) -> str | None:
    m = re.search(rf"x\.{name}=(.*?);(?=x\.|$)", text)
    return m.group(1) if m else None


def read_mat_daq(path) -> DaqRecord:
    path = Path(path)
    d = path.read_bytes()
    if not d.startswith(SIGNATURE):
        raise ValueError(f"{path.name}: not a MATLAB DAQ log file")
    texts = _assignments(d)
    ai = [t for t in texts if "x.SampleRate=" in t]
    if not ai:
        raise ValueError(f"{path.name}: AnalogInput description not found")
    final = ai[-1]  # state after the acquisition
    fs = _num(_field(final, "SampleRate"))
    trig = _field(final, "InitialTriggerTime") or ""
    trig_t = tuple(float(v) for v in re.findall(r"[-\d.eE+]+", trig))
    hw = next((t for t in texts if "x.AdaptorName=" in t), "")
    device = (_field(hw, "DeviceName") or "").strip("[]'")
    bits = int(_num(_field(hw, "Bits") or "0"))
    chans = {}
    for t in texts:
        if "x.HwChannel=" in t and "x.NativeScaling=" in t:
            idx = int(_num(_field(t, "Index")))
            chans[idx] = (_num(_field(t, "NativeScaling")), _num(_field(t, "NativeOffset")))
    n_ch = len(chans)
    if n_ch == 0:
        raise ValueError(f"{path.name}: no channel description found")
    scaling = tuple(chans[k][0] for k in sorted(chans))
    offset = tuple(chans[k][1] for k in sorted(chans))

    # data blocks: find every (type=1, header=48) candidate, keep those whose payload is a
    # whole number of scans, then require that their sample indices tile 0..N exactly.
    blocks = {}
    for m in re.finditer(re.escape(struct.pack("<ii", 1, _HDR)), d):
        p = m.start() - 4
        if p < 0 or p + _HDR > len(d):
            continue
        length = struct.unpack_from("<i", d, p)[0]
        payload = length - _HDR
        if payload <= 0 or payload % (2 * n_ch) or p + length > len(d):
            continue
        i0 = struct.unpack_from("<q", d, p + 32)[0]
        n = payload // (2 * n_ch)
        raw = np.frombuffer(d, dtype="<u2", count=n * n_ch, offset=p + _HDR)
        blocks.setdefault(int(i0), []).append(raw.reshape(n, n_ch))
    if not blocks:
        raise ValueError(f"{path.name}: no data blocks")
    out, k = [], 0
    for i0 in sorted(blocks):
        if len(blocks[i0]) != 1 or i0 != k:
            raise ValueError(f"{path.name}: data blocks do not tile (index {i0}, expected {k})")
        out.append(blocks[i0][0])
        k += blocks[i0][0].shape[0]
    native = np.concatenate(out).T.astype(np.int64)
    acquired = _field(final, "SamplesAcquired")
    if acquired is not None and int(_num(acquired)) != native.shape[1]:
        raise ValueError(f"{path.name}: {native.shape[1]} samples read, header says {acquired}")
    volts = native * np.array(scaling)[:, None] + np.array(offset)[:, None]
    return DaqRecord(
        path=path,
        fs_hz=fs,
        volts=volts,
        native=native,
        scaling=scaling,
        offset=offset,
        trigger_time=trig_t,
        device=device,
        bits=bits,
        info={"n_blocks": len(out)},
    )
