"""Data catalogue (``data/catalog.yaml``): one validated entry per raw file.

The catalogue is *generated* from ``configs/data_sources.yaml`` (static, human-edited
metadata) plus what is read from each file (header, sizes, hashes) and, after the
Phase-1 quality pass, a quality summary and a class (ouro / utilizável / descartar).
Positions of sensors and impact, hammer tip and similar fields stay ``None`` with a
TODO until the researcher provides them — they are never guessed.
"""

from __future__ import annotations

import hashlib
import re
from pathlib import Path
from typing import Literal

import yaml
from pydantic import BaseModel, Field

from dampid.io.ni_csv import read_header, read_ni_csv
from dampid.io.units import unit_info

Quantity = Literal["acceleration", "force", "unknown"]
FileClass = Literal["ouro", "utilizável", "descartar", "pendente"]


class Channel(BaseModel):
    index: int
    name: str  # as written in the file header
    unit: str  # as written in the file header
    quantity_header: str  # quantity implied by the header unit
    role: Literal["accelerometer", "hammer", "unknown"]  # from data_sources.yaml
    sensor_model: str | None = None
    sensitivity_mv: float | None = None  # configured in the DAQ task (nominal)
    full_scale: float | None = None  # in header units, from ±5 V / sensitivity
    location: str | None = None  # TODO researcher (position / direction)


class Entry(BaseModel):
    id: str
    path: str  # relative to data/raw
    sha256: str
    campaign: int
    config_id: str
    specimen: str | None
    material: str | None
    boundary: str | None
    hammer_position: str | None = None  # 1st campaign: from file name
    accel_position: str | None = None
    repetition: int | None
    timestamp: str
    fs_hz: float
    n_samples: int
    duration_s: float
    channels: list[Channel]
    anomalies: list[str] = Field(default_factory=list)
    hammer_tip: str | None = None  # TODO researcher
    impact_location: str | None = None  # TODO researcher
    quality: dict = Field(default_factory=dict)
    file_class: FileClass = "pendente"
    class_reasons: list[str] = Field(default_factory=list)


class Catalog(BaseModel):
    generated_by: str
    sources: str
    n_files: int
    entries: list[Entry]

    def by_config(self) -> dict[str, list[Entry]]:
        out: dict[str, list[Entry]] = {}
        for e in self.entries:
            out.setdefault(e.config_id, []).append(e)
        return out


_EXP5 = re.compile(
    r"^(?P<mat>Aco|Al|Acrilico)_(?P<bc>engastado_livre|biengastado)_m(?P<m>\w+?)_a(?P<a>\w+)$"
)
_MAT1 = {"Aco": "ACO", "Al": "AL", "Acrilico": "ACRILICO"}
_BC1 = {"engastado_livre": "clamped-free", "biengastado": "clamped-clamped"}


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _repetition(stem: str) -> int | None:
    m = re.search(r"_(\d+)$", stem)
    return int(m.group(1)) if m else None


def build_entry(path: Path, raw_dir: Path, sources: dict) -> Entry:
    folder = path.parent.name
    fcfg = sources["folders"][folder]
    hw = sources["hardware"]
    rec = read_ni_csv(path)
    head = read_header(path)
    stem = path.stem
    specimen = fcfg.get("specimen")
    boundary = fcfg.get("boundary")
    config_id = fcfg["config_id"]
    hpos = apos = None
    rep = _repetition(stem)
    if config_id == "from_filename":
        m = _EXP5.match(stem)
        if not m:
            raise ValueError(f"unexpected 1st-campaign file name {stem}")
        specimen = _MAT1[m["mat"]]
        boundary = _BC1[m["bc"]]
        hpos, apos = m["m"], m["a"]
        config_id = f"C1_{specimen}_{'CF' if boundary == 'clamped-free' else 'CC'}_m{hpos}_a{apos}"
        rep = 1
    spec = sources["specimens"].get(specimen, {}) if specimen else {}
    setup = fcfg.get("setup", {}) or {}
    if fcfg["campaign"] == 1 and boundary in setup:
        s1 = setup[boundary]
        sensor_loc = [s1.get("accel"), None, None]
        impact_loc = s1.get("impact")
    else:
        sensor_loc = list(setup.get("sensors", [])) + [None] * 3
        impact_loc = setup.get("impact")
    chans = []
    for k, (name, unit) in enumerate(zip(head["names"], head["units"], strict=True)):
        role = fcfg["channels"][k]
        sensor = hw.get(role) if role in ("accelerometer", "hammer") else None
        chans.append(
            Channel(
                index=k,
                name=name,
                unit=unit,
                quantity_header=unit_info(unit)[2],
                role=role,
                sensor_model=sensor["model"] if sensor else None,
                sensitivity_mv=sensor["sensitivity_mv"] if sensor else None,
                full_scale=(hw["daq"]["range_v"] / (sensor["sensitivity_mv"] * 1e-3))
                if sensor
                else None,
                location=(
                    sensor_loc[k] if role == "accelerometer" and sensor_loc[k] != "TODO" else None
                ),
            )
        )
    return Entry(
        id=f"{config_id}_{rep:02d}" if config_id.startswith("LL_") else config_id,
        path=str(path.relative_to(raw_dir)),
        sha256=sha256(path),
        campaign=fcfg["campaign"],
        config_id=config_id,
        specimen=specimen,
        material=spec.get("material"),
        boundary=boundary,
        hammer_position=hpos,
        accel_position=apos,
        repetition=rep,
        timestamp=rec.timestamp,
        fs_hz=rec.fs_hz,
        n_samples=rec.n_samples,
        duration_s=rec.duration_s,
        channels=chans,
        anomalies=list(rec.anomalies),
        impact_location=impact_loc,
    )


def build_catalog(raw_dir, sources_yaml) -> Catalog:
    raw_dir = Path(raw_dir)
    sources = yaml.safe_load(Path(sources_yaml).read_text(encoding="utf-8"))
    files = sorted(p for p in raw_dir.rglob("*.csv"))
    entries = [build_entry(p, raw_dir, sources) for p in files]
    entries.sort(key=lambda e: (e.campaign, e.config_id, e.repetition or 0))
    return Catalog(
        generated_by="dampid.io.catalog.build_catalog",
        sources=str(sources_yaml),
        n_files=len(entries),
        entries=entries,
    )


def save_catalog(cat: Catalog, path) -> None:
    Path(path).write_text(
        yaml.safe_dump(cat.model_dump(), allow_unicode=True, sort_keys=False, width=100),
        encoding="utf-8",
    )


def load_catalog(path) -> Catalog:
    return Catalog.model_validate(yaml.safe_load(Path(path).read_text(encoding="utf-8")))
