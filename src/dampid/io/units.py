"""Unit conversions at the data boundary (prompt master §6.1).

Every scale factor that turns a quantity read from a file into SI lives here, with a
name. Nothing else in the code base may multiply by 9.81, 9810, etc.
"""

from __future__ import annotations

import numpy as np

G0 = 9.80665  # standard gravity [m/s^2] (CGPM 1901 / ISO 80000-3)

# file unit -> (factor to SI, SI unit, physical quantity)
_UNIT_TABLE: dict[str, tuple[float, str, str]] = {
    "g": (G0, "m/s^2", "acceleration"),
    "m/s^2": (1.0, "m/s^2", "acceleration"),
    "m/s2": (1.0, "m/s^2", "acceleration"),
    "n": (1.0, "N", "force"),
    "lbf": (4.4482216152605, "N", "force"),
    "v": (1.0, "V", "voltage"),
}


def unit_info(unit: str) -> tuple[float, str, str]:
    """(factor_to_si, si_unit, quantity) for a unit string as written in the files."""
    key = unit.strip().strip('"').lower()
    if key not in _UNIT_TABLE:
        raise KeyError(f"unknown unit {unit!r}; add it to dampid.io.units._UNIT_TABLE")
    return _UNIT_TABLE[key]


def to_si(x, unit: str):
    """Convert values written in `unit` to SI. Returns (x_si, si_unit)."""
    factor, si_unit, _ = unit_info(unit)
    return np.asarray(x, dtype=float) * factor, si_unit


def full_scale(sensitivity_mv_per_unit: float, range_v: float = 5.0) -> float:
    """Largest measurable value (in the sensor's engineering unit) for an ADC of
    ±range_v volts and a sensor of the given sensitivity [mV per unit].

    NI 9234: ±5 V nominal (±5.1 V typical) AC full scale (NI-9234 Specifications).
    """
    return range_v / (sensitivity_mv_per_unit * 1e-3)
