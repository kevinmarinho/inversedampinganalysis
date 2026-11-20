"""Impact-test processing pipeline: trigger -> windows -> spectra -> H1/H2/coherence.

One configuration object drives everything, so measured and simulated records go through
exactly the same steps (prompt master §6.2).
"""

from __future__ import annotations

import numpy as np
from pydantic import BaseModel, Field

from dampid.eda.quality import find_pulses, main_pulse
from dampid.signal.frf import FRFEstimate, estimate_frf
from dampid.signal.windows import exponential_window, force_window


class ProcessingConfig(BaseModel):
    pre_trigger_s: float = Field(0.01, description="segment starts this long before the impact")
    length_s: float | None = Field(None, description="segment length; None = longest common")
    force_margin_s: float = Field(1e-3, description="force window extends this past the pulses")
    force_taper_s: float = Field(0.5e-3, description="half-cosine taper of the force window")
    exp_tau_s: float | None = Field(None, description="exponential window tau (D2.2); None = off")
    pulse_rel_threshold: float = 0.05
    multi_impact: str = Field("include", description="'include': force window spans all pulses")


def trigger(force: np.ndarray, fs: float, cfg: ProcessingConfig):
    """Main impact (largest pulse) found on the FORCE channel, plus all pulses."""
    pulses = find_pulses(force, fs, rel_threshold=cfg.pulse_rel_threshold)
    if not pulses:
        raise ValueError("no force pulse found")
    return pulses[main_pulse(pulses)], pulses


def prepare(force, resp, fs: float, cfg: ProcessingConfig, n_seg: int | None = None):
    """Cut, window and transform one record. resp: (n_ch, n) or (n,).

    Returns (F, X, info) with F (n_f,), X (n_ch, n_f) complex spectra (rfft, no scaling:
    the FRF is a ratio), and info with indices, pulses and windows used.
    """
    force = np.nan_to_num(np.asarray(force, dtype=float))
    resp = np.atleast_2d(np.nan_to_num(np.asarray(resp, dtype=float)))
    main, pulses = trigger(force, fs, cfg)
    i0 = main.i_start - int(round(cfg.pre_trigger_s * fs))
    if i0 < 0:
        raise ValueError("impact too close to the start of the record for the pre-trigger")
    avail = force.size - i0
    n = avail if n_seg is None else n_seg
    if n > avail:
        raise ValueError(f"segment of {n} samples does not fit ({avail} available)")
    f_seg = force[i0 : i0 + n].copy()
    x_seg = resp[:, i0 : i0 + n].copy()
    # baseline of the force from the pre-trigger (AC-coupled channel: ~0, but remove it)
    n_pre = max(1, int(round(cfg.pre_trigger_s * fs)) - int(round(cfg.force_margin_s * fs)))
    f_seg -= np.median(f_seg[:n_pre])
    in_seg = [p for p in pulses if 0 <= p.i_start - i0 < n]
    last = max(in_seg, key=lambda p: p.i_end) if cfg.multi_impact == "include" else main
    margin = int(round(cfg.force_margin_s * fs))
    wf = force_window(
        n, main.i_start - i0 - margin, last.i_end - i0 + margin, int(round(cfg.force_taper_s * fs))
    )
    f_seg *= wf
    i_imp = main.i_start - i0
    if cfg.exp_tau_s is not None:
        # window origin = start of the force window: the force is exactly zero before it,
        # which makes the pole shift of -1/tau exact (windows.exponential_window)
        i_exp0 = max(0, main.i_start - i0 - margin - int(round(cfg.force_taper_s * fs)))
        we = exponential_window(n, fs, cfg.exp_tau_s, i_exp0)
        f_seg *= we
        x_seg *= we[None, :]
    info = {
        "i_segment_start": int(i0),
        "n_segment": int(n),
        "i_impact_in_segment": int(i_imp),
        "n_pulses_in_segment": len(in_seg),
        "exp_tau_s": cfg.exp_tau_s,
        "i_exp_origin_in_segment": int(i_exp0) if cfg.exp_tau_s is not None else None,
    }
    return np.fft.rfft(f_seg), np.fft.rfft(x_seg, axis=1), info


def process(records, fs: float, cfg: ProcessingConfig, meta=None) -> FRFEstimate:
    """records: iterable of (force, resp) arrays. Segment length = cfg.length_s or the
    longest length common to all records."""
    records = list(records)
    if cfg.length_s is not None:
        n = int(round(cfg.length_s * fs))
    else:
        n = None
        for force, _ in records:
            main, _ = trigger(np.nan_to_num(np.asarray(force, dtype=float)), fs, cfg)
            avail = len(force) - (main.i_start - int(round(cfg.pre_trigger_s * fs)))
            n = avail if n is None else min(n, avail)
    Fs, Xs, infos = [], [], []
    for force, resp in records:
        F, X, info = prepare(force, resp, fs, cfg, n)
        Fs.append(F)
        Xs.append(X)
        infos.append(info)
    f = np.fft.rfftfreq(n, 1 / fs)
    m = {"config": cfg.model_dump(), "fs_hz": fs, "n_segment": n, "records": infos}
    m.update(meta or {})
    return estimate_frf(f, np.array(Fs), np.array(Xs), m)
