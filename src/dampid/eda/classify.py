"""Record classification (ouro / utilizável / descartar / pendente) — decision D1.2.

Inputs are the Phase-1 tables (survey, per-record consistency, per-config coherence);
thresholds come from ``configs/qc_rules.yaml``. Pure function of those tables, so the
classification is regenerable and auditable: every record carries the list of reasons.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def classify(
    survey: pd.DataFrame,
    consistency: pd.DataFrame,
    config_coherence: dict[tuple[str, int], float],
    rules: dict,
    noise_floor_g: float,
) -> pd.DataFrame:
    """Return one row per record: id, file_class, channel_status (per accelerometer),
    flags (non-blocking) and reasons."""
    cons = consistency.set_index(["id", "canal"])
    med_w = survey.groupby("config").F_width_ms.transform("median")
    out = []
    for i, r in survey.iterrows():
        flags, reasons, ch_status = [], [], {}
        multi = r.second_pulse_ratio >= rules["multi_impact_ratio"] or r.pulses_before_main > 0
        atyp = abs(r.F_width_ms / med_w[i] - 1) > rules["atypical_width_rel"]
        if multi:
            flags.append(f"impacto múltiplo (2º pulso {r.second_pulse_ratio:.0%} do principal)")
        if atyp:
            flags.append(
                f"golpe atípico (largura {r.F_width_ms:.2f} ms vs. mediana {med_w[i]:.2f})"
            )
        for ch in (0, 1):
            tag = f"a{ch}"
            if f"{tag}_n_clipped" not in r or pd.isna(r[f"{tag}_n_clipped"]):
                continue  # channel is not an accelerometer in this campaign
            problems = []
            if r[f"{tag}_n_clipped"] > 0:
                ch_status[ch] = "saturado"
                continue
            key = (r.id, ch)
            if key in cons.index and bool(cons.loc[key, "destoa"]):
                problems.append(
                    f"Acel {ch + 1} destoa da média ({cons.loc[key, 'desvio_dB']:.1f} dB)"
                )
            coh = config_coherence.get((r.config, ch), np.nan)
            if coh < rules["config_coherence_min"]:
                problems.append(
                    f"Acel {ch + 1}: configuração pouco repetível (γ² mediana {coh:.2f})"
                )
            if r[f"{tag}_decay_db"] > rules["truncation_decay_db"]:
                flags.append(f"Acel {ch + 1} truncado (decai só {-r[f'{tag}_decay_db']:.0f} dB)")
            ch_status[ch] = "ok" if not problems else "restrito"
            flags += problems
        pre_db = 20 * np.log10(r.a0_pre_rms / noise_floor_g)
        if pre_db > rules["pre_vibration_db"]:
            flags.append(f"já vibrava antes do golpe ({pre_db:.0f} dB acima do piso)")
        if r.campaign == 1:
            flags.append(rules["campaign1_note"])
        usable = [c for c, s in ch_status.items() if s != "saturado"]
        if not usable:
            cls = "descartar"
            reasons.append(
                "os dois acelerômetros saturaram (±50 g)"
                if len(ch_status) > 1
                else "o acelerômetro saturou (±50 g)"
            )
        elif (
            all(s == "ok" for s in ch_status.values())
            and not multi
            and not atyp
            and r.campaign != 1
        ):
            cls = "ouro"
        else:
            cls = "utilizável"
            for c, s in ch_status.items():
                if s == "saturado" and len(ch_status) > 1:
                    reasons.append(f"Acel {c + 1} saturado — usar só Acel {2 - c}")
            if multi:
                reasons.append("impacto múltiplo")
            if atyp:
                reasons.append("golpe atípico")
            reasons += [f for f in flags if "destoa" in f or "repetível" in f]
            if r.campaign == 1:
                reasons.append(rules["campaign1_note"])
        out.append(
            {
                "id": r.id,
                "file_class": cls,
                "reasons": reasons,
                "flags": sorted(set(flags)),
                "channel_status": ch_status,
            }
        )
    return pd.DataFrame(out)


def infer_tip(width_ms: float, rules: dict) -> str:
    """Hammer tip inferred from the contact time (time above 10 % of the peak).

    Hertzian contact: a softer tip (rubber) gives a longer contact than a metal tip on the
    same structure. The local stiffness of the specimen also matters, so contact times in
    between the two thresholds are left undetermined.
    """
    if width_ms <= rules["tip_metal_max_ms"]:
        return "metal (inferido do contato)"
    if width_ms >= rules["tip_rubber_min_ms"]:
        return "borracha (inferido do contato)"
    return "indeterminado"
