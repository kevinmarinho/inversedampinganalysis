import pandas as pd

from dampid.eda.classify import classify

RULES = {
    "multi_impact_ratio": 0.10,
    "atypical_width_rel": 0.5,
    "config_coherence_min": 0.9,
    "truncation_decay_db": -40.0,
    "pre_vibration_db": 25.0,
}


def _row(rid, **kw):
    base = dict(
        id=rid,
        config="C",
        campaign=2,
        second_pulse_ratio=0.0,
        pulses_before_main=0,
        F_width_ms=5.0,
        a0_n_clipped=0,
        a1_n_clipped=0,
        a0_decay_db=-60.0,
        a1_decay_db=-60.0,
        a0_pre_rms=1e-3,
    )
    base.update(kw)
    return base


def test_classes_follow_the_rules():
    sv = pd.DataFrame(
        [
            _row("good"),
            _row("double", second_pulse_ratio=0.3),
            _row("one_sat", a1_n_clipped=10),
            _row("both_sat", a0_n_clipped=5, a1_n_clipped=3),
            _row("trunc", a0_decay_db=-25.0, a1_decay_db=-25.0),
            _row("c1", campaign=1),
        ]
    )
    cons = pd.DataFrame(
        [{"id": i, "canal": c, "desvio_dB": 0.3, "destoa": False} for i in sv.id for c in (0, 1)]
    )
    out = classify(sv, cons, {("C", 0): 0.99, ("C", 1): 0.99}, RULES, noise_floor_g=5e-4).set_index(
        "id"
    )
    assert out.loc["good", "file_class"] == "ouro"
    assert out.loc["double", "file_class"] == "utilizável"
    assert out.loc["one_sat", "file_class"] == "utilizável"
    assert "usar só Acel 1" in out.loc["one_sat", "reasons"][0]
    assert out.loc["both_sat", "file_class"] == "descartar"
    # truncation is a flag for the Phase-2 decision, not a demotion
    assert out.loc["trunc", "file_class"] == "ouro"
    assert any("truncado" in f for f in out.loc["trunc", "flags"])
    assert out.loc["c1", "file_class"] == "pendente"


def test_unrepeatable_configuration_demotes_channel():
    sv = pd.DataFrame([_row("r1")])
    cons = pd.DataFrame(
        [{"id": "r1", "canal": c, "desvio_dB": 0.3, "destoa": False} for c in (0, 1)]
    )
    out = classify(sv, cons, {("C", 0): 0.4, ("C", 1): 0.95}, RULES, noise_floor_g=5e-4)
    assert out.loc[0, "file_class"] == "utilizável"
    assert out.loc[0, "channel_status"] == {0: "restrito", 1: "ok"}


def test_infer_tip_thresholds():
    from dampid.eda.classify import infer_tip

    r = {"tip_metal_max_ms": 3.5, "tip_rubber_min_ms": 5.0}
    assert infer_tip(2.1, r).startswith("metal")
    assert infer_tip(6.5, r).startswith("borracha")
    assert infer_tip(4.2, r) == "indeterminado"
