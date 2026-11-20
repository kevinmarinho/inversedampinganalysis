# %% [markdown]
# # Fase 1 · Etapa 1.5 — Classificação dos ensaios (D1.2) e catálogo final
#
# > **Bases de dados usadas:** **dados do pesquisador** (156 CSVs), resumidos nas tabelas das
# > Etapas 1.3–1.4 (`data/processed/fase1_*.csv`, regeneráveis). Critério: D1.2,
# > limiares em `configs/qc_rules.yaml`. Nenhuma referência externa.
#
# Regras (D1.2), em uma frase cada:
# * **descartar** — os dois acelerômetros saturaram: a FRF do registro é inválida;
# * **utilizável** — pelo menos um canal sem saturação, mas com impacto múltiplo, golpe
#   atípico, repetição que destoa, ou configuração pouco repetível;
# * **ouro** — impacto único e típico, nenhum canal saturado, repetição coerente com as
#   demais; truncamento e vibração pré-impacto são **sinalizados**, não rebaixam (são
#   propriedade do espécime/protocolo, tratados na Fase 2 — D2.2);
# * **1ª campanha** — no máximo **utilizável** (D1.4): um registro por configuração, sem
#   coerência entre repetições, e escala da força ±15 %. O canal de força é escolhido aqui:
#   Input 2 [N]; onde o Input 2 está sem sinal, Input 1 dividido pela razão In1/In2 mediana.

# %%
import numpy as np
import pandas as pd
import yaml

from dampid.eda.c9 import assign_force_channels, duplicate_hammer_table
from dampid.eda.classify import classify, infer_tip
from dampid.io.catalog import load_catalog, save_catalog

rules = yaml.safe_load(open("configs/qc_rules.yaml", encoding="utf-8"))
sv = pd.read_csv("data/processed/fase1_survey.csv")
cons = pd.read_csv("data/processed/fase1_consistency.csv")
coh_tab = pd.read_csv("docs/reports/fase1/05_consistencia_por_config.csv")
coh = {(r.config, int(r.canal[-1]) - 1): float(r["γ² mediana na banda"]) for _, r in coh_tab.iterrows()}
floor = float(sv[(sv.campaign == 2) & (sv.config != "LL_ACO")].a0_min_strms.median())
cls = classify(sv, cons, coh, rules, noise_floor_g=floor)
print(f"piso de ruído usado: {floor * 1e3:.3f} mg rms")

# %% [markdown]
# ## 1. Contagem por configuração

# %%
m = cls.merge(sv[["id", "config"]], on="id")
counts = pd.crosstab(m.config, m.file_class).reindex(columns=["ouro", "utilizável", "descartar"],
                                                      fill_value=0)
counts["total"] = counts.sum(axis=1)
camp1 = counts.index.str.startswith("C1_")
counts = pd.concat([counts[~camp1], counts[camp1].sum().to_frame("1ª campanha (12 configs)").T])
counts.loc["TOTAL"] = counts.sum()
print(counts.to_string())
counts.to_csv("docs/reports/fase1/06_classificacao_por_config.csv")

# %% [markdown]
# ## 2. Motivos (registros não-ouro)

# %%
nz = m[m.file_class != "ouro"]
with pd.option_context("display.width", 250, "display.max_colwidth", 140):
    print(nz[["id", "file_class", "reasons"]].assign(reasons=nz.reasons.map("; ".join)).to_string(index=False))

# %% [markdown]
# ## 3. Sinalizações não bloqueantes (por configuração)

# %%
fl = m.explode("flags").dropna(subset=["flags"]).reset_index(drop=True)
fl["tipo"] = fl["flags"].str.replace(r"\(.*\)", "", regex=True).str.strip()
print(pd.crosstab(fl.config, fl.tipo).to_string())

# %% [markdown]
# ## 4. Conjunto de desenvolvimento (D1.1: aço e alumínio livre-livre)

# %%
dev = m[m.config.isin(["LL_ACO", "LL_AL"])]
print(pd.crosstab(dev.config, dev.file_class).to_string())
for cfg in ["LL_ACO", "LL_AL"]:
    ids = list(dev[(dev.config == cfg) & (dev.file_class == "ouro")].id)
    print(f"{cfg} ouro ({len(ids)}):", ", ".join(i.split("_")[-1] for i in ids))

# %% [markdown]
# ## 5. Catálogo final: classe, motivos e resumo de qualidade em cada entrada

# %%
cat = load_catalog("data/catalog.yaml")
c9tab = duplicate_hammer_table(cat)
r_med = assign_force_channels(cat, c9tab)
print(c9tab.round(3).to_string(index=False))
print(f"razão In1/In2 mediana usada onde o Input 2 está sem sinal: {r_med:.3f}")
svi = sv.set_index("id")
clsi = cls.set_index("id")
keep = ["n_pulses", "second_pulse_ratio", "F_peak_N", "F_width_ms", "F_impulse_Ns", "F_f10db_hz",
        "F_f20db_hz", "a0_n_clipped", "a1_n_clipped", "a0_decay_db", "a1_decay_db", "t_impact_s"]
for e in cat.entries:
    r = svi.loc[e.id]
    e.quality = {k: (None if pd.isna(r[k]) else round(float(r[k]), 4)) for k in keep}
    e.quality["channel_status"] = {f"Acel {k + 1}": v for k, v in clsi.loc[e.id, "channel_status"].items()}
    e.quality["flags"] = list(clsi.loc[e.id, "flags"])
    e.hammer_tip = infer_tip(float(r["F_width_ms"]), rules)
    e.file_class = clsi.loc[e.id, "file_class"]
    e.class_reasons = list(clsi.loc[e.id, "reasons"])
save_catalog(cat, "data/catalog.yaml")
tips = pd.crosstab(pd.Series([e.config_id if e.campaign == 2 else "1ª campanha" for e in cat.entries], name="config"),
                   pd.Series([e.hammer_tip for e in cat.entries], name="ponteira"))
print(tips.to_string())
cat2 = load_catalog("data/catalog.yaml")
print(pd.DataFrame([{"id": e.id, "canal de força": e.force_channel,
                    "fator p/ N": round(e.force_scale_to_N, 4), "nota": e.force_note}
                   for e in cat.entries if e.campaign == 1]).to_string(index=False))
print("catálogo salvo e revalidado:", cat2.n_files, "entradas;",
      dict(pd.Series([e.file_class for e in cat2.entries]).value_counts()))
