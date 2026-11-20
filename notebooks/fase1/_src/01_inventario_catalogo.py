# %% [markdown]
# # Fase 1 · Etapa 1.1 — Inventário e catálogo (`data/catalog.yaml`)
#
# > **Bases de dados usadas:**
# > * **Dados do pesquisador:** os 156 CSVs de `data/raw/` (1ª campanha `exp_5/`, 12
# >   arquivos; 2ª campanha livre-livre, 9 configurações × 16 repetições), com integridade
# >   conferida por `data/raw/MANIFEST.sha256`; arquivos `.task` do NI (sensibilidades).
# > * **Referências externas (só metadados de hardware):** especificações do NI 9234
# >   (faixa ±5 V, corte AC 0,5 Hz — ni.com), martelo PCB 086C03 (2,25 mV/N, ±2224 N —
# >   pcb.com) e acelerômetro PCB 352C33 (±50 g, 5,8 g — datasheet rev. H).
# > Nenhum dado de terceiros entra nas medições.
#
# O que este notebook faz: lê o cabeçalho e o conteúdo de cada arquivo, confere estrutura
# (canais, unidades, fs, duração, uniformidade do tempo), e gera o catálogo validado por
# pydantic a partir de `configs/data_sources.yaml`.

# %%
from collections import Counter

import pandas as pd

from dampid.io.catalog import build_catalog, load_catalog, save_catalog

cat = build_catalog("data/raw", "configs/data_sources.yaml")
save_catalog(cat, "data/catalog.yaml")
cat = load_catalog("data/catalog.yaml")  # round-trip: o arquivo salvo é válido
print(f"{cat.n_files} arquivos catalogados -> data/catalog.yaml")

# %% [markdown]
# ## 1. Resumo por configuração de ensaio

# %%
rows = []
for cfg, es in cat.by_config().items():
    rows.append({
        "config": cfg,
        "camp.": es[0].campaign,
        "CC": es[0].boundary,
        "n": len(es),
        "datas": ", ".join(sorted({e.timestamp[:10] for e in es})),
        "fs [Hz]": sorted({e.fs_hz for e in es}),
        "duração [s]": sorted({round(e.duration_s, 3) for e in es}),
        "canais": " | ".join(f"{c.name} [{c.unit}]" for c in es[0].channels),
        "anomalias": sum(len(e.anomalies) > 0 for e in es),
    })
inv = pd.DataFrame(rows)
with pd.option_context("display.width", 250, "display.max_colwidth", 70):
    print(inv.to_string(index=False))
inv.to_csv("docs/reports/fase1/01_inventario_configuracoes.csv", index=False)

# %% [markdown]
# ## 2. Conferências estruturais
# * todos os arquivos têm 3 canais, com a mesma unidade por posição (g, g, N)?
# * fs é igual em todos? as colunas de tempo avançam exatamente `dt`?
# * nomes de canal mudam entre arquivos de uma mesma configuração?

# %%
units = Counter(tuple(c.unit for c in e.channels) for e in cat.entries)
fsset = Counter(round(e.fs_hz, 6) for e in cat.entries)
names = Counter((e.campaign, tuple(c.name for c in e.channels)) for e in cat.entries)
print("unidades por posição de canal:", dict(units))
print("taxas de amostragem:", dict(fsset))
print("nomes de canais (campanha, nomes): n arquivos")
for k, v in names.items():
    print("  ", k, v)
dur = Counter((e.config_id.split("_m")[0] if e.campaign == 1 else e.config_id, round(e.duration_s, 3))
              for e in cat.entries)
print("\nduração por configuração:")
for k, v in sorted(dur.items()):
    print("  ", k, v)

# %%
print("Arquivos com anomalias de estrutura:")
for e in cat.entries:
    for a in e.anomalies:
        print(f"  {e.path}: {a}")

# %% [markdown]
# ## 3. Tabela completa por arquivo (exportada para o relatório)

# %%
full = pd.DataFrame([{
    "id": e.id, "arquivo": e.path, "campanha": e.campaign, "config": e.config_id,
    "espécime": e.specimen, "CC": e.boundary, "rep.": e.repetition, "data/hora": e.timestamp,
    "fs [Hz]": e.fs_hz, "amostras": e.n_samples, "duração [s]": round(e.duration_s, 4),
    "martelo (nome do arquivo)": e.hammer_position, "acel. (nome do arquivo)": e.accel_position,
    "anomalias": "; ".join(e.anomalies),
} for e in cat.entries])
full.to_csv("docs/reports/fase1/01_catalogo_arquivos.csv", index=False)
print(full.groupby(["campanha", "config"]).size().to_string())

# %% [markdown]
# ## 4. Campos que continuam `TODO` (só o pesquisador tem)
# Posição e direção de cada acelerômetro, ponto de impacto, ponteira do martelo, cotas de
# fibra de carbono, espessura do tubo de fibra de vidro, cotas e parâmetros de impressão
# do PLA, ligas do aço e do alumínio, e o que diferencia `PLA_LONGARINA` de
# `PLA_LONGARINA_90`. Ficam como `null` no catálogo — nunca são preenchidos por palpite.

# %%
todo = {
    "location (posição/direção do sensor)": sum(c.location is None for e in cat.entries for c in e.channels if c.role == "accelerometer"),
    "impact_location": sum(e.impact_location is None for e in cat.entries),
    "hammer_tip": sum(e.hammer_tip is None for e in cat.entries),
}
print(todo)
