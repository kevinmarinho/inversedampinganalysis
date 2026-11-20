# Diário de sessões

## 2026-10-31 — Fase 0 (setup e auditoria do legado)

**Decisões registradas:** D0.1 uv · D0.2 YAML + pydantic · D0.3 estrutura §8 + notebooks/relatórios por etapa + GitHub · D0.4 código em inglês, docs em pt-BR · D0.5 trilha A primeiro · D0.6 comparação de métodos inversos na Fase 8 (todas pela recomendação, aceitas pelo pesquisador).

**Informação nova do pesquisador:** sensores fixados por rosca numa base impressa em PLA, colada com cera no corpo de prova; posições nas figuras do TCC/relatório. Massa do PCB 352C33: 5,8 g (datasheet rev. H).

**Feito**
- Repositório, ambiente (spec), legado intacto com hashes, dados brutos com manifesto.
- 8 notebooks executados (`notebooks/fase0/00`–`07`) e 8 relatórios (`docs/reports/fase0/`).
- `docs/audit_legacy.md`: 20 achados confirmados, 2 parciais, 1 inconclusivo (C9), 4 novos (N1–N4).
- Pacote: `dampid.analytic.beam` (soluções fechadas) e `dampid.audit` (harness do legado + réguas de referência), 17 testes passando (+1 pulado sem SymPy).

**Resultados principais**
- C2 é mais grave que o previsto: o ζ do código é amplificado ~10⁵–10⁶ ⇒ os "ζ ≈ 10⁻⁶" da monografia são ζ físicos de 4 %–200 %.
- C7 quantificado: 1 % de erro em frequência ⇒ ζ ajustado ×3,7.
- N3: alumínio livre-livre (−17 %) explicável por ≈7,5 g por conjunto sensor nas pontas; aço (+3 %) não — dimensões precisam ser medidas.
- N4: modelo de 6 nós tem metade dos modos acima de Nyquist ⇒ aliasing na simulação temporal.

**Limitações do ambiente:** PyPI bloqueado pela política de rede (proxy 403) — sem `uv.lock`, SymPy, Jupyter, `geneticalgorithm2`. Contornos documentados em `docs/reports/fase0/00_setup.md`.

**Pendências do pesquisador:** ver `docs/audit_legacy.md` §5, mais: criar o repositório GitHub (D0.3), copiar o `prompt_master.md` para `CLAUDE.md`, rodar `uv lock && uv sync && uv run pytest` localmente.

---

## 2026-11-20 — Fase 1 (inventário e EDA)

**Gate da Fase 0:** tratado como aprovado ("continue fazendo o processo", 2026-10-31). Repositório GitHub `kevinmarinho/inversedampinganalysis` criado pelo pesquisador; Fase 0 na `main`, Fase 1 no branch `fase-1`.

**Decisões registradas:** D1.1 conjunto de desenvolvimento = aço + Al livre-livre · D1.2 critério de classificação por canal e gravidade (`configs/qc_rules.yaml`). Ambas pela recomendação, por delegação do pesquisador.

**Feito**
- 96 CSVs que faltavam baixados do Projeto → `data/raw/` com 156 arquivos (manifesto ampliado).
- `dampid.io` (leitor NI CSV robusto, unidades, catálogo pydantic) e `dampid.eda` (qualidade, levantamento, consistência, classificação); 27 testes.
- 6 notebooks (`notebooks/fase1/01`–`06`), 6 relatórios + síntese (`docs/reports/fase1/`), `data/catalog.yaml` com classe e qualidade por arquivo.

**Resultados principais**
- Unidades g/N confirmadas por 3 caminhos físicos (ruído entre sessões, velocidade do martelo, linha de massa κ ≈ 1,5–2,4).
- C9: hipótese "mesmo sinal nas duas entradas" rejeitada; escala de força da 1ª campanha segue desconhecida.
- Saturação dos acelerômetros: fibra de carbono 16/16 (perdida para FRF), acrílico 5/16.
- Banda útil da força nos metais: f₋₁₀ ≈ 150–190 Hz → só o 1º modo elástico bem excitado.
- Aço truncado (decai 25 dB em 5 s; ζ ~ 7·10⁻⁴ grosseiro) e vibrando antes do golpe → D2.2.
- Repetibilidade: aço f₁ = 157,62 ± 0,05 Hz, Al 152,01 ± 0,35 Hz, γ² ≥ 0,99; FV 0° não repetível.
- Classificação: 52 ouro, 69 utilizável, 23 descartar, 12 pendente; dev set: aço 14 ouro, Al 13 ouro.

**Limitações do ambiente:** PyPI continua bloqueado (sem xarray/NetCDF: tabelas processadas em CSV, regeneráveis).


**Adendo (mesmo dia, respostas do pesquisador):** posições/orientação dos sensores (fotos), impacto livre-livre sem ponto fixo, 1ª campanha com impacto em L/4, dimensões da Tabela 1 (D1.3), ponteiras metal/borracha (inferidas pelo contato), FV alinhado com golpe variável e tubo girando. Notebook 07: impacto efetivo perto do meio das barras; linha de massa 1,2–1,4× acima do modelo nominal; **C9: Input 1 e Input 2 são a mesma força** (|In1/In2| = 1,13–1,20, plana) → D1.4 proposta.

---

## Estado atual (§14 do prompt master)

- **Fase atual:** 1 — concluída, **aguardando aprovação do gate** (`docs/reports/fase1/00_sintese_fase1.md`). Branch `fase-1`.
- **Última decisão registrada:** D1.3 (2026-11-20); D1.4 (canal de força da 1ª campanha) aguardando o pesquisador.
- **Dados:** 156 CSVs em `data/raw/`, catalogados e classificados em `data/catalog.yaml`.
- **Pendências do pesquisador** (detalhes na síntese da Fase 1):
  - [ ] Massas: espécimes e um conjunto sensor + base PLA + cera + cabo
  - [ ] D1.4: canal de força da 1ª campanha
  - [ ] 1ª campanha: posição do impacto por arquivo e comprimento livre no engaste
  - [ ] Fibra de carbono: dimensões/layup; repetir o ensaio sem saturar?
  - [ ] FV: espessura/layup; PLA: cotas e parâmetros de impressão
  - [ ] Copiar o `prompt_master.md` para `CLAUDE.md`; rodar `uv lock && uv sync && uv run pytest` localmente
- **Próximo passo:** aprovação do gate da Fase 1 → merge de `fase-1` na `main` → plano da Fase 2 (pipeline de sinais) com D2.1–D2.3.
