# Relatório 1.6 (adendo) — Montagem informada, ponteiras, linha de massa e C9

**Bases de dados:** dados do pesquisador (respostas de 2026-11-20, fotos em `fotos/`, CSVs de aço/Al ouro e de `exp_5/`). Externas, só no modelo de comparação: ρ nominal (aço 7850, Al 2700 kg/m³) e tolerância de sensibilidade do martelo (±15 %, pcb.com).
**Notebook:** `notebooks/fase1/07_adendo_montagem.ipynb` · **Código novo:** `dampid.analytic.beam.free_free_with_point_masses` e `accelerance` (Rayleigh–Ritz testado contra as frequências analíticas, com convergência < 0,1 %); `dampid.eda.classify.infer_tip`

## Informações do pesquisador (registradas em `configs/data_sources.yaml` e no catálogo)
| Item | Resposta |
|---|---|
| Barras livre-livre | dois sensores perto das extremidades, face superior; impacto sem ponto específico |
| Longarina 0 / 90 | sensores com a mesma orientação / a 90° (um na face superior, outro na alma lateral) |
| FV alinhado | o ponto de golpe mudou entre repetições e o tubo girava → confirma a baixa repetibilidade (γ² ≈ 0,35) |
| 1ª campanha | biengastado: acelerômetro em L/2, impacto em L/4; engaste-livre: acelerômetro na ponta, impacto em L/4 do lado do acelerômetro |
| Dimensões | Tabela 1 do TCC (D1.3) |
| Ponteiras | metal e borracha, conforme o material; sem registro por ensaio |
| Ligas / certificados | não identificáveis / não existem; sensibilidades nominais valem |

## Evidência
**1. Ponteira inferida pelo contato** (≤ 3,5 ms metal, ≥ 5 ms borracha):

| Config | contato [ms] | ponteira | f₋₁₀ [Hz] |
|---|---|---|---|
| aço, Al, FV 0°/90°, PLA caverna | 5,6–6,5 | **borracha** | 150–190 |
| acrílico, fibra de carbono | 3,0–3,1 | **metal** | 250–630 |
| PLA longarinas | 2,1 | **metal** | 460–480 |

O catálogo traz `hammer_tip` com o sufixo "(inferido do contato)". Aço 10/11 (contato 2,2 e 0,4 ms) saem como "metal", mas são golpes atípicos, não troca de ponteira. Na fibra de carbono, ponteira metálica + ~40 N explica a saturação em 16/16.

**2. Linha de massa com a montagem real** (dimensões da Tabela 1, ρ nominal, sensores a 0,05 L e 0,95 L, EI ajustado para casar f₁; banda 40–120 Hz):

| Config | impacto | medido/modelo (sensor 0,05 L / 0,95 L), m_s = 6–12 g |
|---|---|---|
| Aço | **0,5 L** | **1,24–1,29 / 1,21–1,25** |
| Aço | 0,4 L | 0,96–1,01 / 1,98–2,02 |
| Al | **0,5 L** | **1,25–1,37 / 1,14–1,27** |
| Al | 0,4 L | 1,02–1,15 / 1,78–1,88 |

- Só o **impacto perto do meio** dá o mesmo nível nos dois sensores, como medido (H₂/H₁ ≈ 1). "Sem ponto específico" foi, na prática, sempre "por volta do meio". Isso combina com γ² ≥ 0,99 entre repetições e com o 2º modo (antissimétrico) quase ausente.
- Com impacto no meio, o nível medido fica **1,2–1,4× acima** do modelo nominal. A suspensão (pulo em ~15–20 Hz, fora do modelo) responde por +3 a +15 %. O resto (~10–25 %) cabe na tolerância das sensibilidades nominais ou numa massa real menor. **Só a pesagem decide.** Unidades seguem confirmadas: um erro de unidade daria um fator ≥ 9,8.

**3. C9 resolvido: Input 1 e Input 2 são a mesma força** (o pesquisador confirmou depois: o martelo estava ligado nas duas entradas). A função de transferência |In1/In2| é **plana** em 5–300 Hz: varia 0,4–4 dB e vale **1,125–1,198** nos 10 arquivos com Input 2 ativo. A exceção é acrílico engaste-livre m-ponta, com 7,8 dB e o Input 0 saturado. Já |In0/In2| varia 20–64 dB, com ressonâncias: Input 0 é a resposta da estrutura.
- Leitura: o sinal do martelo foi registrado duas vezes, com fatores de conversão diferentes em ~13–20 %. No aço engaste-livre o Input 2 estava sem sinal e só o Input 1 tem a força.
- **Correção da Etapa 1.2:** lá a hipótese "mesmo sinal nas duas entradas" foi rejeitada pela razão de picos no tempo. A rejeição foi forte demais, porque a deriva lenta de base do Input 1 distorce essa razão. No domínio da frequência a razão é constante.
- **Para ζ, os dois canais dão o mesmo resultado** (a forma da FRF é a mesma). A escala absoluta fica com ~15 % de incerteza. **D1.4 (aceita):** força = Input 2; no aço engaste-livre, Input 1 ÷ 1,135. A 1ª campanha passa a 10 utilizáveis + 2 descartados.
- Figura para o pesquisador: `figs/07_adendo_montagem_fig02.png` (12 arquivos: canais no tempo e |In1/In2|, |In0/In2|).

## Problemas / surpresas
- O nome dos arquivos da 1ª campanha (`mcentro`, `mponta`, `mforacentro`) não bate exatamente com "impacto em L/4" em todos os casos: é preciso confirmar a posição por arquivo.
- Continuam pendentes: massas, cotas e layup da fibra de carbono, espessura e layup do tubo de FV, e a possibilidade de repetir a fibra de carbono.
