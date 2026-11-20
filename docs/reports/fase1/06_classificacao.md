# Relatório 1.5 — Classificação dos ensaios (D1.2) e catálogo final

**Bases de dados:** dados do pesquisador (156 CSVs), via as tabelas das Etapas 1.3–1.4. Limiares em `configs/qc_rules.yaml`.
**Notebook:** `notebooks/fase1/06_classificacao.ipynb` · **Tabela:** `06_classificacao_por_config.csv`
**Código:** `dampid.eda.classify.classify` (testado em `tests/unit/test_classify.py`)

## Regras (D1.2, ver `docs/decisions/D1.2-criterio-rejeicao.md`)
- **descartar:** todos os acelerômetros do registro saturaram (dois na 2ª campanha, um na 1ª).
- **utilizável:** pelo menos um canal sem saturação, mas com impacto múltiplo, golpe atípico, repetição que destoa ou configuração pouco repetível. Canal saturado = canal excluído.
- **ouro:** impacto único e típico, nenhum canal saturado, coerente com as outras repetições. Truncamento e vibração pré-impacto só **sinalizam**.
- **1ª campanha:** no máximo "utilizável" (D1.4: um registro por configuração, escala da força ±15 %). Força = Input 2; no aço engaste-livre, Input 1 ÷ 1,135.

## Evidência
| Config | ouro | utilizável | descartar | pendente |
|---|---|---|---|---|
| **Aço LL** | **14** | 1 (rep. 10: golpe atípico + 2 pulsos) | 1 (rep. 11) | — |
| **Alumínio LL** | **13** | 2 (reps. 05, 07: impacto múltiplo) | 1 (rep. 13) | — |
| Acrílico LL | 0 | 11 | 5 | — |
| Fibra de carbono LL | 0 | 0 | **16** | — |
| FV 0° | 0 | 16 (não repetível) | 0 | — |
| FV 90° | 8 | 8 | 0 | — |
| PLA caverna | 4 | 12 | 0 | — |
| PLA longarina 0 | 13 | 3 | 0 | — |
| PLA longarina 90 | 0 | 16 (Acel 2 não repetível) | 0 | — |
| 1ª campanha | — | 10 | 2 (acelerômetro saturado: acrílico CF m-ponta, Al CF m-centro) | — |
| **Total** | **52** | **79** | **25** | **0** |

**Sinalizações não bloqueantes:** aço truncado em 15/15 registros não saturados e "já vibrava antes do golpe" em 10/16. Cada entrada do `data/catalog.yaml` traz agora `file_class`, `class_reasons`, `quality.channel_status` (por acelerômetro: ok / restrito / saturado) e `quality.flags`.

**Conjunto de desenvolvimento (D1.1):** aço com 14 registros ouro (todos menos 10 e 11) e alumínio com 13 (todos menos 05, 07 e 13).

## Problemas / surpresas
- A fibra de carbono sai inteira do uso por FRF. É o único compósito de carbono do banco, o que importa para a trilha B (materiais de satélite). Recomendo novo ensaio com golpe mais fraco/ponteira mais macia, se houver acesso ao laboratório.
- Posições, orientação dos sensores e ponto de impacto entraram no catálogo pelo adendo (relatório 07). A ponteira foi inferida pelo tempo de contato.

## Próximo passo
Síntese da Fase 1 e pedido de aprovação do gate.
