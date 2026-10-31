# Relatório 0.2 — Mapa do legado e execução "como está"

**Base de dados:** dados do pesquisador — `exp_5/` (Al e acrílico, engaste-livre) e `AÇO/` (livre-livre).
**Notebook:** `notebooks/fase0/01_tour_legado.ipynb`

## Feito
- Inventário estático das 20 funções (tamanho, quem chama quem) e fluxo de dados de `main.py` e `main_otimizador.py`.
- Execução dos dois scripts sem alteração, via harness.

## Evidência
| Item | Resultado |
|---|---|
| `main.py` | roda em ≈1 s e gera as duas figuras (tempo e FRF) |
| **I1** | **confirmado** — `main_otimizador.py` quebra com `ValueError: too many values to unpack (expected 17)`. A partir da 10ª posição as variáveis ficariam deslocadas (`fase_frf` ← `signal_filtered` …). Os resultados de otimização da monografia vieram de uma versão anterior que não está no repositório. |
| **I2** | **confirmado** — `contorno("livre")` dá `UnboundLocalError`; nenhuma função numérica aplica força em livre-livre. O modelo legado não simula a 2ª campanha. |
| **N2 (novo)** | `pos_experimental_livre_livre(..., max_time=4 s)` quebra no aço: se o trigger cai depois de ~1 s, a FFT daquele ensaio tem 1987 bins em vez de 2000 e a média falha. Funciona com 3,5 s. |
| **I11 (reforço)** | no aço livre-livre, o trigger do legado (diferença de aceleração > 0,05 m/s²) dispara **antes do golpe** em 9 de 16 ensaios, em 4 deles na amostra 0 (até 2 s antes da força). |

Funções nunca chamadas: `calculo_numerico` (Rayleigh), `frf_tempo`, `frf_direta`, `otimizador`, `pos_experimental_livre_livre`.

## Problemas / surpresas
- O trigger antecipado no aço vem da viga ainda vibrando antes do golpe. Qualquer pipeline novo precisa disparar pela **força** (já previsto na Fase 2).

## Próximo passo
Etapa 0.3 — métricas e função objetivo (C1, I15).
