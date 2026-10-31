# Relatório 0.4 — Matriz de amortecimento (C2, C4)

**Base de dados:** nenhum dado experimental; matrizes K e M do próprio legado (6 nós, engaste-livre) e caso sintético com ζ prescrito.
**Notebook:** `notebooks/fase0/03_auditoria_amortecimento.ipynb`

## Feito
Medi o ζ que o modelo legado **realmente** tem, por projeção modal (ζ_r = φ_rᵀCφ_r / 2ω_r) e contando modos superamortecidos no espaço de estados. As duas medidas estão testadas contra a construção correta `C = MΦ diag(2ζω) ΦᵀM`.

## Evidência
| ζ de entrada do legado | ζ real, Al (modos 1/2/4) | ζ real, aço | ζ real, acrílico |
|---|---|---|---|
| 10⁻⁶ | 0,37 / 0,82 / 0,84 | 0,04 / 0,09 / 0,09 | 0,93 / 2,0 / 2,1 |
| 10⁻⁴ | 37 / 82 / 84 | 4 / 9 / 9 | 93 / 200 / 210 |
| 10⁻² | ~4000–8000 | ~400–900 | ~10⁴ |

A construção correta devolve exatamente ζ = 1 % em todos os modos.

| Item | Status |
|---|---|
| **C2** | **confirmado, mais grave do que o previsto**: o ζ de entrada é amplificado por ~10⁵–10⁶, com fator diferente por modo e material (o 3º modo, flexão no plano, quase não é afetado). |
| **C4** | **parcial, reformulado**: os modos 7–30 não ficam com ζ = 0; recebem amortecimento **descontrolado** (≈0 a 9) por vazamento, porque `Φ diag(ζ) Φᵀ` não é diagonal na base modal. |

## O que isso explica
- Os "ζ ótimos ~10⁻⁶" e a frase "1 % ζ = 10⁻⁶" da monografia: 10⁻⁶ é a faixa em que o modelo legado ainda oscila.
- "O modelo superestima o amortecimento": ζ_legado = 10⁻⁶ já é um amortecimento físico de 4 % a 200 %, muito acima do esperado para metais (10⁻⁴–10⁻²). Isso sustenta a hipótese **H1** do projeto.

## Próximo passo
Etapa 0.5 — graus de liberdade e condições de contorno (C3, I12).
