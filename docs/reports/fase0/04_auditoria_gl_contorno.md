# Relatório 0.5 — Graus de liberdade e condições de contorno (C3, I12)

**Base de dados:** nenhum dado experimental; modelo legado e solução analítica de Euler–Bernoulli (Blevins).
**Notebook:** `notebooks/fase0/04_auditoria_gl_contorno.ipynb`

## Evidência
| Item | Status | Evidência |
|---|---|---|
| **C3** | **confirmado** | Engaste-livre: força em `w` do nó da ponta (certo), resposta em `w'` (rotação). A razão rotação/translação na FRF é 5,08 m⁻¹ no 1º modo e 17,64 m⁻¹ no 2º — idêntica à inclinação/deslocamento analítica da ponta, ou seja, um fator diferente por modo que nenhuma constante corrige. |
| **I12** | **parcial** | Posição **confirmada**: no engaste-engaste, com 6 nós, força e resposta caem em x = 0,4 L (não há nó central). Malha: o elemento de Hermite converge rápido — erro de 0,001 % e 0,05 % nos dois modos do Al em balanço abaixo de 500 Hz; ~1 % só no 4º modo. Para metais em balanço a malha não é o problema; para acrílico e livre-livre faltou o estudo de convergência. |

## Próximo passo
Etapa 0.6 — elemento, propriedades e dimensões (I3, I4, I14, C8).
