# Relatório 0.3 — Métricas e função objetivo (C1, I15, N1)

**Base de dados:** dados do pesquisador — `exp_5/Acrilico_engastado_livre_mponta_aponta.csv` (caso do `main_otimizador.py`).
**Notebook:** `notebooks/fase0/02_auditoria_metricas.ipynb`

## Feito
Avaliei as funções objetivo do legado (MAC e RVAC) em função de ζ, usando o próprio código legado.

## Evidência
| Item | Status | Evidência |
|---|---|---|
| **C1** | **confirmado** | O MAC do legado é idêntico até a 15ª casa decimal (0,980924058977776) para ζ de 10⁻⁷ a 0,5: compara **vetores de frequências** (Hz × rad/s), e os autovalores vêm do problema sem amortecimento. Gradiente nulo ⇒ L-BFGS-B "para no chute". O número não tem significado físico: duas listas crescentes sem relação dão MAC = 0,83. |
| **N1 (novo)** | **confirmado** | `frf_direta_modal_comparacao` cria vetores com tamanho fixo 4000, mas calcula só em 1599 frequências; `rvac` quebra com erro de dimensão. A métrica do `main_otimizador.py` não roda com estes dados. |
| **I15** | **confirmado** | Com o tamanho recortado (única mudança), a RVAC varia com ζ e tem máximo em ζ_legado ≈ 1e-4 — **fora** dos limites do GA (0 a 1e-5). O GA terminaria encostado no limite. Melhor RVAC possível ≈ 0,56 (ajuste ruim, porque as frequências não batem). Busca: 8 indivíduos × 20 gerações em 6 dimensões. |

RVAC e TSAC são invariantes à escala: multiplicar a FRF numérica por 1000 não muda a métrica — por isso os fatores 9,81/9810/9810² (C5) nunca apareceram como erro.

## Problemas / surpresas
- Eu esperava a RVAC plana dentro dos limites do GA; **não é** (varia de 0,13 a 0,52). O problema real é o ótimo estar fora da faixa e o parâmetro não ser o ζ físico (C2).
- O GA não pôde ser executado (pacote não instalável aqui, e a métrica quebra antes).

## Próximo passo
Etapa 0.4 — matriz de amortecimento (C2, C4).
