# Relatório 0.7 — Sinais e unidades (C9, C5, C6, I7, I8, I10, I11)

**Bases de dados:** dados do pesquisador (`exp_5/` e livre-livre de aço/Al/acrílico) e casos sintéticos de 1 GL com ζ conhecido (gerados no notebook).
**Notebook:** `notebooks/fase0/06_auditoria_sinais.ipynb`

## Evidência
| Item | Status | Evidência |
|---|---|---|
| **C9** | **inconclusivo (reformulado)** | Nos 2 ensaios de **aço em engaste-livre** (inclui o da Tabela 5 do TCC), o Input 2 [N] só tem ruído e o pulso está no Input 1 — aí o legado provavelmente lê o canal certo, com rótulo "g" errado. Nos outros 10 arquivos, Input 1 e Input 2 têm pulsos quase idênticos. Sem o registro do ensaio não dá para saber o que estava em cada entrada nem a sensibilidade usada → escala de força da 1ª campanha desconhecida. |
| **C5** | **confirmado** | FRF experimental multiplicada por 2/N (1,6·10⁻⁵ para 128 000 amostras); fatores 9,81, 9810 e 9810² ajustados à mão. Escondido porque RVAC/TSAC são invariantes à escala. |
| **C6** | **confirmado** | 1 GL sintético a 150 Hz: com o filtro do código (Wn = 0,5) o ζ por meia-potência dobra para ζ = 0,1 % (+103 %) e fica quase intacto para ζ ≥ 0,5 %; com Wn = 0,02 (valor do TCC) o viés vai de +144 % a +3800 %. Viés maior justamente nos materiais de baixo amortecimento. |
| **I7** | **confirmado** | `frf_direta` põe o pico em 1405 Hz em vez de 35,6 Hz (fator 4π² = 39,5). |
| **I8** | **confirmado** | fase multiplicada por ½ (−1 → 90°). |
| **I10** | **confirmado** | média de fases produz curvas sem sentido (figura); média de magnitudes desalinhadas muda a largura do pico. |
| **I11** | **confirmado** | trigger pela aceleração; impactos múltiplos em **15 de 16** ensaios de acrílico, 2/16 Al, 1/16 aço. |

## Problemas / surpresas
- O duplo impacto quase sistemático no acrílico vai pesar no critério de rejeição (D1.2): rejeitar todos não é opção; a Fase 2 terá de tratar (janela de força no 1º pulso, ou FRF com o espectro de força real, que inclui os dois pulsos).
- C9 não se resolve com os dados: depende de você (registro do ensaio de jan/2023).

## Próximo passo
Etapa 0.8 — simulação no tempo, força, decremento logarítmico e efeito do desajuste de frequência (I5, I6, I9, C7).
