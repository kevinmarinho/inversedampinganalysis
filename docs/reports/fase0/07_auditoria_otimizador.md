# Relatório 0.8 — Simulação no tempo, força, decremento e desajuste de frequência (I6, N4, I5, I9, C7)

**Bases de dados:** dados do pesquisador (`AÇO/` livre-livre: pulso de força e 1º modo; `exp_5/`: α, β do legado); casos sintéticos com parâmetros conhecidos; modelo legado (K, M) com ζ = 0,5 % correto.
**Notebook:** `notebooks/fase0/07_auditoria_otimizador.ipynb`

## Evidência
| Item | Status | Evidência |
|---|---|---|
| **I6** | **confirmado** | 0,5 s de sinal: Radau (tolerância padrão) 1,0 s / 21 696 avaliações; tolerância apertada 5,3 s; solução exata (expm, ZOH) 0,07 s. Mesma diferença de ~2 % na velocidade nos dois Radau — o integrador só adiciona custo. |
| **N4 (novo)** | **confirmado** | 15 dos 30 modos do modelo de 6 nós estão acima de Nyquist (até 244 kHz). Amostrados à taxa do ensaio, dobram para dentro da banda: a aceleração derivada difere **47 %** da exata entre 0–500 Hz; sem os modos > 2 kHz, **0,9 %**. O sensor real tem anti-aliasing; o modelo legado não. |
| **I5** | **confirmado** | Força real (aço): −10 dB em 206 Hz, −20 dB em 312 Hz. Retângulo de 1 ms com o mesmo impulso: −10 dB em 756 Hz, −20 dB em 930 Hz. Energia distribuída de forma completamente diferente na banda. |
| **I9** | **confirmado** | Dados **sintéticos** de ζ conhecido (0,2 % e 0,3 %) passados pela própria `pos_experimental`: ζ₁ = 2,8 % (13,8×) e ζ₂ **negativo**. Nos arquivos de `exp_5/`, o Rayleigh do legado corresponde a ζ físico de 3 % a 1500 % entre 50 e 300 Hz. |
| **C7** | **confirmado** | 1 GL sintético (ζ = 0,2 %): erro de frequência de 0,5 % → ζ ajustado 2×; 1 % → 3,7×; 2 % → 8,5×; 5 % → 43×. Aço real, 1º modo: ζ = 0,0011 com a frequência certa; 0,029 com a do modelo da Tabela 1 (−3 %); 0,45 com a do código (+23 %). |

## Conclusão parcial importante
C2 (escala de C) e C7 (desajuste de frequência) **sozinhos** já explicam a conclusão da monografia de que "o modelo superestima o amortecimento", sem precisar invocar limitações dos otimizadores — sustenta H1 e H2 do projeto e a D0.6.

## Próximo passo
Etapa 0.9 — consolidar `docs/audit_legacy.md` e conferir as afirmações do TCC e do relatório FAPESP.
