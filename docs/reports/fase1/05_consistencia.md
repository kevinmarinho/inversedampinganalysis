# Relatório 1.4 — Consistência entre repetições e entre sensores

**Bases de dados:** dados do pesquisador (144 registros da 2ª campanha). Nenhuma referência externa.
**Notebook:** `notebooks/fase1/05_consistencia.ipynb` · **Tabelas:** `05_consistencia_por_config.csv`, `05_frequencias_pico.csv`; por registro em `data/processed/fase1_consistency.csv`
**Código:** `dampid.eda.consistency` (FRF de um registro por repetição, H1 de conjunto, coerência entre repetições, picos, desvio de log|H|)

Estimador de EDA (janela retangular, sem trigger refinado, impactos múltiplos incluídos): as coerências aqui são **limites inferiores** do que a Fase 2 vai obter.

## Evidência
| Config | γ² mediana na banda (Acel 1 / 2) | fração com γ² ≥ 0,9 | pico principal [Hz] (média ± dp) | desvio típico entre rep. [dB] |
|---|---|---|---|---|
| **Aço** | 0,995 / 0,996 | 0,95 / 0,99 | **157,62 ± 0,05** | 0,35 |
| **Alumínio** | 0,997 / 0,998 | 0,96 / 0,93 | **152,01 ± 0,35** | 0,28 |
| Acrílico | 0,94 / 0,96 | 0,58 / 0,80 | 42,06 ± 0,13; 51,0 ± 0,5; 303 ± 3 | 0,6–0,8 |
| Fibra de carbono* | 0,95 / 0,96 | 0,77 / 0,73 | 280,8 (279,6–281,4) | 0,8–0,9 |
| **FV 0°** | **0,38 / 0,33** | **0,00** | — (sem pico com γ² ≥ 0,8) | **4,5–5,1** |
| FV 90° | 0,96 / 0,98 | 0,86 / 0,89 | 165,03 ± 0,39 | 0,7–1,2 |
| PLA caverna | 0,94 / 0,97 | 0,80 / 0,84 | 193,8 ± 4,2 (Acel 1); 186,4 ± 2,9 (Acel 2) | 0,8–1,1 |
| PLA longarina 0 | 0,995 / 0,992 | 0,92 / 0,94 | 486,15 ± 0,37 | 0,3–0,4 |
| PLA longarina 90 | 0,994 / **0,72** | 0,94 / **0,32** | 493,8 ± 0,5 | 0,3 / 1,9 |

\*Fibra de carbono saturada em todos os registros: a média inclui sinais cortados. Os números são descritivos.

- **Aço e alumínio são muito repetíveis.** γ² ≥ 0,99 de ~40 Hz até o fim da banda da força, e o 1º modo varia 0,03 % (aço) e 0,23 % (Al) entre repetições. A dispersão de 0,23 % do Al é o primeiro número de referência para o gate |Δf|/f < 1 % da Fase 5.
- **Fibra de vidro 0° não é repetível** (γ² ≈ 0,35; repetições até 11–13 dB diferentes). As frequências de pico coincidem (≈ 161 Hz), mas os níveis mudam. É o padrão de **ponto de impacto, direção do golpe ou rotação do tubo mudando entre repetições** (num tubo, os dois planos de flexão têm a mesma frequência). Não dá para fazer média entre essas repetições; cada uma só vale isolada.
- **PLA longarina 90°, Acel 2:** γ² 0,72. O sensor a 90° quase não vê o movimento excitado (nível ~10× menor que o do Acel 1).
- **Acel 1 × Acel 2:** em aço, Al e fibra de carbono, |H₂/H₁| ≈ 1 ± 10 % em toda a banda: os dois sensores estão em **pontos mecanicamente equivalentes** (simétricos ou muito próximos). No FV 90°, |H₂/H₁| ≈ 1,5 constante. Na caverna, a razão varia com a frequência (pontos diferentes, útil para forma modal).
- **Repetições que destoam** (canal não saturado, desvio > max(1,5 dB; 3× a mediana)): acrílico 01 e 12 (Acel 2), PLA caverna 16 (Acel 2). Nenhuma no aço ou no Al.
- **O ponto de impacto parece fixo** em todas as configurações, exceto FV 0°: o sinal de Re(a/F) na linha de massa não troca entre repetições (16/16 iguais). Exceções: FV 0° (15/16 e 13/16), consistente com o item acima, e o Acel 2 da longarina 90° (14/16), canal de pouco sinal.
- Picos de baixa frequência (≈ 15–55 Hz) em quase todas as configurações estão na faixa dos modos de suspensão (barbantes). Não são modos elásticos. As frequências de pico por repetição desses modos variam muito (até 30 %), como esperado de uma suspensão.

## Problemas / surpresas
- O prompt master (§2.3) cita que o relatório fala em impactos aplicados em vários pontos. Os dados indicam **ponto fixo** por configuração, exceto no FV alinhado. **Confirmar com o pesquisador.**
- PLA longarina: o 1º modo elástico forte está em ~486 Hz (0°) e ~494 Hz (90°), justamente na borda da banda útil da força (f₋₁₀ ≈ 460–480 Hz).

## Próximo passo
Etapa 1.5 — classificação (D1.2) e síntese.
