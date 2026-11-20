# Relatório 1.2 — Unidades e mapeamento de canais (inclui C9)

**Bases de dados:** dados do pesquisador (156 CSVs; `.task` de 21/11/2024; Tabela 1 do TCC). Externas, só para ordem de grandeza: massa do martelo PCB 086C03 (0,16 kg, pcb.com), massa do PCB 352C33 (5,8 g, datasheet rev. H), massas específicas nominais (aço 7750–8050, Al 2660–2810, PMMA 1170–1200 kg/m³; `configs/data_sources.yaml`).
**Notebook:** `notebooks/fase1/02_unidades.ipynb`

## Feito
O cabeçalho diz `g, g, N`, mas um cabeçalho só repete a configuração do software. Testei a escala por três caminhos físicos independentes, mais o C9.

## Evidência
**1. O piso de ruído é igual nas duas sessões** (21/11 × 28/11). Uma sensibilidade configurada de forma diferente mudaria o ruído em g ou N na mesma proporção.

| | 21/11 | 28/11 | razão |
|---|---|---|---|
| Acel 1 [mg rms] | 0,510 | 0,497 | 0,97 |
| Acel 2 [mg rms] | 0,475 | 0,504 | 1,06 |
| Força σ [N] | 0,027 | 0,029 | 1,10 |

As sessões de 28/11 (fibra de carbono, PLA) não têm `.task`. O ruído igual sustenta que a configuração foi a mesma.

**2. Impulso do martelo.** J = 0,02–0,09 N·s, o que dá velocidade de impacto de **0,07–0,55 m/s** (J/2m_h a J/m_h, m_h = 0,16 kg): um golpe manual leve, plausível. Com a força em outra escala (por exemplo, volts lidos como N), a velocidade sairia 444× menor. A força máxima chegou a **10 % do fundo de escala** (aço, repetição 11) e ficou abaixo de 3 % no resto: o martelo nunca saturou.

**3. Linha de massa (barras livre-livre).** κ = |a/F|·m tem de ficar entre ~0,5 e 7 para pontos sobre a barra (κ = 1 + 12 x_f x_r/L² na linha média).

| Config | |A| [(m/s²)/N], 30–70 Hz | m nominal + 2 sensores [kg] | κ |
|---|---|---|---|
| Aço, Acel 1 / 2 | 4,67 / 3,98 | 0,372–0,398 | 1,74–1,86 / 1,48–1,58 |
| Alumínio, Acel 1 / 2 | 14,50 / 12,25 | 0,143–0,163 | 2,07–2,36 / 1,75–2,00 |

- Dá praticamente o **mesmo κ nas duas barras** (mesma montagem). Um erro de ×9,81 levaria κ para ~0,2 ou ~18.
- **Unidades g e N confirmadas para o conjunto de desenvolvimento.**
- No acrílico **não há patamar de linha de massa**: f₁ ≈ 42 Hz fica perto da suspensão (15–25 Hz). O número da tabela do notebook não vale como checagem.

**4. Sinal da linha de massa.** Re(a/F) < 0 em aço, Al, FV 90°, PLA caverna e longarinas, e > 0 em acrílico, fibra de carbono e FV 0°. O sinal depende da geometria (x_f·x_r) e da orientação do sensor. É informação sobre a montagem para conferir com as fotos, não uma conclusão.

**5. C9 (1ª campanha).** Input 1 [g] e Input 2 [N] têm pulsos de **mesma forma** (correlação 0,984–0,997, sem defasagem). A razão In1/In2, porém, **varia 23 %** (1,07–1,31): fica em ~1,07–1,16 nos biengastados e em ~1,27–1,31 nos engaste-livre. No aço engaste-livre, o Input 2 só tem ruído (pico 4–6 σ).
- ~~H-a (mesmo sinal ligado nas duas entradas): rejeitada.~~ **Revisto no adendo (relatório 07):** no domínio da frequência, |In1/In2| é plana (1,125–1,198) em 5–300 Hz. Os dois canais carregam a mesma força; a razão de picos no tempo era distorcida pela deriva de base do Input 1.
- **H-b (acelerômetro no martelo + célula de força) e H-c (dois transdutores de força independentes): continuam abertas.** Em H-b, a massa implícita do martelo seria 78–96 g, menor que os 160 g do catálogo e variando com a montagem.
- A escala de força da 1ª campanha **continua desconhecida** → `exp_5/` segue "pendente" (regra do prompt master §13).
- O Input 0 da 1ª campanha satura (≥ 49 g) em `Acrilico_engastado_livre_mponta` e `Al_engastado_livre_mcentro`.

## Problemas / surpresas
- A linha de massa está contaminada em ~10–30 % pela rigidez do 1º modo elástico (a curva sobe para f₁). Para a checagem de unidade basta; para estimar massa (§2.6, Fase 5), o modelo completo vai ter de absorver esse efeito.
- O κ das barras (~1,5–2,4) com dois sensores **praticamente iguais** (H₂/H₁ ≈ 1 ± 10 % em toda a banda; Etapa 1.4) restringe as posições possíveis. Pela teoria, impacto no centro daria κ = 1, o que exigiria uma massa ~2× menor que a nominal; parece improvável, mas só as posições reais resolvem.

## Próximo passo
Etapa 1.3 — qualidade da força e da resposta.
