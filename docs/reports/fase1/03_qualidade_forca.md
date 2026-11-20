# Relatório 1.3a — Qualidade do canal de força

**Bases de dados:** dados do pesquisador (144 registros da 2ª campanha; a 1ª só descritiva). Fundo de escala do martelo: ±5 V / 2,25 mV/N (NI 9234, PCB 086C03).
**Notebook:** `notebooks/fase1/03_qualidade_forca.ipynb` · **Tabela:** `03_forca_por_config.csv`
**Código:** `dampid.eda.quality.find_pulses`, `force_bandwidth` (testados contra um meio-seno analítico: largura a 10 %, impulso 2FT/π e cruzamentos de −10/−20 dB)

## Evidência (medianas por configuração)
| Config | contato [ms] | pico [N] | f₋₁₀ [Hz] | f₋₂₀ [Hz] | impacto múltiplo (≥ 10 %) | 2º pulso: razão / atraso |
|---|---|---|---|---|---|---|
| Aço | 5,6 | 31 | 189 | 300 | 2/16 | — |
| Alumínio | 6,5 | 12 | 158 | 310 | 2/16 | — |
| Acrílico | 3,0 | 19 | 627 | 1062 | **16/16** | 0,90 / 5,6 ms |
| Fibra de carbono | 3,1 | 40 | 248 | 733 | **15/16** | 0,62 / 16 ms |
| FV 0° | 6,5 | 14 | 150 | 300 | 3/16 | — |
| FV 90° | 6,5 | 19 | 151 | 310 | 8/16 | 0,05 / 22 ms |
| PLA caverna | 5,8 | 6,7 | 154 | 192* | 12/16 | 0,27 / 12,5 ms |
| PLA longarina 0 | 2,2 | 32 | 464 | 619 | 2/16 | — |
| PLA longarina 90 | 2,1 | 36 | 480 | 636 | 5/16 | — |

\*Para um pulso quase meio-seno de duração T, o espectro tem um zero em 1,5/T (≈ 260 Hz na caverna) e cruza −20 dB logo antes. Por isso f₋₁₀ é o indicador mais estável.

- **Banda útil da força bem abaixo de 500 Hz nos metais e no FV:** f₋₁₀ ≈ 150–190 Hz e f₋₂₀ ≈ 300 Hz. O 1º modo (152–158 Hz) fica na borda de −10 dB; o 2º modo de flexão do aço (≈ 435 Hz) já está abaixo de −20 dB. A banda útil de identificação, nesses espécimes, provavelmente cobre só o 1º modo elástico (a confirmar pela coerência na D2.3).
- **Três famílias de contato** (5,6–6,8 ms; ~3 ms; ~2,1 ms) sugerem ponteiras diferentes ou a rigidez local da peça. A ponteira não foi registrada: **pendência do pesquisador**.
- **Golpes atípicos no aço:** a repetição 10 tem contato de 2,15 ms e dois pulsos; a 11 tem contato de 0,43 ms, 225 N, e saturou os dois acelerômetros (golpe de quina ou sem ponteira?).
- **Acrílico:** contatos múltiplos quase sistemáticos (2º pulso ≈ 90 % do principal, ~5,6 ms depois). O martelo quica no espécime leve. É impacto múltiplo de verdade, não um único contato com oscilação (a força volta a zero por mais de 0,5 ms entre os pulsos).
- **Fibra de carbono:** pulso principal com duas cristas (~1,5 ms entre elas) e um 2º impacto ~16 ms depois em 15/16 repetições.
- Força nunca saturou (máximo de 10 % do fundo de escala).

## Problemas / surpresas
- Impacto múltiplo **não invalida** a FRF H1 quando todos os pulsos estão no registro, mas cria vales no espectro de força (baixa coerência nessas frequências). Tratamento na Fase 2 (janela de força só no 1º pulso × espectro completo); classificação em D1.2.

## Próximo passo
Etapa 1.3b — canais de aceleração.
