# Relatório 1.3b — Qualidade dos canais de aceleração

**Bases de dados:** dados do pesquisador (144 registros da 2ª campanha; 1ª só na tabela de saturação). Fundo de escala: ±5 V / 100 mV/g = ±50 g (NI 9234; faixa do PCB 352C33 também ±50 g pk).
**Notebook:** `notebooks/fase1/04_qualidade_resposta.ipynb` · **Tabelas:** `04_resposta_por_config.csv`, `04_resposta_por_registro.csv`
**Código:** `dampid.eda.quality` (`clipping`, `response_quality`, `end_rms_split`), com testes sintéticos

## Evidência
**Saturação (|a| ≥ 49 g; o NI 9234 corta em ±51,1–51,2 g):**

| Config | Acel 1 | Acel 2 | os dois | o recorte termina até |
|---|---|---|---|---|
| Aço | 1 (rep. 11) | 1 | 1 | 81 ms após o golpe |
| Alumínio | 1 (rep. 13) | 1 | 1 | 7 ms |
| Acrílico | 5 (reps. 02, 03, 08, 13, 14) | 5 | 5 | 6 ms |
| **Fibra de carbono** | **16/16** | **16/16** | **16** | 2–30 ms |
| FV 90° | 0 | 1 (rep. 16) | 0 | — |
| PLA longarina 0 / 90 | 1 / 1 (Acel 1) | 0 / 0 | 0 | < 1 ms |

A fibra de carbono saturou em todos os ensaios. É um tubo leve, rígido e pouco amortecido, golpeado com ~40 N por uma ponteira de contato curto (3 ms): a aceleração de pico passa de 50 g. **Para a FRF, essa configuração está perdida.** Depois do último recorte (≤ 30 ms), o decaimento livre ainda é linear e pode servir a métodos no tempo (Fase 3). O ideal seria repetir o ensaio com golpe mais fraco ou ponteira mais macia.

**Ruído, SNR, decaimento (medianas, Acel 1):**

| Config | SNR [dB] | decaimento no registro [dB] | fim da banda elástica acima do piso [dB] | pré-impacto acima do piso [dB] |
|---|---|---|---|---|
| Aço | 90 | **−25** | **69** | **33** |
| Alumínio | 87 | −68 | 8 | 17 |
| Acrílico | 84 | −66 | 8 | 16 |
| Fibra de carbono | 94 | −71 | 3 | 16 |
| FV 0° / 90° | 76 / 85 | −53 / −63 | 14 / 11 | 14 / 18 |
| PLA caverna | 81 | −61 | 3 | 15 |
| PLA longarina 0 / 90 | 80 / 81 | −53 / −62 | 10 / 10 | 15 / 17 |

- Piso de ruído: **0,50 mg rms** em banda larga (janelas de 10 ms) e 0,34 mg rms acima de f_c. É igual nos dois canais e nas duas sessões.
- **Aço truncado:** em 5 s a resposta cai só ~25 dB e termina ~69 dB acima do piso na banda do modo elástico. A FFT com janela retangular vai ter vazamento, e isso pesa na D2.2 (janela exponencial com correção × métodos no tempo). Coerente com ζ muito baixo: um decaimento de ~25 dB (2,9 Np) em ~4 s a 158 Hz corresponde a ζ ≈ 2,9/(2π·158·4) ≈ 7·10⁻⁴. É só ordem de grandeza; o valor sai da Fase 3.
- **Aço já vibrava antes do golpe** (33 dB acima do piso, ~57 dB abaixo da resposta): o golpe anterior ainda não tinha morrido. Condição inicial não nula; 10/16 registros sinalizados.
- Nos demais espécimes, o que resta no fim do registro é quase todo **movimento lento da suspensão** (< f_c; 2,5–8 mg). A banda elástica chega a 3–14 dB do piso: sem truncamento relevante.
- **Offset DC** antes do golpe: 1,6–3,4 mg, desprezível (acoplamento AC, corte em 0,5 Hz).
- Amostragem uniforme em todos os arquivos (Etapa 1.1).

## Problemas / surpresas
- A saturação de 16/16 na fibra de carbono é nova e não aparece em nenhum relatório anterior.

## Próximo passo
Etapa 1.4 — consistência entre repetições e sensores.
