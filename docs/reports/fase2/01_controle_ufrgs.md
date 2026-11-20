# Fase 2 · Controle externo — pórtico UFRGS (base externa, não são dados do pesquisador)

**Base:** Barazzutti, Gomes & Drehmer, Mendeley Data v3 (2023), DOI 10.17632/npcn4c5htt.3, CC BY 4.0. Medições 51–100 (as usadas no artigo).

## O que foi feito

1. **Leitor `.daq`** (`src/dampid/io/mat_daq.py`): o formato do MATLAB DAQ Toolbox não tem documentação nem leitor em Python. A estrutura foi deduzida dos arquivos (blocos de 48 bytes de cabeçalho + amostras `uint16` intercaladas; volts = contagem × 10/4095 − 5).
   **Verificação:** as 50 medições reproduzem, amostra a amostra, o bloco correspondente do `Leituras.dat` que o próprio script dos autores escreveu. Erro máximo de 5·10⁻⁵ m/s², que é o arredondamento da impressão (`%12.4e`). Teste: `tests/unit/test_mat_daq.py`.
2. **Reprodução do resultado publicado** (`src/dampid/external/ufrgs.py`): port em Python do `processamento.m` dos autores (soma dos 3 canais → picos da FFT → filtro de ±15 Hz → envelope → ajuste exponencial).

| | f1 | f2 | f3 | f4 |
|---|---|---|---|---|
| f publicada [Hz] (média ± dp) | 36,43 ± 0,12 | 154,16 ± 0,42 | 246,22 ± 0,95 | 294,23 ± 0,67 |
| f nossa leitura + port [Hz] | 36,44 ± 0,12 | 154,19 ± 0,42 | 246,27 ± 0,95 | 294,29 ± 0,67 |
| ζ publicado [%] | 1,26 ± 0,53 | 0,57 ± 0,08 | 0,90 ± 0,15 | 1,17 ± 0,17 |
| ζ nossa leitura + port [%] | 1,26 ± 0,55 | 0,57 ± 0,08 | 0,90 ± 0,15 | 1,17 ± 0,17 |

Por medição: f difere ≤ 0,02 %; ζ dos modos 2–4 difere ≤ 4 % (mediana ≤ 0,6 %). O ζ do modo 1 difere até 23 % em 6 das 50 medições. A causa é a implementação do envelope (MATLAB `envelope(...,'peaks')` × spline nossa) com só ~18 ciclos na janela de 0,5 s.

## O que isso mostra

* Lemos os dados exatamente como os autores (mesmos volts, mesma conversão, mesmos números finais).
* **O ζ do 1º modo publicado depende do detalhe do algoritmo de envelope** na mesma ordem da sua própria dispersão (CoV 42 %). Isso é um exemplo concreto, com dados de terceiros, de que a incerteza de *método* precisa entrar no orçamento de incertezas (Fase 9), e não só a dispersão entre repetições.
* O método dos autores usa um limiar de pico **ajustado à mão por medição** (vetor `level` no script). Nosso método da Fase 3 deve funcionar sem ajuste manual; esta base é um teste para isso.

## Limitações (por que a base não basta)

* **Não há canal de força** (abertas: 0, 25, 51–100; os 101 arquivos têm o mesmo tamanho), embora a descrição cite o martelo Endevco 4416. Sem força não há FRF, e o simulador FE não pode ser comparado em amplitude nem no tempo.
* Escala dos acelerômetros divergente: a descrição dá 970 mV/g; o script usa 3,205 V/g (fator ≈ 3,3). Não afeta f e ζ. Com 970 mV/g, o canal 2 chegaria a ≈ 1,96 g, acima da faixa de ±1,7 g do ADXL203. Não há patamar de saturação no conversor (nenhuma medição tem ≥ 5 amostras no valor extremo).
* Pórtico com juntas parafusadas e cantoneiras de PLA, base presa numa mesa: não é a viga do modelo atual.

## Uso futuro

Referência *output-only* de f e ζ para a Fase 3 (EMA sem ajuste manual, comparada à tabela publicada) e para a discussão de incerteza de método (Fase 9).
