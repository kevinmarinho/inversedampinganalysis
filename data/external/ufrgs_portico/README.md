# Base externa: pórtico de alumínio da UFRGS (D1.5)

**Fonte:** Barazzutti, L. F.; Gomes, H. M.; Drehmer, L. R. C. *Data for: Sensitivity-based model update for estimating generalized proportional damping parameters in finite element model*. Mendeley Data, v3, 2023. DOI 10.17632/npcn4c5htt.3. Licença CC BY 4.0. Artigo: ASCE-ASME J. Risk Uncertainty Eng. Syst., DOI 10.1061/AJRUA6.RUENG-1178.

**Recebido do pesquisador em 2026-11-20** (`pacote_original_mendeley.zip` e, depois, a pasta `ExperimentalData/Conferencia MECSOL` já extraída do `.rar`, compartilhada pelo computador dele). Somente leitura; integridade em `MANIFEST.sha256`. Estes **não são dados do pesquisador**: toda análise com eles é marcada como base externa. Os binários não vão para o git (ver `.gitignore`); quem clonar o repositório baixa o pacote pelo DOI.

## Conteúdo

* `CAD.STEP`, `Detailing - Experimental Model.pdf` (desenho: 3 barras Al 300 × 19 × 5 mm, cantoneiras PLA, pés de aço, parafusos M3), `Materials - FE Model.png`, `Accelerometers.png` (posições), `Figure_1–4.jpg`, `ExperimentalData.rar` (original, não extraível aqui).
* `ExperimentalData/` (extraído pelo pesquisador; só a parte usada no artigo):
  * `datastream51.daq` … `datastream100.daq`: 50 medições, 3 acelerômetros (eixo z), MCC USB-1208FS (12 bits, ±5 V diferencial), 1000,1 Hz, 10 001 amostras (10 s). **Sem canal de força** (excitação manual, não medida).
  * `datastream.daq`: medição de 15:55:34 que não aparece no `Leituras.dat` (origem desconhecida; não usada).
  * `Leituras.dat`: exportação em texto (m/s²) feita pelo script dos autores; 75 blocos (inclui medições sem `.daq` correspondente nesta pasta).
  * `acquisition_3ch.m`, `processamento.m`: scripts originais. Conversão: a [m/s²] = (V − setp)/sens, setp = [2,090 2,211 2,150] V, sens = 3,205/9,81 V/(m/s²).
  * `ResultadosExperimentais.xlsx`, `results.txt`: f e ζ por medição publicados pelos autores (envelope + ajuste exponencial).

## Leitor e validação

`src/dampid/io/mat_daq.py` lê os `.daq` (formato do MATLAB DAQ Toolbox, sem documentação pública; estrutura descrita no módulo). Conferido em 2026-11-20: as 50 medições 51–100 reproduzem o bloco correspondente do `Leituras.dat` amostra a amostra, com erro máximo de 5·10⁻⁵ m/s² (= meia unidade da 4ª casa impressa, `%12.4e`); a correspondência medição ↔ bloco está em `tests/unit/test_mat_daq.py`.
