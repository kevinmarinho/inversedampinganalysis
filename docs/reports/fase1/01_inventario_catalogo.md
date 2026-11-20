# Relatório 1.1 — Inventário e catálogo

**Bases de dados:** dados do pesquisador — 156 CSVs em `data/raw/` (integridade: `MANIFEST.sha256`). Metadados de hardware: arquivos `.task` do NI (pesquisador) e especificações do NI 9234 e dos sensores PCB (externas, só para fundo de escala e massa nominal).
**Notebook:** `notebooks/fase1/01_inventario_catalogo.ipynb` · **Saídas:** `data/catalog.yaml`, `01_inventario_configuracoes.csv`, `01_catalogo_arquivos.csv`
**Código:** `dampid.io.ni_csv` (leitor), `dampid.io.catalog` (modelos pydantic + gerador), `configs/data_sources.yaml` (metadados estáticos por pasta)

## Feito
- Baixei do Projeto os 96 CSVs que faltavam: fibra de carbono, fibra de vidro alinhada e a 90°, PLA caverna e PLA longarina nas duas configurações. O `data/raw/` agora tem os **156 arquivos** e o manifesto foi ampliado.
- Leitor novo (`read_ni_csv`): lê cada canal com a sua coluna de tempo, confere se o tempo avança exatamente `dt` e registra qualquer quebra.
- Gerei o catálogo: uma entrada validada por arquivo, com campanha, configuração, espécime, condição de contorno, repetição, data/hora, fs, amostras, duração, canais (nome, unidade, papel, sensor, sensibilidade configurada, fundo de escala) e hash SHA-256.

## Evidência
| Conferência | Resultado |
|---|---|
| Unidades por posição de canal | `g, g, N` em **156/156** arquivos |
| Taxa de amostragem | 25 600 Hz em 156/156 |
| Duração | 2ª campanha: 5,0 s, exceto `PLA_LONGARINA*` (15 de 16 com **3,0 s**; só a repetição 1 tem 5 s). 1ª campanha: 2–6 s, conforme o ensaio |
| Datas | 2ª campanha em duas sessões: **21/11/2024** (aço, Al, acrílico, FV 0° e 90°) e **28/11/2024** (fibra de carbono, PLA caverna, PLA longarinas). 1ª campanha: 25/01/2023 |
| Nomes de canal | 2ª campanha: `Acelerômetro 1 / Acelerometro_1 / Acelerometro 1`, `Acelerometro 2`, `Martelo`. 1ª campanha: `Input 0/1/2` |
| Uniformidade do tempo | desvio máximo de `dt` < 10⁻¹⁰ em todos os arquivos, menos um (abaixo) |

## Problemas / surpresas
- **`exp_5/Al_engastado_livre_mponta_aponta.csv`**: depois das 102 400 amostras (4,0 s) há **2 193 linhas de um segundo bloco colado** no fim. O tempo volta para 3,914 s, e a primeira linha desse bloco tem tempo 2,76·10¹⁴. O leitor corta o arquivo na quebra e registra o descarte em `anomalies`. O legado lia essas linhas como se fossem dados.
- **`exp_5/Aco_engastado_livre_mponta_aponta.csv`** tem carimbo de data 31/12/1903: o relógio não estava configurado. Não afeta os dados.
- **PLA longarina**: só a primeira repetição tem 5 s. As outras 15 têm 3 s e, como o golpe ocorre por volta de 1 s, sobram ~1,9 s de resposta (Δf ≈ 0,5 Hz).
- A 2ª campanha **não registra** posição dos sensores, ponto de impacto nem ponteira. Esses campos ficam `null` no catálogo, nunca preenchidos por palpite.

## Próximo passo
Etapa 1.2 — unidades.
