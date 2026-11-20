# Fase 1 — Inventário e exploração dos dados: síntese

**Bases de dados usadas na fase:** só **dados do pesquisador** (156 CSVs das duas campanhas, `.task` do NI, Tabela 1 do TCC). Referências externas apenas para metadados de hardware (NI 9234: ±5 V e corte AC em 0,5 Hz; PCB 086C03: 2,25 mV/N, ±2224 N, 0,16 kg; PCB 352C33: ±50 g, 5,8 g) e massas específicas nominais da literatura, estas só para a checagem de ordem de grandeza da Etapa 1.2.

| Etapa | Notebook | Relatório |
|---|---|---|
| 1.1 Inventário e catálogo | `notebooks/fase1/01_inventario_catalogo.ipynb` | `01_inventario_catalogo.md` |
| 1.2 Unidades e C9 | `02_unidades.ipynb` | `02_unidades.md` |
| 1.3a Força | `03_qualidade_forca.ipynb` | `03_qualidade_forca.md` |
| 1.3b Aceleração | `04_qualidade_resposta.ipynb` | `04_qualidade_resposta.md` |
| 1.4 Consistência | `05_consistencia.ipynb` | `05_consistencia.md` |
| 1.5 Classificação | `06_classificacao.ipynb` | `06_classificacao.md` |
| 1.6 Adendo: montagem, ponteiras, linha de massa, C9 | `07_adendo_montagem.ipynb` | `07_adendo_montagem.md` |

Código novo: `dampid.io` (leitor, unidades, catálogo), `dampid.eda` (qualidade, levantamento, consistência, classificação); 27 testes (+1 pulado sem SymPy). Decisões: D1.1, D1.2, D1.3 (dimensões = Tabela 1 do TCC), D1.4 (canal de força da 1ª campanha).
Regenerar tudo, em ordem: `python scripts/build_notebooks.py fase1`. O 01 recria o catálogo, o 02 gera `data/processed/fase1_survey.csv`, o 04 e o 05 geram as tabelas por registro e o 06 classifica e atualiza o catálogo.

## Resultados principais
1. **Catálogo completo**: 156 arquivos em `data/catalog.yaml`, validado por pydantic, com hash, metadados, qualidade e classe. Unidades `g, g, N` e fs = 25 600 Hz em todos.
2. **Unidades confirmadas** para o conjunto de desenvolvimento por três caminhos físicos: piso de ruído igual entre sessões, velocidade do martelo de 0,07–0,55 m/s e linha de massa com κ ≈ 1,5–2,4 nas duas barras metálicas.
3. **C9 resolvido:** o martelo foi ligado em duas entradas (pesquisador). Input 1 e Input 2 são a mesma força, com |In1/In2| plana em 5–300 Hz e valor 1,13–1,20; Input 0 é a resposta. D1.4: força = Input 2; no aço engaste-livre, Input 1 ÷ 1,135. Escala absoluta ±15 %; ζ não é afetado.
4. **Saturação dos acelerômetros:** 16/16 na fibra de carbono, 5/16 no acrílico, 1/16 no aço e no Al. A fibra de carbono está perdida para FRF.
5. **Banda útil da força** f₋₁₀ ≈ 150–190 Hz nos metais e no FV: só o 1º modo elástico (152–158 Hz) fica bem excitado. Longarinas e acrílico têm banda mais larga (460–630 Hz).
6. **Aço**: muito pouco amortecido (ζ ~ 7·10⁻⁴ na estimativa grosseira); truncado em 5 s (decai só 25 dB) e já vibrando antes de cada golpe → D2.2.
7. **Montagem real das barras livre-livre:** sensores perto das extremidades, impacto perto do meio (inferido; o pesquisador não fixou o ponto). O nível medido fica 1,2–1,4× acima do modelo com massa nominal → pesar os espécimes.
8. **Repetibilidade**: aço e Al excelentes (γ² ≥ 0,99; f₁ ± 0,03 % e ± 0,23 %). FV 0° não repetível (ponto/direção do golpe ou rotação do tubo variando).
9. **Classificação (D1.2 + D1.4)**: 52 ouro, 79 utilizável, 25 descartar. Conjunto de desenvolvimento (D1.1): **aço 14 ouro, Al 13 ouro**. 1ª campanha: 10 utilizáveis, 2 descartados (acelerômetro saturado).

## Gate da Fase 1
> "Catálogo completo e unidades confirmadas para o conjunto de desenvolvimento."

- Catálogo: ✅ (156/156; campos que só o pesquisador tem ficam `null`, sem palpite).
- Unidades no conjunto de desenvolvimento: ✅ (seção 2 acima).
- **Aguardando a aprovação do pesquisador.** O trabalho está no branch `fase-1` do GitHub. O merge na `main` acontece depois da aprovação.

## Medições físicas e informações pendentes (pesquisador)
Respondidas em 2026-11-20 (registradas no catálogo e em `configs/data_sources.yaml`): posições e orientação dos sensores, ponto de impacto (livre-livre: sem ponto fixo; 1ª campanha: L/4), dimensões (Tabela 1), ponteiras (metal/borracha → inferidas pelo contato), FV alinhado (golpe e rotação variaram), longarina 0/90 (orientação dos sensores), ligas e certificados (não aplicáveis; sensibilidades nominais).

Ainda pendentes:
1. **Massas** de cada espécime e de um conjunto sensor + base PLA + cera + cabo. É o item de maior impacto: decide se a diferença de 1,2–1,4× na linha de massa vem da massa ou das sensibilidades.
2. Posição do impacto **por arquivo** da 1ª campanha (os nomes `mcentro`/`mponta`/`mforacentro` não batem todos com "L/4") e comprimento livre no engaste.
3. Fibra de carbono: dimensões e layup; dá para repetir sem saturar?
4. Tubo de FV: espessura de parede e layup. PLA: cotas e parâmetros de impressão.
5. Qual acelerômetro (1 ou 2) está em cada extremidade nas barras (só importa se a simetria for quebrada).

## O que muda nas próximas fases
- **Fase 2** lê `file_class`, `force_channel`/`force_scale_to_N` e `quality.channel_status` do catálogo (canais saturados nunca entram). D2.2 (janela exponencial) é obrigatória para o aço. D2.3 (banda útil) parte de f₋₁₀ por configuração.
- **Fase 3**: só o 1º modo elástico de aço e Al está bem excitado, então o ζ de referência será de um modo por barra, a menos que a coerência mostre mais.
