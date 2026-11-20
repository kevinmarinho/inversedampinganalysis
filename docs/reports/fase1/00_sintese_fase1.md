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

Código novo: `dampid.io` (leitor, unidades, catálogo), `dampid.eda` (qualidade, levantamento, consistência, classificação); 27 testes (+1 pulado sem SymPy). Decisões: D1.1, D1.2.
Regenerar tudo, em ordem: `python scripts/build_notebooks.py fase1`. O 01 recria o catálogo, o 02 gera `data/processed/fase1_survey.csv`, o 04 e o 05 geram as tabelas por registro e o 06 classifica e atualiza o catálogo.

## Resultados principais
1. **Catálogo completo**: 156 arquivos em `data/catalog.yaml`, validado por pydantic, com hash, metadados, qualidade e classe. Unidades `g, g, N` e fs = 25 600 Hz em todos.
2. **Unidades confirmadas** para o conjunto de desenvolvimento por três caminhos físicos: piso de ruído igual entre sessões, velocidade do martelo de 0,07–0,55 m/s e linha de massa com κ ≈ 1,5–2,4 nas duas barras metálicas.
3. **C9 continua aberto.** A hipótese "mesmo sinal nas duas entradas" foi rejeitada (razão In1/In2 varia 23 %); a escala de força da 1ª campanha segue desconhecida → `exp_5/` pendente.
4. **Saturação dos acelerômetros:** 16/16 na fibra de carbono, 5/16 no acrílico, 1/16 no aço e no Al. A fibra de carbono está perdida para FRF.
5. **Banda útil da força** f₋₁₀ ≈ 150–190 Hz nos metais e no FV: só o 1º modo elástico (152–158 Hz) fica bem excitado. Longarinas e acrílico têm banda mais larga (460–630 Hz).
6. **Aço**: muito pouco amortecido (ζ ~ 7·10⁻⁴ na estimativa grosseira); truncado em 5 s (decai só 25 dB) e já vibrando antes de cada golpe → D2.2.
7. **Repetibilidade**: aço e Al excelentes (γ² ≥ 0,99; f₁ ± 0,03 % e ± 0,23 %). FV 0° não repetível (ponto/direção do golpe ou rotação do tubo variando).
8. **Classificação (D1.2)**: 52 ouro, 69 utilizável, 23 descartar, 12 pendente. Conjunto de desenvolvimento (D1.1): **aço 14 ouro, Al 13 ouro**.

## Gate da Fase 1
> "Catálogo completo e unidades confirmadas para o conjunto de desenvolvimento."

- Catálogo: ✅ (156/156; campos que só o pesquisador tem ficam `null`, sem palpite).
- Unidades no conjunto de desenvolvimento: ✅ (seção 2 acima).
- **Aguardando a aprovação do pesquisador.** O trabalho está no branch `fase-1` do GitHub. O merge na `main` acontece depois da aprovação.

## Medições físicas e informações pendentes (pesquisador)
Ordem = impacto no projeto.
1. **Posições** de cada acelerômetro e do ponto de impacto, por espécime, e direção do golpe (fotos da montagem bastam). Explicam κ ≈ 2 e H₂/H₁ ≈ 1 nas barras; necessárias na Fase 4.
2. **Massas**: de cada espécime (balança de 0,01 g) e de um conjunto sensor + base PLA + cera + trecho de cabo. Sem elas, a Fase 5 usa só o prior nominal (§2.6).
3. **Dimensões medidas** (paquímetro, vários pontos) do aço e do Al. A Tabela 1 do TCC e o código legado divergem (I14).
4. **Ponteira do martelo** por espécime. Os contatos de 6 ms, 3 ms e 2 ms sugerem pelo menos duas ponteiras.
5. **FV alinhado (0°)**: o golpe/ponto de impacto mudou entre repetições? O tubo podia girar na suspensão?
6. **Fibra de carbono**: dimensões (comprimento, diâmetros, espessura, layup). Dá para repetir o ensaio com golpe mais fraco/ponteira macia?
7. **PLA**: o que diferencia `PLA_LONGARINA` de `PLA_LONGARINA_90`? Cotas e parâmetros de impressão (preenchimento, orientação) da caverna e da longarina.
8. **FV**: espessura de parede e layup do tubo.
9. **C9**: o que estava ligado em cada entrada do NI 9234 na 1ª campanha (jan/2023), com qual sensibilidade e por que o Input 2 do aço engaste-livre está sem sinal.
10. Ligas exatas do aço e do Al; certificados de calibração dos sensores, se existirem.

## O que muda nas próximas fases
- **Fase 2** lê `file_class` e `quality.channel_status` do catálogo (canais saturados nunca entram). D2.2 (janela exponencial) é obrigatória para o aço. D2.3 (banda útil) parte de f₋₁₀ por configuração.
- **Fase 3**: só o 1º modo elástico de aço e Al está bem excitado, então o ζ de referência será de um modo por barra, a menos que a coerência mostre mais.
