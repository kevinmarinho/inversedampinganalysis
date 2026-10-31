# Relatório 0.6 — Elemento, propriedades e dimensões (I3, I4, I14, C8)

**Bases de dados:** dados do pesquisador (CSVs livre-livre de aço, Al e acrílico, só para identificar frequências); Tabela 1 do TCC (dimensões); datasheet PCB 352C33 (5,8 g, fonte externa do fabricante); faixa de PMMA do prompt master (**a confirmar**); soluções analíticas (Blevins; Timoshenko & Goodier).
**Notebook:** `notebooks/fase0/05_auditoria_elemento.ipynb`

## Evidência
| Item | Status | Evidência |
|---|---|---|
| **I14** | **confirmado** | L do código é 29–33 mm menor que a Tabela 1 nos metais (provável comprimento livre da 1ª campanha). Espessura do aço: 3,04 mm no código × 2,90 mm no TCC (4,8 %, ≈ 5 % em frequência). Nenhuma das fontes é medição documentada. |
| **I3** | **confirmado** | J = Iy+Iz ≈ 68–74 × J_t de Saint-Venant. 1ª torção livre-livre: legado 3,3–5,8 kHz; correto ≈ 700 Hz (Al), 640 Hz (aço), 400 Hz (acrílico, dentro da banda). |
| **I4** | **confirmado, desprezível** | termo de inércia rotacional errado (forma e unidade), mas Δf ≤ 0,08 %. |
| **C8** | **confirmado** | acrílico com E = 20 GPa / ρ = 2000 → f₁ livre-livre ≈ 108 Hz, sem pico medido perto; com PMMA típico, 48–57 Hz, na ordem dos picos medidos (42 e 51 Hz). |

**Frequências medidas (H1 simples, 16 ensaios):** aço 157,5 e 439,6 Hz (razão 2,79 ≈ 2,757 teórica dos dois primeiros modos de flexão); alumínio 151,7 Hz; acrílico 42,1 / 50,8 / 304,6 Hz.

## Achado novo N3 (aberto)
- **Alumínio:** previsto sem sensores 182 Hz × medido 151,7 Hz (−17 %). Dois conjuntos sensor+base de ≈ 7,5 g nas pontas explicariam — compatível com 5,8 g do acelerômetro + ~2 g de base de PLA, cera e cabo. Hipótese, depende das posições reais.
- **Aço:** previsto sem sensores 153 Hz (espessura 2,90 mm) × medido 157,5 Hz (+3 %). Massa adicionada só reduziria a frequência: **não explica**. A espessura (ou E) do aço precisa ser medida.

## Pendências físicas que este notebook tornou concretas
- Massa de **um conjunto** acelerômetro + base de PLA + cera (+ cabo próximo); posições dos dois conjuntos e do ponto de impacto (as figuras do Apêndice A não vieram no texto do PDF).
- Comprimento, largura e espessura **medidos** de cada barra (vários pontos); liga do aço.

## Próximo passo
Etapa 0.7 — sinais e unidades (C5, C6, C9, I7, I8, I10, I11).
