# Auditoria do código legado e da monografia — Fase 0

**Data:** 2026-10-31 · **Status:** aguardando revisão do pesquisador (gate da Fase 0)

**Bases de dados usadas na auditoria (todas identificadas em cada notebook):**
- **Dados do pesquisador:** 60 CSVs — `exp_5/` (1ª campanha, 12) e livre-livre `AÇO/`, `Al/`, `ACRILICO/` (16 cada); arquivo `Analog Input_2.task` (configuração do DAQ, 2ª campanha).
- **Documentos do pesquisador:** TCC (2025) e Relatório Parcial FAPESP (jul/2025) — texto extraído (as figuras não vieram).
- **Fonte externa:** datasheet PCB 352C33 rev. H (massa 5,8 g).
- **Casos sintéticos** com parâmetros conhecidos e **soluções analíticas** (Blevins; Timoshenko & Goodier; Ewins).
- Nenhuma base experimental externa foi usada.

**Como o legado foi executado:** sem alteração, via `src/dampid/audit/legacy_harness.py`. Única substituição: `matrixes()` (SymPy, não instalável no ambiente) por reprodução numérica linha a linha, verificada contra as entradas em forma fechada do elemento (`tests/legacy/`).

---

## 1. Resumo

| Status | Itens |
|---|---|
| **Confirmados** | C1, C2 (mais grave do que o previsto), C3, C5, C6, C7, C8, I1, I2, I3, I4 (efeito desprezível), I5, I6, I7, I8, I9, I10, I11, I14, I15 |
| **Parciais / reformulados** | C4 (modos altos não ficam com ζ = 0: recebem ζ descontrolado), I12 (posição confirmada; malha suficiente para metais em balanço) |
| **Inconclusivo** | C9 (mapeamento de canais da 1ª campanha — precisa do registro do ensaio) |
| **Não avaliado aqui** | I13 (massa dos sensores no modelo) — é modelagem da Fase 4; a sensibilidade foi quantificada em N3 |
| **Novos** | N1 (RVAC quebra por tamanho de vetor), N2 (livre-livre quebra com max_time = 4 s), N3 (alumínio −17 % explicável por massa; aço +3 % não), N4 (aliasing de modos acima de Nyquist) |

**Conclusão principal.** Os resultados de otimização da 1ª campanha não podem ser usados: a função objetivo MAC era constante (C1), o parâmetro "ζ" do código era amplificado ~10⁵–10⁶ vezes (C2), a resposta numérica era uma rotação (C3), e qualquer desajuste de frequência infla o ζ ajustado (C7). Isso explica, **sem invocar limitações dos otimizadores**, a conclusão da monografia de que "o modelo superestima o amortecimento" — sustenta as hipóteses H1 e H2 e a decisão D0.6.

---

## 2. Achados — evidência, impacto e ação

| # | Achado | Status | Evidência (notebook) | Impacto | Ação (fase) |
|---|---|---|---|---|---|
| C1 | MAC entre vetores de frequências, sem dependência de ζ | confirmado | MAC = 0,980924058977776 para ζ de 10⁻⁷ a 0,5; listas sem relação dão 0,83 (02) | Gradiente nulo; Tabela 4 do TCC ("parou no chute") é consequência disso, não de mínimos locais | MAC só para parear formas modais reais (F4/F5) |
| C2 | `C = Φ diag(ζ) Φᵀ` | **confirmado, grave** | ζ de entrada 10⁻⁶ → ζ real 0,04–2; 10⁻² → modelo superamortecido; correto recupera 1 % (03) | Todos os ζ do legado sem significado físico; "1 % ζ = 10⁻⁶" | `C = MΦ diag(2ζω) ΦᵀM` com teste (F6) |
| C3 | Resposta em `w'` (rotação) no engaste-livre | confirmado | GL 34 = nó 5, `w'`; razão rot./transl. 5,08 e 17,64 m⁻¹ = analítica (04) | Comparação com acelerômetro inválida, fator diferente por modo | Seleção por posição/direção (`dofs.py`, F4) |
| C4 | Modos além do 6º sem amortecimento | parcial | Recebem ζ de ~0 a ≈ 9 por vazamento (03) | Amortecimento descontrolado nos modos altos | Termos residuais / truncamento consciente (F7) |
| C5 | Unidades e fatores 2/N, 9,81, 9810, 9810² | confirmado | 2/N = 1,6·10⁻⁵ na FRF; fatores manuais (06) | Escala da FRF errada; escondida por métricas invariantes | Conversões nomeadas, acelerância em m/s²/N (F1/F2) |
| C6 | Butterworth na magnitude da FRF | confirmado | ζ = 0,1 %: +103 % (Wn 0,5) e +3800 % (Wn 0,02 do TCC) (06) | Viés para cima, maior em baixo amortecimento | Proibido; médias de espectros + coerência (F2) |
| C7 | Rigidez/massa não atualizadas antes de ζ | confirmado | 1 % de erro em f → ζ × 3,7; aço real: 0,0011 → 0,45 com f do código (07) | Explica "superestima o amortecimento" | Fase 5 antes da Fase 8 (H2) |
| C8 | Acrílico com E = 20 GPa, ρ = 2000 | confirmado | f₁ livre-livre 108 Hz vs. picos medidos 42/51 Hz; PMMA típico 48–57 Hz (05) | Modelo do acrílico inválido | Propriedades com fonte + atualização (F5) |
| C9 | Canal da força na 1ª campanha | **inconclusivo** | Aço engaste-livre: Input 2 [N] só ruído, pulso no Input 1; outros 10: pulso nos dois (06) | Escala de força da 1ª campanha desconhecida | Registro do ensaio (pesquisador); `exp_5/` fora de qualquer uso de amplitude |
| I1 | 18 retornos × 17 variáveis | confirmado | `ValueError` ao rodar (01) | Script de otimização não roda; resultados vieram de outra versão | — (legado congelado) |
| I2 | Sem livre-livre no modelo | confirmado | `UnboundLocalError` (01) | 2ª campanha sem modelo | CC livre-livre (F4) |
| I3 | J = Iy + Iz na rigidez torcional | confirmado | J/J_t ≈ 68–74; 1ª torção 3–6 kHz → correto 400–700 Hz (05) | Sem modos de torção na banda | J_t de Saint-Venant (F4) |
| I4 | ρI∫NNᵀ em vez de ρI∫N′N′ᵀ | confirmado, desprezível | Δf ≤ 0,08 % (05) | Nenhum prático | Corrigir na F4 |
| I5 | Pulso retangular de 1 ms | confirmado | −20 dB: real 312 Hz × retângulo 930 Hz (07) | Excitação errada no tempo | Força medida como entrada (F7) |
| I6 | Radau no laço | confirmado | 15–80× mais lento que a solução exata, sem ganho (07) | Custo | Discretização exata / superposição modal (F7) |
| I7 | ω = f/2π em `frf_direta` | confirmado | pico em 39,5 × fₙ (06) | Função inutilizável | — |
| I8 | Fase ×½ | confirmado | −1 → 90° (06) | Fases erradas | — |
| I9 | Decremento logarítmico | confirmado | sintético ζ 0,2 % → 2,8 %; 0,3 % → negativo (07) | Estimativas iniciais inválidas | EMA (F3) |
| I10 | Média de magnitudes e fases | confirmado | figura de fases (06) | FRF média sem sentido físico | H1/H2 + coerência (F2) |
| I11 | Trigger pela aceleração; duplo impacto | confirmado | trigger até 2 s antes do golpe no aço; duplo impacto em 15/16 do acrílico (01, 06) | Cortes inconsistentes, FRFs contaminadas | Trigger pela força, detector, D1.2 (F1/F2) |
| I12 | Malha de 6 nós; centro não é nó | parcial | posição 0,4 L no engaste-engaste; erro de malha ≤ 0,05 % para metais em balanço na banda (04) | Ponto errado; malha OK para esse caso | Nós nos pontos medidos + convergência (F4) |
| I13 | Massa dos sensores não modelada | não avaliado (é F4) | sensibilidade em N3 | ver N3 | Massas concentradas (F4) |
| I14 | Dimensões código × TCC | confirmado | L 29–33 mm menor; espessura do aço 3,04 × 2,90 mm (05) | ≈5 % em frequência só pela espessura | Medir dimensões (pendência) |
| I15 | GA (0, 1e-5), 8 × 20 | confirmado | ótimo da RVAC em ζ_legado ≈ 1e-4, fora dos limites (02) | Busca não alcança o ótimo da própria métrica | D0.6 (F8) |
| N1 | RVAC: vetores 4000 × 1599 | novo, confirmado | erro de dimensão (02) | Métrica do `main_otimizador` não roda | — |
| N2 | `pos_experimental_livre_livre` com 4 s | novo, confirmado | FFT de 1987 × 2000 bins (01) | Quebra com trigger tardio | Pipeline novo (F2) |
| N3 | Frequências livre-livre × modelo | novo, aberto | Al: 182 → 151,7 Hz explicável por ≈ 7,5 g por conjunto nas pontas; aço: 153 → 157,5 Hz **não** explicável por massa (05) | Dimensões/massas decidem o gate da F5 | Medir; F5 |
| N4 | Modos acima de Nyquist na simulação temporal | novo, confirmado | 15/30 modos > 12,8 kHz; erro 47 % → 0,9 % ao truncar (07) | Comparação no tempo contaminada por aliasing | Domínio da frequência / truncamento + residuais (F7) |

---

## 3. Afirmações do TCC e do relatório FAPESP a corrigir no paper

Conferidas no texto extraído dos dois documentos (os trechos são idênticos nos dois).

| Afirmação no texto | O que a evidência mostra |
|---|---|
| SHELL181 é "de quatro nós" e a malha usa "elementos finitos quadráticos" (4.5.1) | Contraditório: SHELL181 é o elemento **linear** de 4 nós. |
| "amortecimento modal com um coeficiente uniforme de 1% ζ = 10⁻⁶" (4.3.4.2) | 1 % = 10⁻²; e o número do código não é ζ físico (C2). |
| "os métodos BDF foram adotados como a principal abordagem" (4.3.1.1) | O código usa `method='Radau'`. |
| RK e **HHT** avaliados "da biblioteca scipy.integrate" (4.3.1.1) | HHT não existe no `scipy.integrate`. |
| MAC calculado "a partir dos autovalores da matriz A amortecida" (4.4) | O código usa `eigh(K, M)` **sem** amortecimento (C1). |
| Seção 4.2.1 (engaste-livre do alumínio) diz "para este ensaio específico da viga de alumínio livre-livre" | O caso descrito é engaste-livre. |
| L-BFGS-B "convergiu para mínimos locais" (Tabela 4); GA limitado por hiperparâmetros (Tabela 5) | A função MAC era constante (C1); o ζ do código é amplificado ~10⁵ (C2): os ζ ≈ 3·10⁻⁷–10⁻⁵ da Tabela 5 (aço) correspondem a ζ físico de ~1 % a ~90 % nos modos de flexão fora do plano. Conclusão sobre os otimizadores deve ser revista. |
| Taxa de amostragem "25.560 Hz (intervalo de 39,12 µs)" (4.2.1.1) | Arquivos e `.task`: **25 600 Hz** (39,0625 µs). |
| Martelo "com sensibilidade de 100 mV/N" (4.2.1.1) | `.task` da 2ª campanha: **2,25 mV/N** configurados no DAQ. |
| Filtro do alumínio com `Wn = 0,02` (Algoritmo 4) | Código: `Wn = 1/2` para todos os materiais. |
| Eq. (6): termo de torção escrito como `GJ ∫ δu′u′` e J = Iyy + Izz | O termo deveria ser no ângulo de torção, e J (rigidez) é a constante de Saint-Venant (I3). |
| "Os impactos foram aplicados em múltiplos pontos da viga" (4.2.1.1) | A repetibilidade de f₁ e da linha de massa nos 16 ensaios sugere ponto fixo por espécime — **confirmar com o pesquisador**. |
| Objeto descrito como "placa" e modelado como viga | Coerente com a escolha declarada, mas o paper deve explicitar a hipótese de viga e verificar efeitos de placa (torção, largura/comprimento ≈ 1/6). |

---

## 4. Configuração confirmada do DAQ (2ª campanha, `Analog Input_2.task`, alumínio)
Input 0 e 1: acelerômetros, IEPE 2 mA, acoplamento AC, 100 mV/g. Input 2: martelo ("DynamicForce"), IEPE 2 mA, AC, **2,25 mV/N**, faixa configurada ±50 N. Taxa 25 600 Hz, sem trigger. (A faixa de ±50 N é só a escala configurada; o limite físico do NI 9234 a 2,25 mV/N é muito maior — saturação será checada na Fase 1.)

---

## 5. Pendências que dependem do pesquisador

1. **C9:** o que estava ligado em cada entrada do NI 9234 na 1ª campanha (jan/2023) e com qual sensibilidade.
2. **Montagem dos sensores:** massa de um conjunto acelerômetro + base de PLA + cera (+ trecho de cabo); posições dos dois conjuntos e do ponto de impacto em cada espécime (as figuras do Apêndice A não vieram no texto do PDF — anexar os PDFs de novo, ou fotos/desenhos).
3. **Dimensões medidas** de cada barra (comprimento, largura, espessura em vários pontos); confirmar se o L do código é o comprimento livre da 1ª campanha.
4. **Liga do aço e do alumínio**; procedência do acrílico.
5. **Ponto de impacto** fixo ou variável entre as 16 repetições; ponteira do martelo por espécime.
6. Arquivos originais do legado (para fechar a cadeia de custódia dos hashes — ver `legacy/README.md`).
