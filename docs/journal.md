# Diário de sessões

## 2026-10-31 — Fase 0 (setup e auditoria do legado)

**Decisões registradas:** D0.1 uv · D0.2 YAML + pydantic · D0.3 estrutura §8 + notebooks/relatórios por etapa + GitHub · D0.4 código em inglês, docs em pt-BR · D0.5 trilha A primeiro · D0.6 comparação de métodos inversos na Fase 8 (todas pela recomendação, aceitas pelo pesquisador).

**Informação nova do pesquisador:** sensores fixados por rosca numa base impressa em PLA, colada com cera no corpo de prova; posições nas figuras do TCC/relatório. Massa do PCB 352C33: 5,8 g (datasheet rev. H).

**Feito**
- Repositório, ambiente (spec), legado intacto com hashes, dados brutos com manifesto.
- 8 notebooks executados (`notebooks/fase0/00`–`07`) e 8 relatórios (`docs/reports/fase0/`).
- `docs/audit_legacy.md`: 20 achados confirmados, 2 parciais, 1 inconclusivo (C9), 4 novos (N1–N4).
- Pacote: `dampid.analytic.beam` (soluções fechadas) e `dampid.audit` (harness do legado + réguas de referência), 17 testes passando (+1 pulado sem SymPy).

**Resultados principais**
- C2 é mais grave que o previsto: o ζ do código é amplificado ~10⁵–10⁶ ⇒ os "ζ ≈ 10⁻⁶" da monografia são ζ físicos de 4 %–200 %.
- C7 quantificado: 1 % de erro em frequência ⇒ ζ ajustado ×3,7.
- N3: alumínio livre-livre (−17 %) explicável por ≈7,5 g por conjunto sensor nas pontas; aço (+3 %) não — dimensões precisam ser medidas.
- N4: modelo de 6 nós tem metade dos modos acima de Nyquist ⇒ aliasing na simulação temporal.

**Limitações do ambiente:** PyPI bloqueado pela política de rede (proxy 403) — sem `uv.lock`, SymPy, Jupyter, `geneticalgorithm2`. Contornos documentados em `docs/reports/fase0/00_setup.md`.

**Pendências do pesquisador:** ver `docs/audit_legacy.md` §5, mais: criar o repositório GitHub (D0.3), copiar o `prompt_master.md` para `CLAUDE.md`, rodar `uv lock && uv sync && uv run pytest` localmente.

---

## Estado atual (§14 do prompt master)

- **Fase atual:** 0 — concluída, **aguardando aprovação do gate** (revisão de `docs/audit_legacy.md`).
- **Última decisão registrada:** D0.6 (2026-10-31).
- **Dados:** 60 CSVs em `data/raw/` (exp_5, AÇO, Al, ACRILICO). Faltam no ambiente: fibra de carbono, fibra de vidro (0°/90°), PLA caverna/longarina (entram na Fase 1).
- **Pendências do pesquisador:**
  - [ ] C9: o que estava em cada entrada do DAQ na 1ª campanha e com qual sensibilidade
  - [ ] Massa de um conjunto acelerômetro + base PLA + cera; posições dos conjuntos e do impacto (ou reanexar os PDFs com figuras)
  - [ ] Dimensões medidas das barras; liga do aço e do Al; procedência do acrílico
  - [ ] Ponto de impacto fixo ou variável entre repetições; ponteira por espécime
  - [ ] Criar repositório GitHub privado e vazio (D0.3) e informar `usuário/nome`
- **Próximo passo:** aprovação do gate da Fase 0 → plano da Fase 1 (inventário e EDA) com as decisões D1.1 e D1.2.
