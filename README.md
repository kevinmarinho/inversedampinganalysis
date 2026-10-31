# DampID — identificação inversa de amortecimento estrutural

Código e documentação do projeto de identificação de amortecimento a partir de ensaios
de impacto (UNESP – São João da Boa Vista; FAPESP 2023/14176-5). O roteiro completo, as
regras de colaboração e o diagnóstico herdado estão no prompt master (`CLAUDE.md`).

## Estado
Fase 0 — setup e auditoria do código legado. Ver `docs/journal.md` e `docs/audit_legacy.md`.

## Como usar
```bash
uv sync                      # cria .venv com todas as dependências (D0.1)
uv run pytest                # testes de verificação e de caracterização do legado
uv run python scripts/build_notebooks.py fase0   # re-executa os notebooks da Fase 0
```
Os dados brutos não estão no git: copie as pastas de `ultima versao/` para `data/raw/`
(ver `data/raw/README.md`) e confira com `sha256sum -c data/raw/MANIFEST.sha256`.

## Organização
| Pasta | Conteúdo |
|---|---|
| `legacy/` | código original da monografia, intacto (somente leitura) |
| `src/dampid/` | pacote novo (`analytic/`, `audit/` na Fase 0; demais módulos nas fases seguintes) |
| `tests/` | `unit/` (verificação contra soluções fechadas) e `legacy/` (caracterização) |
| `notebooks/<fase>/` | notebooks executados, um por etapa; fontes em `_src/` |
| `docs/reports/<fase>/` | relatório curto por etapa (+ figuras) |
| `docs/decisions/` | registro de decisões (D0.1 …) |
| `data/raw/` | dados brutos (fora do git, com manifesto de hashes) |

Cada notebook e relatório informa no topo **qual base de dados usou** (dados do
pesquisador, casos sintéticos analíticos ou bases externas).
