# Relatório 0.1 — Setup do repositório e integridade

**Base de dados:** dados do pesquisador (60 CSVs: `exp_5/`, `AÇO/`, `Al/`, `ACRILICO/`), só leitura.
**Notebook:** `notebooks/fase0/00_setup.ipynb`

## Feito
- Repositório git com a estrutura da §8 (+ `analytic/`, `audit/`), `pyproject.toml` (uv; grupos `dev` e `legacy`), `ruff` e `pytest` configurados.
- Legado copiado para `legacy/` com SHA-256 registrado; dados brutos em `data/raw/` sem permissão de escrita, com `MANIFEST.sha256`.
- Decisões D0.1–D0.6 registradas em `docs/decisions/`.
- Ferramenta `scripts/build_notebooks.py`: executa os notebooks e embute saídas e figuras.

## Evidência
- Hashes do legado e dos 60 CSVs conferidos no notebook (todos íntegros).
- 12 testes passam (raízes de Euler–Bernoulli contra tabelas; constante de torção; amortecimento modal recuperado; meia potência; H1; elemento legado contra entradas em forma fechada); 1 teste pulado (comparação com SymPy, que não está instalado).

## Problemas / surpresas
1. **PyPI bloqueado pela política de rede do ambiente** (proxy respondeu 403 para `pypi.org`). Consequências: `uv.lock` não foi gerado aqui; SymPy, Jupyter, `geneticalgorithm2` e `npTDMS` não estão disponíveis. Contornos: notebooks executados por script próprio; `matrixes()` do legado substituída por reprodução numérica verificada (ver `legacy/README.md`). **Ação na sua máquina:** `uv lock && uv sync && uv run pytest` — o teste pulado compara a reprodução com o SymPy original.
2. **1ª campanha:** o cabeçalho de `exp_5/` diz `Input 0 (g) | Input 1 (g) | Input 2 (N)`. A força é o 3º canal, mas o legado lê o 2º → reforça C9 (notebook 06).
3. Um arquivo de `exp_5/` (`Aco_engastado_livre_mponta_aponta.csv`) tem data `31/12/1903` no cabeçalho: relógio do DAQ não ajustado; sem impacto nos sinais.
4. Durações da 1ª campanha variam de 2 a 6 s por arquivo (2ª campanha: 5 s fixos).

## Próximo passo
Etapa 0.2 — mapa do legado e execução dos scripts como estão.
