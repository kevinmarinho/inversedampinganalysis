# data/raw — dados brutos (somente leitura)

Os CSVs **não** são versionados no git (≈2 GB). A integridade é garantida pelo
`MANIFEST.sha256` (gerado com `sha256sum`). Para reconstruir: copie as pastas originais
de `ultima versao/` para cá, **mantendo os nomes** (`exp_5/`, `AÇO/`, `Al/`, `ACRILICO/`,
`FIBRA_CARBONO/`, `Fibra vidro_ alinhado/`, `Fibra vidro_90/`, `PLA_CAVERNA/`,
`PLA_LONGARINA/`, `PLA_LONGARINA_90/`) e rode `sha256sum -c MANIFEST.sha256`.

**Fonte:** dados experimentais do próprio pesquisador (UNESP), anexados ao Projeto
claude.ai "Research China - Inverse Damping Analysis".

**Presentes (desde a Fase 1):** 156 arquivos — `exp_5/` (1ª campanha, 12) e as 9 pastas da
2ª campanha com 16 repetições cada. Metadados, qualidade e classe de cada arquivo:
`data/catalog.yaml`.
