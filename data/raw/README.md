# data/raw — dados brutos (somente leitura)

Os CSVs **não** são versionados no git (≈2 GB). A integridade é garantida pelo
`MANIFEST.sha256` (gerado com `sha256sum`). Para reconstruir: copie as pastas originais
de `ultima versao/` para cá, **mantendo os nomes** (`exp_5/`, `AÇO/`, `Al/`, `ACRILICO/`,
`FIBRA_CARBONO/`, `Fibra vidro_ alinhado/`, `Fibra vidro_90/`, `PLA_CAVERNA/`,
`PLA_LONGARINA/`, `PLA_LONGARINA_90/`) e rode `sha256sum -c MANIFEST.sha256`.

**Fonte:** dados experimentais do próprio pesquisador (UNESP), anexados ao Projeto
claude.ai "Research China - Inverse Damping Analysis".

**Presentes nesta cópia (Fase 0):** 60 arquivos — `exp_5/` (1ª campanha, 12),
`AÇO/` (16), `Al/` (16), `ACRILICO/` (16). As demais pastas da 2ª campanha entram na Fase 1.
