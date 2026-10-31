# legacy/ — código original (somente leitura)

Cópia do código da monografia (TCC) e do relatório FAPESP, **sem nenhuma alteração de
lógica**. Serve como referência histórica e como alvo dos testes de caracterização da
Fase 0 (`tests/legacy/`, `notebooks/fase0/`). **Nunca edite estes arquivos.**

| Arquivo | Origem | SHA-256 (desta cópia) |
|---|---|---|
| `Biblioteca.py` | Projeto claude.ai, `ultima versao/Biblioteca.py` | `7dc816954be367691f89077d42a856cfa0ce1fad3e9263047b5d9779727294b1` |
| `main.py` | Projeto claude.ai, `ultima versao/main.py` | `a7c313a623329c89315ba2c1bc73ffafd7839214db62ab667b45af319513cae1` |
| `main_otimizador.py` | Projeto claude.ai, `ultima versao/main_otimizador.py` | `3e7073c1bbe45bda9f1878cd7473135bcce1cdf317657cfa4a1699ba0d32404f` |

**Nota de fidelidade.** Os arquivos foram transcritos a partir do texto lido no Projeto.
Diferenças em relação aos originais, todas só de espaço em branco: finais de linha
CRLF → LF, um caractere de tabulação no fim da linha `nu = 0.25` do acrílico removido e
uma quebra de linha final adicionada. Por isso o hash acima é o **desta cópia**, não o
do original. Para fechar a cadeia de custódia, substitua estes três arquivos pelos
originais do seu computador (o conteúdo de código é o mesmo) e atualize os hashes.

**Como o legado é executado na auditoria.** Via `src/dampid/audit/legacy_harness.py`:
os caminhos relativos (`exp_5/…`, `AÇO/…`) são resolvidos por links para `data/raw/`;
e, como o SymPy não pôde ser instalado no ambiente da auditoria, `matrixes()` é
substituída por uma reprodução numérica *linha a linha* das mesmas integrais,
verificada em `tests/legacy/test_legacy_element.py`. Todas as outras funções rodam
como estão.
