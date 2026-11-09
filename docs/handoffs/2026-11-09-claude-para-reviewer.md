# claude → reviewer · controle positivo e rótulo "raiz certificada" (2026-11-09)

## Entregue
Worktree `.claude/worktrees/fase3-certificacao`, diff ainda sem commit sobre `feat/fase3-certificacao`, nos arquivos:
- `chronosim/provenance.py`
- `chronosim/certify.py`
- `chronosim/phase3.py`
- `tests/test_certify.py`
- `tests/test_certified_origin.py`
- `docs/phase3_certification.md`
- `docs/phase3_plan.md`

Decisões: `docs/phase3_plan.md`, blocos de 2026-11-08 e 2026-11-09.

## O que falta
Revisar, somente leitura:
1. Os critérios `passed` e `complete_near_solutions` de `control_summary` são fiéis à decisão (b)? Uma solução ou cópia em caixa excluída tem que falhar o controle.
2. A mudança de escala (`scale`) do controle mantém "a mesma máquina" do caso vazio?
3. Nenhum caminho promove origem nem cria CERTIFIED sem o rótulo nomeando o sistema.
4. O texto do relatório não afirma nada além do que foi calculado. O texto das tabelas está em `_control_lines`.

## Riscos
- A escala do controle mudou de "relativa ao pedaço" para "1/8 do intervalo de s", como no caso vazio. Antes disso, os pedaços colados à boca da Fig. 7(b) explodiam.
