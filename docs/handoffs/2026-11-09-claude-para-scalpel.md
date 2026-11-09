# claude → scalpel · controle positivo e rótulo "raiz certificada" (2026-11-09)

## Entregue
Worktree `.claude/worktrees/fase3-certificacao`, branch `feat/fase3-certificacao`. As mudanças ainda não têm commit.
- `chronosim/provenance.py`: `CERT_ROOT_PREFIX`, `CERT_ROOT_MODEL` e `CERT_LABELS` com quatro rótulos.
- `chronosim/certify.py`:
  - `default_ranges`, `periodic_copies`;
  - `exclude(..., scale=None)`;
  - o `watch` mantém sem normalizar um ponto que já esteja dentro da caixa.
- `chronosim/phase3.py`:
  - `control_piece`: escala do caso vazio, cópias periódicas vigiadas, Krawczyk só nas caixas examinadas, volume não examinado;
  - `control_summary`: por solução, CERTIFIED com `CERT_ROOT_MODEL`;
  - `_control_lines`: duas tabelas;
  - nova heurística `ekt_cert_control_sol_budget`.
- Testes novos em `tests/test_certify.py` e `tests/test_certified_origin.py`.
- Decisões que governam a tarefa: `docs/phase3_plan.md`, bloco "Decisões do usuário sobre o relatório do controle positivo (2026-11-09)".

## O que falta
1. Cobrir:
   - `exclude` com `scale`: escala inválida, e que `scale` igual aos lados da caixa não muda nada;
   - `periodic_copies` nas bordas: ν exatamente 0 e mira exatamente −π;
   - `control_summary` quando uma cópia está numa caixa não examinada;
   - a recusa de rótulos "raiz certificada" sem o nome do sistema.
2. Rodar `python -m pytest` inteiro com o `.venv` e relatar a saída.

## Riscos
- `control_summary` agrupa as cópias pela chave normalizada da solução (comparação de floats).
- A fração do volume examinado supõe que os pedaços cobrem a caixa sem sobreposição.
