# claude → reviewer · 2026-11-10

## Entregue
Branch `feat/fase3-faixa-boca`, worktree `.claude/worktrees/fase3-faixa-boca`, mudanças não commitadas (sobre `de09722`).
- Decisão do usuário de 2026-11-10: as raízes de controle da faixa recebem o rótulo "raiz certificada das equações da máquina por direção" (`chronosim/provenance.py:46-55`, `CERT_ROOT_DIR`), usado em `chronosim/phase3.py` `_band_lines`; frase explícita de que o controle da faixa testa a máquina contra raízes das próprias equações.
- Testes: `tests/test_certified_origin.py`, `tests/test_certify_scalpel_2.py`, `tests/test_certify_band.py`.
- Análise do item (i) (EKT p. 1081) em `docs/phase3_borda.md`, seção 4; nota em `docs/phase3_certification.md:11`.

## O que falta
- Revisar o diff (`git diff`): rótulo, conjunto CERT_LABELS, testes, e a seção 4 (números conferidos contra o código: `billiards.verify` e `certify.conditions_dir`).
- `--certify` rodando; pytest 930 passed, execução padrão 394 s.

## Riscos
- O rótulo do controle do intervalo principal (CERT_ROOT_MODEL) não mudou: aquelas raízes passam no verificador.
