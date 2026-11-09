# Handoff · claude → scalpel (item B, certificação)

Worktree: `/home/fasn98/chronosim/.claude/worktrees/fase3-certificacao`, branch `feat/fase3-certificacao`, a partir de `main` (0ed21ec). Mudanças ainda sem commit. Python: sempre `.venv/bin/python` da worktree.

## Entregue

- **N-a:** a busca da mira passa a registrar os colchetes descartados.
  - `chronosim/billiards.py`: `roots_by_bracketing` (parâmetros `log` e `where`), `_aim_roots` (tipos `aim_bracket_lost` e `aim_bracket_jump`), `_nu_intervals`, `middle_roots` e `refine_nu`.
  - `chronosim/phase3.py`: `AIM_KINDS`, `item_span` (fatorado de `unresolved_fraction`), `refine_unresolved` e os textos do relatório.
- **Origem CERTIFIED:** `chronosim/provenance.py`, com os rótulos, `certified()`, `physical_solution()` e as regras de `combine()` e `degrade()`. Testes em `tests/test_certified_origin.py`, escritos antes da implementação.
- **Exclusão intervalar:** `chronosim/certify.py`; em `chronosim/phase3.py`, `certify_items`, `certify_chunk`, as tarefas `r8r` e `r8c` e `_certification_lines`. Testes em `tests/test_certify.py`. O enunciado está em `docs/phase3_certification.md`.
- **Decisões que governam:** `docs/phase3_plan.md`, segunda rodada, item 1, e o pedido B do usuário: rótulos, testes antes do uso e a redação "provado, em aritmética intervalar…".

## O que falta (para você)

1. Testes que tentem quebrar a exclusão, por exemplo uma caixa que contém uma solução conhecida e mesmo assim aparece como excluída. Atenção à dependência intervalar (`x*x` contra `x**2`) e às pontas de ν, da mira e de s.
2. Conferir que `conditions` é o mesmo modelo de `billiards._loop_state`, inclusive sinais e fatores positivos.
3. Testes de regressão para N-a, `item_span` e `unresolved_fraction` (os números não podem mudar no caso vazio).
4. Rodar o pytest. O `run --full` está rodando em paralelo; não rode outro.

## Riscos

- O registro de N-a fez a contagem do caso vazio, depois do refinamento, subir de 30 para 427 itens (397 colchetes da mira).
- O `CERT_EXCLUSION` é um rótulo proposto, ainda sem aprovação.
- `combine()` com CERTIFIED e CITED dá COMPUTED: é uma escolha deste branch.
