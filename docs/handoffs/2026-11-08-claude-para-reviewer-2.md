# Handoff: claude → reviewer (fase 3: revisão curta de B1-a e B1-b)

Branch `feat/fase3-bilhar`, 234550d mais o diff sem commit. Origem: `docs/handoffs/2026-11-08-reviewer-para-claude.md`.

## Entregue

- **B1-a:** `chronosim/phase3.py`, `_refine_s`.
  - Cada ponto novo de s chama `middle_roots(..., log=)`.
  - Uma folha com item registrado numa ponta nova volta não resolvida, com o motivo.
  - Os itens nas pontas originais são refinados à parte, por `refine_nu`.
  - Teste: `tests/test_ekt_reproduction.py::test_refine_s_keeps_a_leaf_whose_new_end_logged_nu_items`.
- **B1-b:** `chronosim/phase3.py`, `unresolved_fraction`.
  - Um item em ν num ponto da grade conta, na medida conservadora, com os dois intervalos vizinhos.
  - A docstring e o texto do caso vazio dizem que a cota só vale sob o critério de bracketing e na grade declarada.
  - Teste: `test_unresolved_fraction_on_a_hand_made_grid`.
- **Notas aplicadas:**
  - N-b: `off_grid` é impresso.
  - N-d: a frase da regra de sinal só sai quando há soluções.
- **Saídas regeneradas** com `.venv/bin/python`: `--full` e depois a execução padrão. Caso vazio em `outputs/phase3_full.md:97-104`.

## O que falta

1. Revisar só B1-a e B1-b.
2. Classificar como bloqueante ou nota, e dar o veredito.
3. Resposta final: o conteúdo de `docs/handoffs/2026-11-08-reviewer-para-claude-2.md`.

## Riscos

- O critério de "resolvido" continua sendo o do bracketing.
