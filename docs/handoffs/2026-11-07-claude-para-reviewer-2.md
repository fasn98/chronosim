# Handoff: claude → reviewer (fase 2, segunda rodada)

Branch `feat/fase2-kim-thorne`. Escopo: só os commits desde `ae1e7fe`, inclusive, ou seja, `git diff ae1e7fe~1 HEAD`.

## Entregue

- **`ae1e7fe`:** decisões do usuário sobre os achados da primeira rodada (`docs/handoffs/2026-11-07-reviewer-para-claude.md`):
  - I3: a referência de Unruh passou a HEURISTIC (`chronosim/unruh.py`).
  - I4: `ValueError` para entradas fora do domínio e para toda soma vazia (`chronosim/kim_thorne.py`, `chronosim/thermal.py`, `tests/test_input_validation.py`).
  - I6: as leituras de F7 ficaram marcadas como pendentes.
  - Notas: vírgula decimal, "Até aqui:" e o caso-limite de rótulo.
- **Commit seguinte:** registra a verificação humana das leituras de F7 (`docs/equation_transcription.md`, `chronosim/estimates.py`) e a redação do usuário sobre a eq. 14 (`chronosim/run.py`, `outputs/report.md`).
- **Estado:** `python -m pytest` com 409 testes passando, e `python -m chronosim.run` executado.

## O que falta

Revisão somente-leitura desse diff, em especial:

1. Algum caminho válido de `run.py` ou dos testes agora levanta `ValueError` sem querer?
2. A mudança de `_series_origin` mantém a regra de que a origem só piora?
3. O texto do relatório usa a redação pedida, sem a linguagem vetada?

## Riscos

- `thermal_stress` agora exige `n_terms` ≥ 2, o que é mais restritivo que o pedido do usuário. Está justificado no código.
