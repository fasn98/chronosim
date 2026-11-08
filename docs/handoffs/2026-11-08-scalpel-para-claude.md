# Handoff: scalpel → claude (fase 3: correções B1/B2 e divisão run/--full)

Branch `feat/fase3-bilhar`, árvore principal. Nada de `chronosim/` foi alterado, nada foi escrito em `outputs/`, nada foi commitado. A `--full` em andamento não foi tocada. Python usado: `.venv/bin/python` (o python do PATH não tem pytest).

## Resultado real do pytest

1. Suíte existente (inclui os 3 testes novos de `tests/test_ekt_reproduction.py`), `python -m pytest -q`:
   `439 passed in 472.91s (0:07:52)`. Rodou com a `--full` competindo por CPU; uma primeira tentativa com `timeout 580` foi cortada aos ~50% e não conta.
2. Arquivos novos, rodados depois (não estavam nos 439):
   - `tests/test_ekt_phase3_run.py`: `17 passed in 14.04s`.
   - `tests/test_ekt_gaps.py`: `3 failed, 37 passed in 37.08s` (as 3 falhas estão em "Bloqueante?" abaixo).

## Verificação de que os testes pegam B1 e B2

Feita em cópia fora do repositório (scratchpad), com uma mutação por vez em `billiards.py`/`phase3.py`:

| Mutação (reverte a correção) | Teste que falha |
| :- | :- |
| `nu_state_lost` não registrado (B1) | `test_middle_roots_logs_nu_state_lost_when_the_exchange_is_undefined_at_a_bracket_end` |
| `nu_count_change` não registrado (B1) | `test_nu_interval_with_a_count_change_is_logged_not_skipped` |
| `s_bisect_lost` não registrado (B1) | `test_crossings_reports_a_bisection_that_lost_its_branch_as_s_bisect_lost` |
| rejeição por impulso ≤ 0 não registrada (B1) | `test_crossings_discards_a_crossing_that_needs_a_negative_kick` |
| `dev_a6` removido das linhas de R8 (B2) | `test_r8_rows_carry_the_eq_A6_and_sign_rule_checks` (KeyError) |

Nota: `test_search_reports_every_discard_with_its_kind_or_reason` (já existente) **passou nas mutações de B1**; no caso testado a lista pode estar vazia ou os tipos serem legítimos, então ele não detecta a regressão sozinho. Os testes acima é que protegem B1.

## Bloqueante?

Nada nas mudanças sem commit está quebrado pelas correções. Três testes de `tests/test_ekt_gaps.py` estão vermelhos; deixei vermelhos, sem afrouxar. Os três exercitam código que o diff não toca (`roots_by_bracketing`, `crossings`/`max_split`), então são falhas **preexistentes** a este diff, achadas pela rodada anterior (marcadas "Measured RED" nas docstrings):

- **B-a** `test_roots_by_bracketing_does_not_report_a_zero_on_an_interior_node_twice`: para `2 - x` em `[1, 2, 3]`, devolve `[1.9999999999995453, 2.0]`. O zero num nó interior é contado duas vezes.
- **B-b** `test_roots_by_bracketing_finds_a_root_sitting_exactly_on_the_last_grid_node`: falha na mesma família (raiz exatamente no último nó da grade).
- **B-c** `test_crossings_with_max_split_zero_reports_the_whole_interval_without_halving`: `max_split = max_split or MAX_SPLIT.value` (`billiards.py`, em `crossings`) trata 0 como ausente; o intervalo é bisseccionado 12 vezes e relatado com 1/4096 da largura (obtido `('s_count_change', 0.69995..., 0.70019..., 1, 0)`).

Decisão para você: se são bloqueantes depende de a grade real de R1–R11 pôr um zero num nó interior (B-a/B-b) ou de alguém passar `max_split=0` (B-c). Não verifiquei se a grade real cai nesse caso; não afirmo impacto nos resultados da fase 3.

## Notas (não bloqueantes)

- **N-1** Lacuna: `run.main` com `--full` **não tem teste**. `main` roda as fases 1 e 2 inteiras e escreve em `outputs/`; sem poder tocar em `outputs/`, não executei. Só verifiquei por leitura (`run.py:870-932`): `--full` escolhe `report_section` e escreve `phase3_full.md`, CSV e fig. 5; sem ele, `report_section_quick`. Uma refatoração que separasse a escolha da seção em função pequena tornaria isso testável.
- **N-2** `report_section_quick` foi testado só com resultados montados à mão (formato, caso sem soluções, contagem de `sign_rule`, aviso de "não recalculados"). Não valida números científicos.
- **N-3** O risco do handoff (lista `log` compartilhada entre varredura e bisseções contando a mesma ocorrência mais de uma vez) **não** tem teste: não defini qual contagem seria a correta. Cobri só que o `log` passado a `crossings` é a mesma lista devolvida em `unresolved` e que entradas prévias são preservadas. Se a duplicação importa para o relatório, é decisão de especificação.
- **N-4** O caminho mpmath de `solutions` (k=12, B=0,01) é testado só por contrato (devolve par, tipos válidos, `mp.dps` restaurado), em grade minúscula. Não achei nele nenhuma solução verificada; não afirmo correção numérica ali.
- **N-5** `_circ` tem teste para mpf (fecha em 1e-35 com ε=1e-30; com `math.pi` falharia) e para float.
- **N-6** `test_collect_shuts_the_pool_down_even_when_a_task_fails` usa `ex.submit` após `collect` para checar o desligamento; depende do `RuntimeError` do `ProcessPoolExecutor` (comportamento de stdlib).
- **N-7** Não portei `tests/test_ekt_drivers.py` da worktree antiga: ~4–5 min extras, e os drivers R3–R11 já são cobertos em parte pelo `--full`. Está em `/home/fasn98/chronosim/.claude/worktrees/agent-a7568da6636b87b7c/tests/` se quiser.
- **N-8** A suíte inteira leva ~8 min com a `--full` rodando em paralelo; vale saber para a regra "pytest antes de cada commit".

## Arquivos de teste

- Criado: `tests/test_ekt_phase3_run.py` (17 testes).
- Criado (portado da worktree antiga e ajustado a `crossings` com 3 retornos, `middle_roots(log=)` e tuplas etiquetadas): `tests/test_ekt_gaps.py` (40 testes, 3 vermelhos).
- Não alterados: `tests/test_ekt_reproduction.py`, `tests/test_anti_fabrication.py`.
