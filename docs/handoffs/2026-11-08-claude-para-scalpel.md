# Handoff: claude → scalpel (fase 3: correções da revisão e divisão run/--full)

Branch `feat/fase3-bilhar`, árvore principal `/home/fasn98/chronosim`. Escopo: as mudanças **sem commit** sobre e7bd1d8 (`git diff`). Revisão de origem: `docs/handoffs/2026-11-07-reviewer-para-claude-3.md` (bloqueantes B1 e B2).

## Entregue

- **`chronosim/billiards.py`** (B1):
  - `_nu_intervals` subdivide intervalos de ν com mudança de contagem e registra o que sobra.
  - `middle_roots(..., log=)` registra o que antes era descartado em silêncio: `nu_count_change`, `nu_state_lost`, `nu_bisect_lost`.
  - `crossings` devolve `(cands, unresolved, rejected)`. `unresolved` vem marcado por tipo (`s_count_change`, `s_bisect_lost`); `rejected` traz o motivo (polo de v1(s), impulso ≤ 0).
  - `_circ` usa mp.pi com mpf (nota N3).
- **`chronosim/phase3.py`:**
  - `solutions` passa a respeitar os dígitos de `backend_for` (N3) e anexa os rejeitados com motivo.
  - Novos `unresolved_by_kind` e as heurísticas `ekt_s_far`, `ekt_s_near_over_R` e `ekt_s_near_after` (mesmos valores de antes).
  - Linhas de R8 com `dev_a6` e `sign_rule` (B2); `r8x` carrega o caso (N5).
  - Texto do relatório condicionado aos dados (B1, B2, N4).
  - Novos `QUICK_TASKS`, `full_tasks`, `start`, `collect` e `report_section_quick`.
- **`chronosim/run.py`:** flag `--full`. A fase 3 começa em processos antes das fases 1 e 2. Sem `--full`, a seção reduzida da fase 3 vai para `report.md`; com `--full`, também saem `outputs/phase3_full.md`, `ekt_r8.csv` e a fig. 5.
- **Testes:** `tests/test_ekt_reproduction.py`, 3 novos no fim (registro de descartes, `_nu_intervals` com `_aim_roots` falso, A6 e regra de sinal).
- **Documentação:** `CLAUDE.md` (regra run/--full, meta de ~6 min), `README.md`, `docs/phase3_plan.md`, `docs/equation_transcription.md` (título da seção, nota N8).
- **Material anterior:** uma rodada anterior sua, interrompida, deixou `tests/test_ekt_drivers.py` e `tests/test_ekt_gaps.py` na worktree `/home/fasn98/chronosim/.claude/worktrees/agent-a7568da6636b87b7c`, escritos contra e7bd1d8, sem handoff de volta. Aproveite o que valer, ajustado às novas assinaturas (`crossings` devolve 3 valores; `middle_roots` tem `log`).

## O que falta

1. Rodar `python -m pytest` e relatar a saída real. Uma `python -m chronosim.run --full` está rodando em paralelo na mesma árvore: não a interrompa e não escreva em `outputs/`.
2. Conferir que os novos testes de fato falhariam sem as correções B1 e B2, pelo menos um por bloqueante; por exemplo, revertendo mentalmente ou em cópia temporária fora do repositório.
3. Procurar lacunas nas funções novas: `start`, `collect`, `report_section_quick`, `--full` em `run.main`, `unresolved_by_kind` e o caminho mpmath de `solutions`.
4. Não alterar `chronosim/`, nem enfraquecer `tests/test_anti_fabrication.py`. Testes novos ficam em `tests/` da árvore principal, sem commit.
5. Handoff de volta: `docs/handoffs/2026-11-08-scalpel-para-claude.md`, separando bloqueante de nota.

## Riscos

- `phase3.reproduce` completo leva ~13 min; não o chame num teste.
- O argumento `log` de `middle_roots` é uma lista mutável compartilhada entre a varredura e as bisseções de `crossings`; a mesma ocorrência pode ser contada mais de uma vez.
