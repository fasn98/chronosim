# Handoff: claude → scalpel (artigo ReScience C, branch paper/rescience-kim-thorne)

Worktree: `/home/fasn98/chronosim/.claude/worktrees/paper-rescience`. Escopo: `git -C /home/fasn98/chronosim/.claude/worktrees/paper-rescience log main..HEAD` (ae66a33, 9424490, e61c2a4). Use sempre `/home/fasn98/chronosim/.venv/bin/python`.

## Entregue

- **`chronosim/run.py`:**
  - `write_phase2_tables` grava em CSV, com origem, as tabelas que só existiam no `report.md`;
  - `figures` e `figures_phase2` desenham cada figura em português e em inglês (`_en`) a partir dos mesmos dados.
- **Testes:** `tests/test_phase2_tables.py` e `tests/test_figures_lang.py`.
- **`paper/build_tables.py`:** gera `paper/generated/numbers.tex` e as tabelas a partir de `outputs/`, dos docs e do git. Tem a regra da tag `v0.2.0` e a contagem de trailers.
- **Estado:** 524 testes passando; execução padrão concluída.

## O que falta

1. Rodar o pytest e relatar a saída real.
2. Conferir que os testes novos falhariam diante de uma regressão. Exemplos:
   - um rótulo de CSV que diverge do relatório;
   - uma figura `_en` com dados diferentes da versão em português.
3. **Testes para `paper/build_tables.py`**, que hoje não tem nenhum:
   - `tex()` com números de 20 dígitos;
   - `literature_counts` e `transcribed_equations` quebrando quando o doc muda;
   - `release()` com a tag ausente, com HEAD diferente da tag e com a árvore suja (num repositório git temporário);
   - `trailers()`.
4. Não alterar `chronosim/` nem `paper/build_tables.py`. Testes novos ficam na worktree, sem commit.
5. Handoff de volta: `/home/fasn98/chronosim/docs/handoffs/2026-11-08-scalpel-para-claude-2.md`, separando bloqueante de nota.

## Riscos

- `build_tables.py` muda `mp.mp.dps` globalmente ao ser importado; num teste, isso pode vazar para outros testes.
- `tests/test_figures_lang.py` leva ~20 s.
