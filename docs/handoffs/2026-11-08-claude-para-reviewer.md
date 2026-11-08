# Handoff: claude → reviewer (fase 3: revisão curta dos bloqueantes B1 e B2)

Branch `feat/fase3-bilhar`, árvore principal `/home/fasn98/chronosim`. Pedido seu, de `docs/handoffs/2026-11-07-reviewer-para-claude-3.md`, "O que falta", item 4: uma revisão curta **só** de B1 e B2.

## Entregue

- **Commit 234550d:** as correções de B1 e B2 e a divisão run/--full. Ver `git show 234550d --stat`.
- **Sem commit, sobre 234550d** (`git diff`):
  - **Refinamento adaptativo dos intervalos não resolvidos do caso vazio:**
    - `billiards.refine_nu`; em `phase3.py`, `_refine_s`, `refine_unresolved` e a tarefa `r8r`;
    - a heurística `ekt_refine_max_split`.
  - **Fração do espaço de busca coberta pelos não resolvidos:** `phase3.unresolved_fraction`, com uma medida fina e uma conservadora. Entra no caso vazio e nas linhas do R7.
  - **Nova redação do caso vazio**, pedida pelo usuário: "nenhuma solução encontrada na região resolvida; N intervalos não resolvidos, cobrindo X% do espaço de busca".
  - **Versões dos pacotes:** `tests/test_report_format.py` confere as versões instaladas contra o `requirements.txt` e a formatação de números das fases 1 e 2. A causa da formatação diferente era rodar com o mpmath 1.3.0 do miniconda em vez do 1.4.1 do `.venv`.
- **Saídas:** `outputs/report.md` (execução padrão) e `outputs/phase3_full.md` (`--full`), regenerados com `.venv/bin/python`.

## O que falta (o que revisar)

1. **B1:**
   - Nada que o bracketing não consegue seguir é descartado em silêncio: `billiards.py` (`_nu_intervals`, `middle_roots`, `crossings`, `refine_nu`).
   - O texto do caso vazio corresponde ao que o código faz, com o critério e o espaço de busca declarados: `phase3.report_section`, o trecho do caso vazio, e `outputs/phase3_full.md`.
   - A definição das duas medidas de fração em `unresolved_fraction` é honesta. Em particular, confira se a "conservadora" é de fato uma cota superior da região em que uma solução poderia ter escapado. Os itens fora da grade estão contados.
2. **B2:**
   - As linhas de R8 trazem `dev_a6` e `sign_rule`, e o texto do árbitro está condicionado a eles.
   - "Estabelece" restringe o limite B → 0 a 3.1, 3.2 e 4.12–4.13.
3. **Veredito:** classificar como bloqueante ou nota. Gravo o seu relato em `docs/handoffs/2026-11-08-reviewer-para-claude.md`.

## Riscos

- O critério de "resolvido" no refinamento é o do bracketing: uma grandeza com um único sinal nas pontas. Isso não é prova. Ele pode errar se a grandeza trocar de sinal duas vezes dentro de uma folha.
- No caso vazio sobram 21 raízes do meio em ν que a busca principal não seguiu. O refinamento as detecta, mas não as incorpora à varredura em s.
