# Handoff: claude → reviewer (fase 3, Echeverria, Klinkhammer & Thorne 1991: reprodução R1–R11)

Branch `feat/fase3-bilhar`, escopo `main..e7bd1d8`. A árvore de trabalho pode mudar durante a revisão (preparação da divisão `run` / `run --full`); revise o commit, com `git show e7bd1d8:<caminho>` ou `git diff main..e7bd1d8`.

## Entregue

- **Código:** `chronosim/billiards.py`, `chronosim/ekt.py`, `chronosim/phase3.py`, `chronosim/run.py`, `chronosim/sources.py`.
- **Testes:** `tests/test_ekt_transcription.py`, `tests/test_ekt_reproduction.py`.
- **Documentação:** `docs/equation_transcription.md` (parte da fase 3), `docs/phase3_plan.md`, `outputs/report.md`, `README.md`.
- **Governança:** regras do projeto em `CLAUDE.md`; decisões do usuário em `docs/phase3_plan.md`, seção "Decisões do usuário (2026-11-07)".

## O que falta (o que revisar)

1. **Honestidade científica:** toda saída é `Quantity` com origem correta. Nada do modelo geométrico pode estar acima de COMPUTED, e nada transcrito pode estar acima de DERIVED. Toda fórmula de `ekt.py` deve apontar para a equação exata. Constantes sem fonte devem passar por `heuristic()`.
2. **O texto do relatório** (`phase3.report_section`, `outputs/report.md`) diz só o que os números mostram. Nenhum "paradoxo" é declarado: o caso sem solução deve aparecer como "nenhuma solução encontrada sob [critério]". Achados sobre o artigo devem usar linguagem neutra ("inconsistência"); nunca "erro de impressão" ou "erro tipográfico".
3. **Correção do modelo:**
   - regra de travessia e itinerário (`billiards.traverse`, `_gap`, `itinerary_windows`);
   - verificador de história (`billiards.history`, `verify`): ele precisa rejeitar colisões extras;
   - filtro de raízes das eqs. 4.2–4.3 (eq. 4.6 e distâncias do Apêndice A).
4. **Robustez numérica:** guardas de bisseção, mistura float/mpmath, contexto `mp.dps` vazando entre tarefas ou processos.
5. **Anti-RNG:** nada alcança um gerador aleatório; as grades são determinísticas.
6. Classificar cada achado como **bloqueante** ou **nota**. Concluir "nada bloqueante" é aceitável. Devolver em `docs/handoffs/2026-11-07-reviewer-para-claude-3.md`.

## Riscos

- O termo em B das duas formas da eq. B20 difere por um fator 4. A geometria segue a segunda forma; a conferência humana dessa leitura está pendente.
- Minha eliminação para a eq. B30 dava um termo em R diferente. A geometria favorece a forma impressa, por isso não aponto inconsistência no artigo.
- A citação 12 (p. 1088) está pendente de conferência humana.
