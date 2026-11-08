# Handoff: claude → scalpel (fase 3, Echeverria, Klinkhammer & Thorne 1991: reprodução R1–R11)

Branch `feat/fase3-bilhar`, escopo `main..e7bd1d8`. Regras do projeto em `CLAUDE.md` e `README.md`. Plano aprovado em `docs/phase3_plan.md` (seção "Decisões do usuário"). Não há ADR: o portão foi a aprovação do plano pelo usuário.

## Entregue

- **Código:** `chronosim/billiards.py` (modelo geométrico independente, COMPUTED), `chronosim/ekt.py` (equações impressas, DERIVED), `chronosim/phase3.py` (R1–R11 e seção do relatório), `chronosim/run.py` (chamada da fase 3, figura 5), `chronosim/sources.py` (`EKT_1991`).
- **Testes:** `tests/test_ekt_transcription.py` (13), `tests/test_ekt_reproduction.py` (12).
- **Documentação:** `docs/equation_transcription.md` (parte da fase 3), `docs/phase3_plan.md`, `outputs/report.md` (seção da fase 3), `outputs/ekt_r8.csv`.
- **Estado:** 436 testes passando (`python -m pytest`, cerca de 3,8 min); `python -m chronosim.run` concluído (cerca de 18 min, 13 na fase 3).

## O que falta

1. Rodar `python -m pytest` e relatar a saída real.
2. Conferir, teste a teste, que cada tolerância tem o valor medido ao lado e que esse valor está correto: meça, não confie no comentário. Atenção a `tests/test_ekt_reproduction.py:127` (larguras das janelas de itinerário, cotas largas) e `:50` (razão d1/d2 com folga 0,1).
3. Procurar tolerâncias frouxas demais para detectar uma regressão e lacunas de cobertura: funções públicas de `billiards.py`, `ekt.py` e `phase3.py` sem teste; ramos de erro sem teste (por exemplo `ekt.appendix_b_coefficients`, `ekt.b11_v1`, `phase3.r11_after_traversal`, `billiards.crossings` com `max_split`).
4. Verificar que nenhum teste da fase 3 altera `mp.mp.dps` globalmente e vaza para outros testes.
5. Não alterar código do pacote nem enfraquecer `tests/test_anti_fabrication.py`. Testes novos ficam sem commit; liste-os no handoff de volta (`docs/handoffs/2026-11-07-scalpel-para-claude-2.md`, na árvore principal `/home/fasn98/chronosim`).

## Riscos

- A busca de raízes por bracketing aninhado tem intervalos "não resolvidos" contados no relatório; um teste que exija contagem exata pode ser frágil se a grade mudar.
- `phase3.reproduce` é lento (≈13 min com processos paralelos); não o chame inteiro num teste.
- Testes com `mp.workdps(40)` dependem do backend escolhido por `billiards.backend_for`; o limiar é a heurística `ekt_float_depth`.
