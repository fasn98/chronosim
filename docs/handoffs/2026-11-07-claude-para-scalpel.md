# Handoff: claude → scalpel (fase 2, Kim & Thorne 1991)

Branch `feat/fase2-kim-thorne`. Regras do projeto em `CLAUDE.md` e `README.md`. Plano aprovado em `docs/phase2_plan.md`. Não há ADR para esta fase: o portão foi a aprovação do plano pelo usuário.

## Entregue

- **Código:** `chronosim/kim_thorne.py` (fórmulas transcritas, DERIVED), `chronosim/wormhole.py` (numérica, COMPUTED), `chronosim/thermal.py` (árbitro térmico), `chronosim/unruh.py` (árbitro secundário), `chronosim/estimates.py` (F7), `chronosim/sources.py`, `chronosim/run.py`.
- **Testes:** `tests/test_kim_thorne_arbiter.py`, `tests/test_kim_thorne_geometry.py` (F1), `tests/test_kim_thorne_f3.py`, `tests/test_kim_thorne_f4.py`, `tests/test_kim_thorne_f2_f5.py`, `tests/test_kim_thorne_f6.py`, `tests/test_estimates.py`, `tests/test_thermal.py`, `tests/test_unruh.py`, `tests/test_literature_factor2.py`, `tests/test_sources.py`.
- **Documentação:** `docs/equation_transcription.md`, `docs/thermal_images.md`, `docs/literature_check_factor2.md`, `outputs/report.md`.
- **Estado:** a última execução teve 371 testes passando (`python -m pytest`, cerca de 2,3 min) e `python -m chronosim.run` concluído (cerca de 4,5 min).

## O que falta

1. Rodar `python -m pytest` e relatar a saída real.
2. Conferir, teste a teste, que cada tolerância tem ao lado o valor medido e que esse valor está correto: rode medições, não confie no comentário.
3. Procurar tolerâncias frouxas demais para detectar uma regressão. Exemplos a olhar: `tests/test_kim_thorne_f6.py` (eq. 54 com `dev[1] < 0.01`), `tests/test_literature_factor2.py` (1e-6 no segundo teste) e `tests/test_kim_thorne_f6.py::test_eq54_and_eq68_are_the_same_expression`, que só testa o sinal.
4. Procurar lacunas de cobertura: funções públicas sem teste, ramos de erro sem teste.
5. Não alterar código do pacote nem enfraquecer `tests/test_anti_fabrication.py`. Se escrever testes novos, deixe-os sem commit e liste-os no handoff de volta.

## Riscos

- `tests/test_thermal.py::test_density_has_dimension_length_to_minus_4` é o único teste cuja tolerância foi ajustada depois de falhar; ele traz uma estimativa a priori do piso de erro.
- `tests/test_kim_thorne_f6.py` leva cerca de 2 minutos (quadraturas).
- Alguns testes definem `mp.mp.dps = 80` ou `50` globalmente (`test_literature_factor2.py`, `test_kim_thorne_f4.py::test_eq52_chain_rule_identity`); verifique se isso vaza para outros testes.
