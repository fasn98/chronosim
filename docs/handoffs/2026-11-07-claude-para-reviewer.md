# Handoff: claude → reviewer (fase 2, Kim & Thorne 1991)

Branch `feat/fase2-kim-thorne`, contra `main`. Regras em `CLAUDE.md` (o do projeto e o global) e em `README.md`. Plano aprovado em `docs/phase2_plan.md`. A revisão de testes do scalpel está em `docs/handoffs/2026-11-07-scalpel-para-claude.md`: 388 testes passando, nenhum achado bloqueante.

## Entregue

- **Código:** `chronosim/kim_thorne.py`, `chronosim/wormhole.py`, `chronosim/thermal.py`, `chronosim/unruh.py`, `chronosim/estimates.py`, `chronosim/sources.py` e `chronosim/run.py` (diff contra `main`: `git diff main...HEAD -- chronosim/`).
- **Relatório e documentação:** `outputs/report.md`, `docs/equation_transcription.md`, `docs/thermal_images.md`, `docs/literature_check_factor2.md`.

## O que falta

Revisão somente-leitura, apontando arquivo:linha. Pontos a olhar em especial:

1. **Proveniência:** algum valor COMPUTED rotulado como DERIVED ou CITED? Alguma promoção de origem? Alguma constante sem fonte fora de `heuristic()`?
2. **Fórmulas transcritas** (`kim_thorne.py`, `estimates.py`) contra as transcrições de `docs/equation_transcription.md`.
3. **Independência:** os cálculos ditos independentes (`wormhole.py`, `thermal.py`) usam, sem querer, a fórmula que pretendem testar?
4. **Afirmações do relatório** (`run.py`, função `report`) que excedam o que os números mostram, ou que usem linguagem vetada pelo usuário: "erro de impressão", "erro tipográfico", "erro de Kim & Thorne".
5. **Guarda anti-fabricação:** nenhum módulo pode alcançar um gerador de números aleatórios.
6. **Achado 7 do scalpel:** falta de validação de entrada em `doppler_xi`, `series_terms` e `momentum_horizon`. Classifique a gravidade; a decisão fica com o usuário.

## Riscos

- `run.py` ficou longo, com relatório construído por concatenação de strings.
- As quadraturas de F6 usam 15 dígitos (`kt_quad_dps`), e todo o resto usa 80.
- `wormhole.eq4_stress` duplica de propósito a lógica de diferenças finitas de `pointsplit`, para medir o fator 2 em vez de impô-lo.
