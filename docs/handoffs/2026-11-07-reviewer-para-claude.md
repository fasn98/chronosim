# Handoff: reviewer para claude (fase 2, revisão de código, proveniência e segurança)

Gravado por claude a partir do relato do reviewer, que é somente-leitura. Veredito do reviewer: **BLOQUEADO por 1 achado (B1), de rótulo de proveniência, sem mudança de valores.** Nada mais é bloqueante.

Branch `feat/fase2-kim-thorne`, comparada com `main` por `git diff main...HEAD`, até o commit 30be426. A revisão foi só de leitura: o reviewer não rodou pytest nem `python -m chronosim.run`, e os números citados vêm de `outputs/report.md` e `outputs/kim_thorne_f3.csv`, commitados. Não existe `docs/adr/`; o documento de governo é `docs/phase2_plan.md`, aprovado.

## Entregue

1. **Proveniência:** 1 bloqueante (B1) e 2 importantes (I1, I3).
2. **Fórmulas contra `docs/equation_transcription.md`:** conferidas uma por uma as eqs. 3, 17, 29 a 33, 36 a 38, 42, 43, 44a, 46 a 48, 50 a 54, 57, 61 a 64, 66 a 69 (`chronosim/kim_thorne.py`) e 9c, 11, 13 e 14 (`chronosim/estimates.py`). Todas batem. O zero do denominador da eq. 67 é a eq. 62, conferido pela álgebra.
3. **Independência:** sem achado. `chronosim/wormhole.py:28-55` é escrito à parte de `pointsplit._tensor_at`, com separação em y. F3 deriva a eq. 66 sem usar a 67, F4 deriva a eq. 50 sem usar a 51, F2 usa matrizes ABCD, F5 acha raízes das eqs. 46 e 61, F6 integra as eqs. 51 e 67 já validadas, e o árbitro térmico usa a soma KMS contra π²T⁴/30. O lapso de F2 (`wormhole.py:191-198`) usa as eqs. 29 e 32, que são definições de coordenadas e não a eq. 57 sob teste.
4. **Afirmações do relatório:** a linguagem vetada não aparece em `chronosim/`, `docs/`, `tests/`, `outputs/` nem no `README.md`. Há uma contradição interna (I2) e uma frase forte demais (I5).
5. **Guarda anti-fabricação:** sem achado. `tests/test_anti_fabrication.py`, `tests/conftest.py`, `chronosim/provenance.py`, `pointsplit.py`, `misner.py` e `config.py` estão iguais a `main`, e nenhum módulo novo alcança `random`, `secrets`, `uuid`, `numpy.random`, `getattr`, `importlib`, `eval`, `exec` ou `subprocess`.
6. **Achado 7 do scalpel:** importante, não bloqueante (I4).

**Segurança:** nada encontrado. Não há rede, shell, eval nem segredo. O único caminho de escrita é `--out` da CLI local (`chronosim/run.py:806-811`). Nenhum PDF está no índice nem no histórico do branch, `.gitignore` exclui `papers/`, e `requirements.txt` não mudou.

### Bloqueante

- **B1. Somas truncadas rotuladas DERIVED** (`chronosim/kim_thorne.py:116`, `:184`, `:200`, `:270`). `stress_between_kk` (eq. 67), `stress_throat` (eq. 51), `klinkhammer_kk_on_axis` e `momentum_horizon` (eqs. 54 e 68) somam uma série até `n_terms`, que vem de `kt_series_tol` (HEURISTIC), e rotulavam DERIVED. Pela definição de `chronosim/provenance.py:8-10` e do `README.md`, isso é COMPUTED. Com `orders=[N]`, um termo só, DERIVED está correto. A diferença numérica é de cerca de 1e-30 relativo. Propagação: `chronosim/run.py:407`, `:644`, `:683` e `outputs/report.md`.

### Importante

- **I1.** `chronosim/estimates.py:336-340`: `hawking_cutoff_dt` usava `findroot` e rotulava DERIVED, embora a forma fechada seja Δt = l_P²/D.
- **I2.** `chronosim/run.py:788-790` e `chronosim/kim_thorne.py:138`: diziam que a eq. 52 "aguarda conferência humana", em contradição com `run.py:709` e `docs/equation_transcription.md:81`.
- **I3.** `chronosim/unruh.py:39`: a referência de Unruh é DERIVED, com fonte não lida. `unruh.py:18-25` não cita fonte. Atenuado por ser declarado e por ser só o árbitro secundário. **Decisão do usuário.**
- **I4.** `chronosim/kim_thorne.py:24-26`, `:75-78`, `:265-270`: sem validação de entrada. `series_terms` devolve um inteiro ≤ 0 quando b·ξ/2D > 1, e as somas dão 0 em silêncio. `doppler_xi` só faz sentido para 0 < β < 1. Nenhum caminho de `run.py` sai do domínio. **Decisão do usuário:** transformar em `ValueError`?
- **I5.** `docs/literature_check_factor2.md:47` e `:5`: as frases "ninguém apontou" e "confirma de forma independente" excedem a busca, que deixou 44 trabalhos sem leitura e teve transcrições de leitura única.
- **I6.** `chronosim/estimates.py:295-304`: os números impressos de F7 não estão em `docs/equation_transcription.md`. A discrepância entre 10⁴⁹ e 10⁵⁹ depende de uma leitura do scan. **Recomendação:** transcrição dupla ou leitura humana antes de publicar.
- **I7.** `chronosim/thermal.py:40` (`workdps(30)`), `:43` (31 subintervalos) e `chronosim/run.py:43` (`MODE_N = 400`) estavam fora de `heuristic()`.

### Notas

- `chronosim/run.py:359` e `:361` reimplementavam inline as eqs. 68 e 54.
- `chronosim/wormhole.py:42-55`: faltava em `eq4_stress` a checagem `h < delta`.
- `chronosim/run.py:656-657` e `:695` tinham parênteses fixos, não calculados dos dados.
- `chronosim/run.py:754-755`: a independência de ρ foi testada só em ρ/X ∈ {0, 1e-3, 1e-2}.
- Textos desatualizados: `chronosim/run.py:2`, `:521` e `:805`, `README.md:88` e `README.md:13`.

## O que falta

1. Corrigir B1 e ajustar legendas e texto.
2. Corrigir I2.
3. Decisões do usuário: I3, I4 e I6.
4. Opcional: I1, I5, I7 e as notas.
5. Rodar `python -m pytest` e `python -m chronosim.run`, com `outputs/` regenerado.
6. Uma rodada curta de revisão do diff das correções.

## Riscos

- O reviewer não rodou testes nem o pipeline e não leu os PDFs. Conferiu o código contra `docs/equation_transcription.md`.
- A guarda anti-fabricação é sintática (AST).
- B1 é de rótulo, não de valor.
- F6 roda com `kt_quad_dps = 15`.
