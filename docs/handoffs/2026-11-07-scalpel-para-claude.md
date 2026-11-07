# Handoff: scalpel para claude (fase 2, revisão de testes)

Branch `feat/fase2-kim-thorne`. Nada foi commitado. O pacote `chronosim/` e `tests/test_anti_fabrication.py` não foram tocados.

## Entregue

- **pytest, estado recebido:** `.venv/bin/python -m pytest`, 371 passed em 152,5 s, 0 falhas.
- **pytest, com meus testes:** 388 passed em 204,9 s, 0 falhas.
- **Arquivo novo, sem commit:** `tests/test_phase2_gaps.py` (17 testes): `series_terms` (limite exato, 9 casos), `unruh_temperature`, `metric_perturbation`, `hawking_cutoff_dt` como raiz da eq. 9c, soma sobre N da eq. 54 (um termo, soma explícita, escala 1/delta_T), eq. 68 numérica com ordens 1 a 3 contra a soma da eq. 54 (medido 5,3e-11, tolerância 1e-9), e a integral da garganta em delta_T = 1e-3 (medido 1,02e-3, tolerância 1,1e-3).
- **Medições** (rodadas por script em scratchpad, valores deterministas):
  - Os comentários "medido" conferem nos casos medidos: eq. 54 (0,0787 / 0,00907 / 0,00102 para delta_T 0,1 / 0,01 / 0,001), eq. 68 (5,5e-11), eq. 53 (1,25e-9, comentário diz 1,2e-9), eq. 69 (máximo 2,0e-8), Klinkhammer diagonal (máximo 4,0e-20).
  - Não medi o restante da suíte (F2, F3, F4, geometria, térmico, misner): lá só li os comentários.
  - Pytest não mostrou vazamento de `mp.mp.dps` (ver Riscos).

## O que falta (achados; nenhum é bloqueante)

1. `tests/test_kim_thorne_f6.py:43` (eq. 54, `dev[1] < 0.01`): medido 0,0091, folga de 10%. A tolerância não detecta erro de 0,1% no coeficiente, porque o desvio de truncamento (0,9%) é maior que o erro que se quer detectar. Mitigado pelo meu teste em `tests/test_phase2_gaps.py` (delta_T = 1e-3). Não bloqueante.
2. `tests/test_kim_thorne_f6.py:35` (`test_eq54_and_eq68_are_the_same_expression`): só verifica `< 0`. O nome promete identidade e o teste não a verifica. A soma sobre N da eq. 54 nunca era comparada a nada. Coberto agora em `tests/test_phase2_gaps.py` (testes da eq. 54 e eq. 68 ordens 1 a 3). Sugiro renomear ou apagar o teste antigo, decisão sua.
3. `tests/test_literature_factor2.py:62` (1e-6): medido 3,27e-8 (o comentário diz 6,5e-8 "medido em F3"; 6,5e-8 é o valor de eq4a contra 2 vezes Klinkhammer, e o desvio do with_half medido aqui é 3,3e-8). Folga de 30x. Como o alvo é distinguir o fator 2 (erro O(1)), a tolerância basta para o objetivo. Seria possível apertar para 1e-7. Não bloqueante. A identidade algébrica da linha 61 (1e-70, medido 1,7e-76) está certa.
4. `tests/test_literature_factor2.py:50`: só compara a diagonal de Klinkhammer com a do operador. O fora-da-diagonal medido é 0 em K, mas não há assert para isso. Lacuna pequena.
5. `tests/test_kim_thorne_f6.py:23` (eq. 69, 1e-7): único limite para todos os casos. N=1 mede entre 2e-16 e 3e-14; só N=3 chega a 2e-8. Para N=1 a tolerância tem 7 ordens de folga e deixa passar regressão pequena. Sugestão: tolerância por N.
6. `tests/test_kim_thorne_f6.py:53` (eq. 53): só o coeficiente do logaritmo é testado. O valor absoluto diverge 0,56% da fórmula em a = 1e-6 b (esperado, é termo dominante), mas isso não aparece em teste nenhum.
7. Sem validação de entrada (código do pacote, só relato): `kt.doppler_xi(1)` devolve 0, `doppler_xi(1.5)` devolve `mpc` (imaginário puro) e `doppler_xi(-0.5)` devolve 1,73 (xi > 1, série divergente). `kt.series_terms` com b·xi/2D >= 1 devolve inteiro negativo ou zero sem erro. `momentum_horizon(delta_T=0)` levanta ZeroDivisionError e `n_terms=0` devolve 0 sem aviso. Não há teste desses ramos porque não há comportamento especificado; decidir se viram `ValueError`.
8. Funções públicas sem teste direto, ainda: `pointsplit.image_pairs` (só indireto via `image_part`).
9. Os testes de análise dimensional (tolerância 1e-75 com dps 80) validam homogeneidade, não valores; não pegam erro numérico de constante.

## Riscos

- **Vazamento de precisão:** `tests/test_literature_factor2.py:44`, `:54` e `tests/test_kim_thorne_f4.py:140` atribuem `mp.mp.dps` direto. O autouse `tests/conftest.py:9` usa `mp.workdps`, que restaura a precisão ao sair, então o vazamento não ocorre, e a suíte inteira passou em ordem normal. Depende do fixture continuar autouse; se alguém o tornar opt-in, esses três testes contaminam os seguintes.
- **Tempo:** `test_kim_thorne_f6.py` leva cerca de 2 min; a suíte inteira com o meu arquivo, cerca de 3,4 min. Meu teste de delta_T = 1e-3 acrescenta cerca de 50 s e o de ordens 1 a 3 cerca de 20 s.
- **Tolerância ajustada após falha:** `tests/test_thermal.py::test_density_has_dimension_length_to_minus_4` (tolerância 1e-58, piso estimado em cerca de 2e-59, `assert tolerance > 3 * floor`). A justificativa a priori é coerente; não a reproduzi, só li o argumento.
- **O que esta revisão não estabelece:** não rodei mutação no pacote (proibido alterar código). A afirmação "tolerância frouxa" nos itens 1, 3 e 5 vem da comparação entre tolerância e valor medido, não de uma regressão injetada.
