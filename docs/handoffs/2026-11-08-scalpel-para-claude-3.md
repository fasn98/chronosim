# Handoff · scalpel → claude (item B, certificação)

Worktree: `/home/fasn98/chronosim/.claude/worktrees/fase3-certificacao`, branch `feat/fase3-certificacao`, mudanças sem commit. Nenhum arquivo sob `chronosim/` foi editado por mim; `outputs/` intocado; `chronosim.run` não foi executado.

## Entregue

- Arquivo novo: `tests/test_certify_scalpel.py` (56 funções de teste, 146 casos com parametrização).
- Suíte completa (`.venv/bin/python -m pytest`, com o `run --full` rodando em paralelo): **759 passed, 0 failed**, 315 s.
- Bloqueador na exclusão intervalar: **nenhum encontrado**. Nenhuma caixa que contém uma solução verificada foi excluída, em nenhum teste.

Cobertura, por risco:

1. **Falsa exclusão** (`tests/test_certify_scalpel.py`, testes `test_a_box_containing_a_verified_solution_is_never_excluded`, `test_budget_exhaustion_keeps_the_solution_in_the_leftover`, `test_default_ranges_with_a_narrow_s_window_keep_every_solution`).
   - 7 soluções verificadas por `solve_mp` e `verify` dentro do próprio teste: 3 do caso central, 1 de R8_CASES[1], 1 de [3] (com ν − 2π = −0,0098, dentro da folga NU_PAD), 2 de [4].
   - A solução é colocada no centro, a 1e-7 de cada uma das 6 faces e num canto; o invariante é "alguma caixa indecidida contém a solução", com orçamentos de 1 a 1200 caixas.
   - Com as faixas padrão de ν e da mira, s estreito, e todas as representações 2π-equivalentes dentro da faixa.
2. **Pontas e partição** (`test_splitting_tiles_the_root_box`, `test_a_point_on_a_corner_or_on_a_split_plane_survives_a_pointwise_exact_test`, `test_complete_means_no_box_left_even_when_the_budget_ends_on_the_last_box`): a bissecção ladrilha a caixa raiz sem perder volume; um teste pontual exato no lugar do teste intervalar mantém pontos de canto e de plano de corte; `complete` com orçamento 1.
3. **Dependência intervalar e exceções** (`certify.py:74`, `certify.py:104`): invólucro dos 6 valores de `conditions` contém o valor de ponto em mpmath, em 480 caixas aleatórias (largura de 1e-7 a 2,5), 4 casos; caixa em que X − P contém 0 (L0 com divisão por intervalo que contém zero) não levanta exceção e não gera NaN; `ivprec` restaura a precisão após exceção e `exclude` não altera `iv.prec` global.
4. **Entradas** (`certify.py:61`): B, R, v1, p0 e e são invólucros dos decimais exatos, comparados com mpmath calculado de forma independente de `case_inputs` (comparar com a própria função deixaria a mutação por `float()` passar).
5. **Arredondamento externo de s** (`certify.py:156`): estrito para extremos decimais inexatos, iguais e mpf; `exclusion_quantity` só certifica se `complete`.
6. **Cola em phase3** (`phase3.py:248`, `phase3.py:272`): os pedaços de `certify_chunk` se encostam sem lacuna nem sobreposição, de n = 1, 3, 8, e cobrem o intervalo declarado; `certify_items` usa um `s` externo ao `item_span` de cada item, compartilha execuções de spans iguais e conta itens sem span como indecididos.
7. **Regressão N-a / `item_span` / `unresolved_fraction`** (`phase3.py:103`, `phase3.py:125`): comparação diferencial com uma cópia da função de `main` (0ed21ec) em 8 sementes aleatórias, tipos S e ν, grades de 2 a 15 pontos, no grid e fora dele: resultado idêntico, e lista vazia dá zeros. Também: tipo desconhecido, extremos da grade, grade de 1 ponto, área em ν só para tipos ν, item da mira igual ao item de ν na medida conservadora, e o log de `roots_by_bracketing` não muda as raízes.
8. **Regras de CERTIFIED** (`provenance.py`): `Quantity` direto com o rótulo físico e sem provas (no momento em que testei pela primeira vez, **não levantava**; a versão atual de `provenance.py`, que já valida no construtor, passa); provas como uma string única não destravam; provas em gerador; variações do nome da prova não contam; `degrade` de qualquer rótulo para COMPUTED/HEURISTIC solta o rótulo e não promove a CITED/DERIVED; nenhuma origem vira CERTIFIED por `degrade`; `combine` nunca devolve CERTIFIED para nenhuma mistura (3 rótulos × 4 origens × 1 a 3 repetições); `at_least=CERTIFIED` é recusado.

Validação por mutação (cópia do `certify.py` mutada e carregada por plugin do pytest em diretório temporário, original intocado; cada mutação derrubou pelo menos um teste novo, ou existente quando indicado):

- `dn.a>=0` → `.b>=0`: derrubado só por `test_default_ranges_with_a_narrow_s_window_keep_every_solution` (o mais lento, ~15 s); idem `along.b<=0` → `.a<=0`.
- `kick`, teste do zero por ponto médio, descartar uma metade na divisão, descartar a folha indecidida (`test_splitting_tiles_the_root_box`), descartar a pilha ao esgotar o orçamento, sinal em V, sinal em miss, sinal em C, refletir sem o fator 2, W = L0 − along/L0: derrubados pelo teste de solução central ou de ladrilhamento.
- AIM_PAD = 3,0 e NU_PAD = 0: `test_default_ranges_cover_a_full_turn_of_nu_and_of_aim`.
- `d[0]*d[0]` no lugar de `d[0]**2`: derrubado só por `test_the_box_count_adds_up`, por exceção do mpmath, não por uma asserção sobre o sinal; é uma detecção acidental.
- Sem `nextafter`: `test_outward_s_range_contains_the_exact_decimals_strictly`. `ivprec` sem restaurar: `test_ivprec_restores_after_an_exception`. `float(x)` nas entradas: `test_interval_inputs_enclose_the_decimal_strings` (e `tests/test_certify.py`).

## O que falta

1. As mutações em `phase3.py` (pedaços de `certify_chunk`, spans de `certify_items`) **não** foram executadas; a validação desses testes é só por raciocínio (um intervalo com lacuna ou sem arredondamento externo viola as asserções de igualdade e de desigualdade estrita).
2. A condição `conditions` foi conferida contra `billiards._loop_state` pelo teste já existente em `tests/test_certify.py` (fatores positivos L0 e W, em 2 casos); não encontrei divergência de sinal. Não há ainda um teste desse tipo para um segundo ponto de s próximo de 0 nem para k ≠ 1 (o modelo da exclusão só trata k = 1, e isso está escrito no enunciado).
3. Decisão sua: `CERT_EXCLUSION` continua um rótulo proposto, sem aprovação (como no seu handoff).

## Riscos

- `exclude` com uma faixa degenerada (extremos iguais em qualquer um dos três eixos) divide por zero em `certify.py:127` (largura relativa). Não é alcançável pelos chamadores atuais (`s_range_outward` sempre dá largura positiva), então não é bloqueador; fica como armadilha se alguém passar `nu_range` ou `aim_range` degenerado.
- `item_span` (`phase3.py:121`): um item fora da grade e depois do último nó dá um intervalo de largura zero. Os itens logados ficam dentro da grade, então não ocorre hoje; não tem teste que prove que isso não acontece.
- A exclusão é sobre o s ∈ [−S_FAR, −R·S_NEAR_OVER_R] e só para uma colisão e uma travessia (hipóteses em `certify.py:162`); soluções fora dessa família não são cobertas, e os testes não pretendem cobri-las.
- Os testes de solução são lentos (cerca de 2 s cada, ~60 s no total) porque dependem do iv do mpmath; `test_a_tiny_box_around_a_solution_ends_undecided_not_excluded` foi limitado a 4000 caixas (com 100 mil levava 115 s).
- Durante meu trabalho `certify.py` e `provenance.py` mudaram (AIM_PAD e NU_PAD viraram `heuristic()`, `Quantity` ganhou `proofs`, `exclusion_quantity` passou a dar valor 0 e a contar subcaixas na fonte). Os testes foram ajustados para essa versão e usam `.value` quando existe.
