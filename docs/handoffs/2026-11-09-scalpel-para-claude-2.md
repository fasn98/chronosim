# scalpel → claude · testes da faixa junto à boca, máquina por direção, --certify (2026-11-09)

## Entregue
- `tests/test_certify_band_scalpel.py` (novo, 64 testes, ~13 s). Nenhum arquivo de produção ou teste existente foi alterado.
- Resultado da suíte inteira com `.venv/bin/python -m pytest`: **924 passed in 1299.86s (0:21:39)**, 0 falhas. Esse número inclui os 64 novos. Os 15 de `tests/test_certify_band.py` já passavam antes.
- Cobertura, por item do handoff:
  - `conditions_dir` em iv: 12 caixas em torno da quina tangente (α = φ ± π/2) em s = 0 e ao lado, caixas que cruzam a quina, caixa degenerada em s = 0, grade 5³ de pontos da família contra o envoltório. Mais três verificações: finitude e largura limitada em s = 0, ausência de explosão em caixa pequena, intervalo fino em torno do valor pontual. Além disso rodei fora da suíte um fuzz de 74 524 pontos da família em caixas aleatórias de 3 casos, com 0 violações do envoltório.
  - `_meet` (`certify.py:161`): sobreposto, aninhado, tangente, único, semirreta, e disjunto (devolve o casco, não um intervalo invertido). `_nonneg` também.
  - `_reason_dir` (`certify.py:204`): um teste por motivo alcançável. Mais um teste de solidez: em 5 mil caixas da varredura, todo ponto da família amostrado viola a condição alegada.
  - `krawczyk(machine="dir")`: 5 caixas sem raiz dão `none`; caixa que contém a raiz nunca dá `none` (4 larguras, 4 deslocamentos); a máquina chega ao Jacobiano (K difere do da máquina pela mira); máquina desconhecida levanta `KeyError`.
  - `exclude(machine="dir")` com `watch`: ponto em caixa excluída; ponto fora da caixa-raiz não é normalizado; faixa padrão de α é `alpha_range`; orçamento de caixas esgotado dá `complete=False`; caso sem `pc` levanta erro.
  - `start(mode="certify")`: executor síncrono falso e `_task` stub. Confere a lista de tarefas, os `r8p` de seguimento do `r8pk`, que as fases 1–11 não rodam, que `quick`/padrão não rodam certificação, que `mode` explícito vence `quick`, e `KeyError` para modo desconhecido. A tarefa `r8rc` roda `certify_items` e a `r8r` não.
  - `run.main(["--certify"])`: simulado, porque o real leva ~100 min. `start` e `collect` foram trocados por um resultado enlatado. Escreve só `phase3_certify.md`, deixa `report.md` intacto, não chama `report`, imprime a frase de prova só quando completa, e lida com `res = {}`. `--certify` e `--full` são exclusivos nas duas ordens.

## O que falta
- Nenhum defeito de produção encontrado (arquivo:linha): nada a corrigir.
- Não coberto, com o motivo:
  - Motivos `dn ≥ 0` e `V ≠ 0` de `_reason_dir`: não alcançáveis em caixa real do caso vazio (varredura de ~8 mil caixas). Foram testados só em `_reason_from` (`certify.py:122`) com enlatados, não ponta a ponta.
  - `main()` real com os workers: não rodei, por causa do `--certify` em segundo plano e do custo. O fluxo real de `start` com `ProcessPoolExecutor` (pickling dos resultados das novas tarefas `r8b`/`r8bx`/`r8bp`) fica sem teste.
  - Caixas com s > 0: fora da família, ficam de fora de propósito (só o motivo "P dentro da boca" é testado).
  - Soundness do rigor fica por amostragem de pontos, não por prova.
- O teste de solidez e a quina tangente só valem para o caso vazio e os controles de `BAND_CONTROLS`; outros (θ, φ) não foram varridos.

## Riscos
- `test_every_reason_claimed_on_a_scan...` leva ~10 s; os demais são rápidos. A suíte inteira levou 21 min, mais lenta que o usual por competir com o `--certify` em segundo plano; este não foi tocado e nada foi escrito em `outputs/`.
- `_meet` com intervalos disjuntos devolve o casco, em silêncio. É conservador e seguro, mas esconderia um bug de envoltório em vez de acusá-lo. O teste fixa só que o resultado não é invertido e contém ambos.
- Os testes novos importam `tests/test_certify_band.py` (`_band_res`) por `sys.path`; se aquele helper mudar de nome, quebram.
