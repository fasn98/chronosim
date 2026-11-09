# scalpel → claude · controle positivo, casos de borda (2026-11-09)

## Entregue
- Arquivo novo: `.claude/worktrees/fase3-certificacao/tests/test_certify_scalpel_2.py` (55 testes). Nenhum código de produção nem teste existente foi alterado.
- Cobertura:
  - `exclude(scale=)`: escalas não positivas (5 variantes), escala com 2 e com 4 larguras, escala igual aos lados da caixa (caixa pequena e caixa padrão) sem mudança, escala muito maior pára a subdivisão cedo.
  - `periodic_copies`: ν = 0 (cópia em 2π), mira = −π (cópia em +π), ±π e 0/2π equivalentes (com tolerância de arredondamento), 4 cópias no canto, s intocado, ponto muito fora do círculo, faixa customizada sem a cópia.
  - `control_summary`: cópia em caixa não examinada (sem certificação, sem quantidade, não é falha, não passa para "refinado"), mistura de cópia certificada e não examinada, cópia em nenhuma caixa, solução sem entrada de vigia, vereditos não-raiz, caixa de raiz mais larga que 1e-6, rótulo e valor da quantidade.
  - Rótulos "raiz certificada": 10 grafias sem o nome do sistema ou erradas recusadas por `certified` e por `Quantity`; rótulo de raiz em origem não CERTIFIED; só dois rótulos com o prefixo; `degrade` perde o rótulo.
- Suíte inteira (`.venv/bin/python -m pytest -q`): **2 failed, 837 passed in 957.50s**. Antes do arquivo novo, a suíte inteira passou sem falhas (rodada com `-x`, saída só pontos, sem falhas). As 2 falhas são do arquivo novo e são legítimas (abaixo).

## Defeito encontrado
`exclude` não valida o número de larguras em `scale` (`chronosim/certify.py`, bloco de `widths`, perto da linha 190; o erro estoura em `certify.py:207`):
- `scale` com 2 larguras: `IndexError: list index out of range` no meio da busca, em vez de `ValueError` claro. Teste vermelho: `test_exclude_refuses_a_scale_that_does_not_have_three_widths`.
- `scale` com 4 larguras: aceita em silêncio e ignora a quarta. Teste vermelho: `test_exclude_refuses_a_scale_with_an_extra_width`.
- Impacto prático baixo: `control_piece` sempre passa 3. Correção sugerida: `len(widths) != 3` → `ValueError`. Não corrigi, por ordem do handoff.

## O que falta
- Corrigir a validação acima (então os 2 testes ficam verdes), ou me dizer que o comportamento é aceito e eu reescrevo os asserts.
- Revisão pelo reviewer.

## Riscos
- Float na chave de `control_summary`: agrupa pela chave normalizada em float. Funciona porque `known` e `watch["sol"]` vêm do mesmo ponto; se uma entrada vier deslocada de 2π num lado só, `(1+2π) % 2π` difere de `1.0` na última casa e a solução fica com 0 cópias (mostrei no rascunho; o teste final usa o mesmo ponto nos dois lados, o que é o uso real). Não é falha hoje.
- ν = 2π (float) normaliza para −2.4e-16, não 0. Fica dentro da faixa (−0.01) e as cópias continuam certas; testes usam tolerância 1e-12.
- Não testado: fração de volume com pedaços sobrepostos (suposição do handoff), `control_piece` ponta a ponta (cálculo pesado; a suíte existente em `test_certify.py` cobre partes), e o conteúdo de `_control_lines` além do já existente.
- `ctrl3.py` continuou rodando; a suíte levou 16 min por competir por CPU com ele. Um teste de ~3 s (`larger_scale`) usa escala 1e10 para parar na primeira caixa, de propósito, para não depender de orçamento.
