# reviewer → claude · controle positivo e rótulo "raiz certificada" (2026-11-09)

## Entregue
Revisão somente leitura do diff sem commit em `.claude/worktrees/fase3-certificacao` (branch `feat/fase3-certificacao`). Arquivos lidos: `chronosim/provenance.py`, `chronosim/certify.py`, `chronosim/phase3.py`, `tests/test_certify.py`, `tests/test_certified_origin.py`, `docs/phase3_certification.md`, `docs/phase3_plan.md`, `outputs/heuristics.csv`. Decisões aplicadas: `docs/phase3_plan.md`, blocos de 2026-11-08 e 2026-11-09.

Respostas:
1. `passed` falha com uma solução ou cópia em caixa excluída (`phase3.py:408`).
2. O `scale` reproduz a resolução do caso vazio (`certify_chunk`, `phase3.py:272-282`).
3. Não há caminho de promoção, e "raiz certificada" sem sistema é recusado em `Quantity.__post_init__`.
4. Há um texto falso no cabeçalho (`phase3.py:1050`) e um rótulo de coluna impreciso (`phase3.py:963`).

Nada bloqueante. Achados IMPORTANTES:
- `phase3.py:401/408`: um caso de controle sem soluções conhecidas passa sem testar nada.
- `phase3.py:401`: uma cópia em caixa nenhuma fica escondida por outra cópia não examinada.
- `phase3.py:394`: o source do CERTIFIED não diz a caixa nem as hipóteses, e o valor é o ponto COMPUTED.
- `phase3.py:1050`: o cabeçalho "Todas COMPUTED, exceto a exclusão" está falso com o controle.
- `phase3.py:963` e `:377-378`: "outras caixas" conta também as caixas com solução, e "none" numa caixa com solução não é sinalizado.

NOTAS:
- `phase3.py:391/394`: 1e-6 não registrado como heurística.
- `phase3.py:989`: "duas vezes" (até 4), e os intervalos estão fixos no texto.
- `phase3.py:985`: "caixas" conta cópias.
- `phase3.py:289`: constantes da Fig. 7(b) duplicadas.
- `certify.py:202-204`: a vigilância usa o ponto aproximado, não a raiz.

## O que falta
1. Corrigir os 5 achados IMPORTANTES.
2. pytest, que eu não pude rodar.
3. Regerar `outputs/`.
4. Rodar `--full`.
5. Registrar `run: --full` no commit.

## Riscos
- O Krawczyk só foi verificado na leitura.
- O encadeamento `start()`/`collect()` só será testado de verdade no `--full`.
- R8 #0, #3 e #4 podem ficar incompletos por orçamento, e então a frase "Provado" fica bloqueada.

APROVADO COM RESSALVAS
