# Handoff: reviewer → claude (fase 3: segunda revisão curta de B1-a e B1-b)

Gravado por claude a partir do relato do reviewer, que é somente-leitura. **Veredito: APROVADO. B1-a e B1-b resolvidos; nada bloqueante encontrado.**

## Entregue

- **Escopo:** só B1-a e B1-b, sobre 234550d mais o diff sem commit. O reviewer não executou nada.
- **B1-a resolvido** (`phase3.py`, `_refine_s`):
  - a folha com item registrado numa ponta nova não passa como resolvida em nenhum dos três caminhos;
  - os itens das pontas originais já estão em `case["unresolved"]`;
  - o teste distingue o caso;
  - efeito medido: as 8 folhas em s que sobram têm todas o motivo de B1-a.
- **B1-b resolvido** (`phase3.py`, `unresolved_fraction` e o texto):
  - um item na grade conta com os dois intervalos vizinhos;
  - a cota está restrita ao critério de bracketing e à grade declarada, com o que ela não cobre dito;
  - a tolerância de "na grade" só pode aumentar a medida conservadora.

### Notas (nenhuma bloqueante)

- **`phase3.py:150`:** o corte `min(fine, conservadora)` não está na docstring.
- **`phase3.py:140-141`:** itens fora da grade além das pontas são latentes; uma asserção os tornaria explícitos.
- **`phase3.py:166-168`:** os itens em ν dos pontos novos não passam por `refine_nu`. Isso é conservador; passá-los poderia reduzir as 8 folhas que sobram.
- **N-c:** em `r8x`, com `after` e `refined=True`, a grade da fração seria a errada. Latente.
- **Cobertura:** o teste da fração não exercita `nu_state_lost`, `nu_bisect_lost` nem a ponta da grade.

## O que falta

1. Confirmar o pytest.
2. Commit com `run: --full`.
3. Tarefas para as notas. N-a continua recomendada antes do mapa.

## Riscos

- A cota de 2,5% só vale sob o critério de bracketing e na grade declarada.
- 21 dos 30 itens que sobram são raízes do meio que a busca principal não seguiu: são os candidatos concretos de Q1.
