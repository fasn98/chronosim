# Handoff: reviewer → claude (fase 3: revisão curta de B1 e B2)

Gravado por claude a partir do relato do reviewer, que é somente-leitura. **Veredito: BLOQUEADO: 2 achados bloqueantes (B1-a, B1-b). B2 está resolvido.**

## Entregue

- **Escopo:** só B1 e B2 de `docs/handoffs/2026-11-07-reviewer-para-claude-3.md`, em 234550d mais o diff sem commit, e as saídas `outputs/phase3_full.md` e `outputs/report.md`. O reviewer não executou nada.
- **B1 conferido sem achado:**
  - `_nu_intervals`, `middle_roots` e `crossings` registram tudo o que não seguem;
  - `max_split = 0` é respeitado;
  - os intervalos de s estão declarados como heurísticas e impressos;
  - `refine_nu` devolve as folhas com motivo;
  - em `refine_unresolved`, as bisseções perdidas continuam contadas;
  - o texto do caso vazio depende dos dados.
- **B2 conferido sem achado:**
  - `sign_rule` compara o sinal de β com o sinal de d medido, e `dev_a6` = |a(A6) + s|;
  - os dois entram na tabela de R8, no resumo e na execução padrão;
  - a frase do árbitro depende dos dados;
  - "Estabelece" restringe o limite B → 0 a 3.1, 3.2 e 4.12–4.13.

### BLOQUEANTE B1-a: `phase3.py` (`_refine_s`)

- **O quê:** `_refine_s` chama `middle_roots` sem `log`. Nos pontos novos de s, os itens em ν que não puderam ser seguidos são descartados em silêncio, e a folha pode passar como resolvida justamente por faltar um ramo.
- **Correção:** passar `log` e tratar como não resolvida a folha que tiver item registrado numa ponta nova.

### BLOQUEANTE B1-b: `phase3.py` (texto do caso vazio e docstring de `unresolved_fraction`)

- **O quê:** a afirmação de que a medida conservadora "é a cota superior" não está garantida, por três motivos:
  - a célula dual de um ponto da grade é estreita demais, porque a ligação de v₁(s) é feita nos dois intervalos vizinhos;
  - a cota só vale sob o critério de bracketing e na grade declarada;
  - o descarte de B1-a também entra na medida.
- **Correção:** para itens na grade, usar a união dos dois intervalos vizinhos e qualificar a frase.

### Notas

- **N-a:** `roots_by_bracketing`, na busca da mira, ainda descarta sem registro os colchetes com f indefinida e os saltos rejeitados. As docstrings "Nothing is dropped silently" são fortes demais. O problema já existia antes.
- **N-b:** `off_grid` é calculado, mas nunca impresso.
- **N-c:** em `r8x` com `kind == "after"` e `refined=True`, a grade usada na fração não seria a grade da busca. Hoje é latente.
- **N-d:** sem nenhuma solução, a frase sai "vale em todas as 0 soluções".
- **N-e:** "Estabelece" é texto fixo, e o relatório não imprime a ordem medida para 4.13 (classe II).
- **N-f:** valeria calcular v₁ nas 21 raízes do meio que a busca principal não seguiu.

## O que falta

1. Corrigir B1-a e B1-b.
2. Rodar `pytest` e `run --full`.
3. Nova revisão curta, só de B1-a e B1-b.
4. N-a e N-b merecem tarefa antes do mapa.

## Riscos

- O reviewer não executou testes nem o run.
- B1-a é um descarte demonstrado no código; o efeito sobre o resultado não foi medido.
- O cenário de B1-b.1 é um contraexemplo possível, não observado.
