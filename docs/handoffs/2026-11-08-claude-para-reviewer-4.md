# Handoff · claude → reviewer (item B, certificação)

Worktree: `/home/fasn98/chronosim/.claude/worktrees/fase3-certificacao`, branch `feat/fase3-certificacao`. As mudanças estão sem commit; veja `git diff` e os arquivos novos `chronosim/certify.py`, `docs/phase3_certification.md`, `tests/test_certify.py` e `tests/test_certified_origin.py`. Somente leitura.

## Entregue

- **N-a:** registro dos descartes da busca da mira (`billiards.py`, `phase3.py`).
- **CERTIFIED:** em `provenance.py`, conforme `docs/phase3_plan.md`, segunda rodada, item 1.
- **Exclusão intervalar** do caso vazio (θ, φ) = (0; 0,4), B = 0,3, R/B = 0,45, para uma travessia e colisão antes da boca 2. Fiz duas exclusões: a dos 427 itens não resolvidos e a de Ω inteiro, que deu 656.220 subcaixas, todas excluídas.

## O que falta (para você)

1. **A prova é válida?**
   - As condições e desigualdades de `certify.conditions` são **necessárias** para uma solução do modelo de `billiards.py`?
   - A caixa cobre de fato o espaço declarado, com as pontas arredondadas para fora?
   - O critério de exclusão está correto?
   - Algum caminho do código exclui uma caixa sem prova (exceção engolida, divisão por intervalo que contém 0, `sqrt` de intervalo negativo)?
2. **A redação** do relatório e do documento promete mais do que foi provado?
3. **As regras de CERTIFIED** respeitam a decisão: "origem só piora", nada promovido, combinação com COMPUTED dá COMPUTED.
4. **N-a:** sobrou algum descarte silencioso?

## Riscos

- `CERT_EXCLUSION` é um rótulo proposto, à espera de aprovação do usuário.
- Krawczyk não foi implementado (não foi preciso no caso vazio).
