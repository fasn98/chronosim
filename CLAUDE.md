# chronosim · instruções para agentes

Leia o README.md antes de qualquer mudança. Princípio central: **honestidade científica estrutural**.

- Nunca invente resultados nem "arredonde" números. Relate o que a execução produziu, inclusive falhas.
- Toda saída científica é um `Quantity` com origem (`CITED`, `DERIVED`, `COMPUTED`, `HEURISTIC`). A origem só piora; nunca promova.
- Toda fórmula transcrita aponta para a fonte e a equação exata (`sources.py`). Se não tiver o artigo em mãos, pare e peça; não reconstrua fórmulas de memória.
- Constantes sem fonte são registradas com `heuristic()` e aparecem em toda execução.
- Nenhum módulo do pacote pode alcançar um gerador de números aleatórios (`tests/test_anti_fabrication.py`). Não enfraqueça essa guarda.
- Não "conserte" uma referência até ela concordar com o código: concordância obtida assim não é validação.
- Antes de cada commit: `python -m pytest` passando e `python -m chronosim.run` (execução padrão: fases 1 e 2 completas e um subconjunto pequeno da fase 3, meta de até ~6 min) executado, com `outputs/` atualizado.
- `python -m chronosim.run --full` (inclui R1–R11 completos da fase 3, escreve `outputs/phase3_full.md`, `outputs/ekt_r8.csv` e a figura 5) é obrigatória antes de qualquer merge em main e sempre que o código da fase 3 mudar (`billiards.py`, `ekt.py`, `phase3.py`).
- A mensagem de cada commit registra qual execução foi rodada: `run: padrão` ou `run: --full`.
- Cada fase declara o que estabelece e o que **não** estabelece.
