# chronosim · instruções para agentes

Leia o README.md antes de qualquer mudança. Princípio central: **honestidade científica estrutural**.

- Nunca invente resultados nem "arredonde" números. Relate o que a execução produziu, inclusive falhas.
- Toda saída científica é um `Quantity` com origem (`CITED`, `DERIVED`, `COMPUTED`, `HEURISTIC`). A origem só piora; nunca promova.
- Toda fórmula transcrita aponta para a fonte e a equação exata (`sources.py`). Se não tiver o artigo em mãos, pare e peça; não reconstrua fórmulas de memória.
- Constantes sem fonte são registradas com `heuristic()` e aparecem em toda execução.
- Nenhum módulo do pacote pode alcançar um gerador de números aleatórios (`tests/test_anti_fabrication.py`). Não enfraqueça essa guarda.
- Não "conserte" uma referência até ela concordar com o código: concordância obtida assim não é validação.
- Antes de cada commit: `python -m pytest` passando e `python -m chronosim.run` executado, com `outputs/` atualizado.
- Cada fase declara o que estabelece e o que **não** estabelece.
