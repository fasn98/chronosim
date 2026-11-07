# Handoff: reviewer → claude (fase 2, segunda rodada)

Gravado por claude a partir do relato do reviewer, que é somente-leitura. Veredito: **APROVADO COM RESSALVAS**. Não há achado bloqueante nem importante; há cinco notas.

## Entregue

- **Escopo:** revisão somente-leitura de `git diff ae1e7fe~1 HEAD` na branch `feat/fase2-kim-thorne`, conferida contra o `CLAUDE.md` do projeto. Cobre `chronosim/kim_thorne.py`, `chronosim/thermal.py`, `chronosim/unruh.py`, `chronosim/estimates.py`, `chronosim/run.py`, `outputs/report.md`, `docs/equation_transcription.md`, `tests/test_input_validation.py` e `README.md`.
- **Pergunta 1:** nenhum caminho válido de `run.py` ou dos testes levanta `ValueError`. Os valores de β são sempre 0,1, 0,5 ou 0,9, `thermal_stress` recebe o padrão ou 12 termos, e `momentum_horizon` recebe de 1 a 4 termos.
- **Pergunta 2:** `_series_origin` (`chronosim/kim_thorne.py:94`) só rebaixa a origem, e nenhum caminho foi promovido.
- **Pergunta 3:** o texto da eq. 14 é o mesmo em `outputs/report.md:250` e em `docs/equation_transcription.md:93`, sem a linguagem vetada. O reviewer não tinha o original do usuário para conferir palavra por palavra.

## O que falta

1. **Opcional:** alinhar a regra de origem de `momentum_horizon` (`chronosim/kim_thorne.py:310`, DERIVED com n_terms = 1) com a de `stress_between_kk`, ou documentar a diferença na interface.
2. **Opcional:** aplicar `_br` também aos desvios em `chronosim/run.py:663`, para que a linha 159 de `outputs/report.md` não misture vírgula e ponto decimal.
3. **Opcional:** reescrever a frase com dois "e" em `chronosim/run.py:779-781`.
4. **Decisão do usuário:** manter `PLANCK_LENGTH_CM` (`chronosim/estimates.py:18`) como CITED com a ressalva atual, ou rebaixar a origem.
5. **Opcional:** a fonte de `chronosim/unruh.py:43` está em português e as demais em inglês. Se mudar, ajustar o teste em `tests/test_input_validation.py:61`.
6. **Antes do merge:** rodar `python -m pytest` e `python -m chronosim.run` de novo se houver mais commits.

## Riscos

- O reviewer não executou os testes, por causa do hook somente-leitura.
- `thermal_stress` exige `n_terms` ≥ 2 (`chronosim/thermal.py:75`). Isso é mais restritivo que o pedido do usuário, mas está justificado no código, e nenhuma chamada é afetada.
- Na assimetria de `momentum_horizon`, um `series_terms(...)` igual a 1 resultaria em DERIVED. O código atual não faz essa chamada.
