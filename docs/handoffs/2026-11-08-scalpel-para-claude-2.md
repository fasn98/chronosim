# Handoff: scalpel → claude (artigo ReScience C, branch paper/rescience-kim-thorne)

Worktree: `/home/fasn98/chronosim/.claude/worktrees/paper-rescience`. Python: `/home/fasn98/chronosim/.venv/bin/python`.

## Entregue

- Novo, sem commit: `tests/test_build_tables.py` (49 testes). Nada em `chronosim/` nem em `paper/build_tables.py` foi alterado por mim.
- Resultado real do pytest:
  - Linha de base, suíte da branch antes dos meus testes: `524 passed in 262.97s (0:04:22)`.
  - Suíte completa, com 43 dos meus testes: `567 passed in 257.44s (0:04:17)`.
  - Depois acrescentei 6 testes (ver Nota 1) e rodei só `tests/test_build_tables.py`: `49 passed in 0.38s`. A suíte completa **não** foi repetida depois desses 6; o total esperado é 573, não medido.
- Nada vermelho, nenhum teste desabilitado ou afrouxado.

## Item 2: os testes novos falhariam diante de regressão?

Verificado em cópia descartável da worktree (em `scratchpad`, não na worktree). Cada mutação falhou pelo motivo certo:

| Mutação | Teste que falhou |
| :- | :- |
| Rótulo de `kim_thorne_f6.csv` divergente do `report.md` (0.1 → 0.2, 0.001 → 0.002) | `test_phase2_tables.py::test_f6_labels_match_the_report` |
| Figura `_en` fig2 com dados diferentes (×1,1 só em `en`) | `test_figures_lang.py::test_english_and_portuguese_figures_plot_the_same_values[fig2]` |
| `release()` sem a checagem de árvore suja | 2 testes de `release`, árvore modificada e staged |
| `release()` sem a comparação HEAD x tag | 2 testes de `release`, HEAD além da tag |
| `trailers()` sem exigir "Claude" | 4 testes de `trailers` |
| `trailers()` com critério de merge errado | 2 testes de `trailers` com merge |
| `transcribed_equations()` sem a contagem | 3 testes |
| `mp.mp.dps = 15` no script | 3 testes de `tex()` e precisão |

## O que `tests/test_build_tables.py` cobre

- `tex()` com 20 dígitos, 19 noves arredondando para 2, zero, potência de 10, sinal, vai-um do arredondamento e perda de dígitos em dps 15.
- `Macros`: nome só com letras e nome duplicado.
- `literature_counts` e `transcribed_equations` quebrando quando o doc muda. Cada falha por frase alterada, por aritmética inconsistente, por linha apagada e por linha sem "sim". Também, linhas depois de `# Echeverria` não contam.
- `release()` em repositório git temporário: tag ausente (com e sem `require`), HEAD na tag, HEAD além da tag, árvore modificada, staged, untracked ignorado, tag anotada, diretório sem git.
- `trailers()`: contagem, coautor humano, "Claude" só no corpo, merge, merge com trailer, `|` no assunto, repositório vazio.
- `rows()` com origem desconhecida, `computed()` e a exclusão de `paper/generated` em `release()`.

## Bloqueante

Nenhum encontrado. Nenhum teste vermelho.

## Notas

1. **`paper/build_tables.py` e `paper/bibliography.bib` mudaram na worktree enquanto eu testava** (14:03, não fui eu; ao começar o `git status` estava limpo). A mudança: `rows()` valida colunas de origem, nova `computed()`, e `release()` ignora `paper/generated`. Acrescentei 6 testes para isso. Os 43 primeiros foram escritos contra a versão commitada e continuam passando contra a nova. A suíte de 567 testes rodou durante a edição, então pode ter misturado as duas versões do script. **Antes do commit, rode a suíte uma vez de novo.**
2. **Vazamento de `mp.mp.dps`:** importar `build_tables.py` fixa `dps = 40` no processo inteiro. Meu fixture importa e restaura o valor. Quem importar o script em outro teste sem isso muda a precisão de tudo o que roda depois. Não escrevi teste que afirme o vazamento; ele só documentaria um defeito.
3. **Fragilidade em `test_f6_labels_match_the_report`:** o CSV arredonda a 3 dígitos e o relatório a 2 (`run.py:420` e `run.py:915`). Hoje os valores são 0.1, 0.01 e 0.001, então coincidem. Uma grade como 0.0125 faria o teste falhar sem haver divergência real. Esse teste lê `outputs/` já gravado, não regenera. Ele detecta divergência entre arquivos, não uma regressão em `run.py`.
4. **`test_figures_lang.py`** compara só `ax.lines`. Dados em coleções (scatter, barras) divergentes entre pt e en não seriam vistos. Hoje as figuras só usam linhas.
5. **`trailers()`:** a checagem é por substring `'Claude' in co`. Um coautor humano chamado "Claude" contaria como assistência de IA. Teste não cobre.
6. **`trailers()` em repositório vazio** devolve zeros em vez de falhar. O teste afirma esse comportamento. Se preferir que falhe, é decisão sua.
7. Não testados: `main()` de ponta a ponta, as tabelas `.tex` geradas e `versions()`. Dependem de `outputs/` completo e do estado do git real.
