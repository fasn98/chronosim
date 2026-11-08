# Handoff: claude → reviewer (artigo ReScience C, branch paper/rescience-kim-thorne)

Worktree: `/home/fasn98/chronosim/.claude/worktrees/paper-rescience`. Escopo: `git -C /home/fasn98/chronosim/.claude/worktrees/paper-rescience diff main..HEAD` (ae66a33, 9424490, e61c2a4). Somente leitura.

## Entregue

- **Rascunho do artigo:** `paper/content.tex`, `paper/metadata.yaml`, `paper/bibliography.bib`.
- **Gerador:** `paper/build_tables.py`, que produz `paper/generated/*.tex`.
- **Notas:** `paper/NOTES.md` (exigências da ReScience) e `paper/outline.md`.
- **`chronosim/run.py`:** CSVs das tabelas e figuras `_en`.
- **`LICENSE`:** MIT.
- **Contexto:** a revisão do usuário já pediu as correções do commit 9424490.

## O que falta (o que revisar)

1. **Honestidade:**
   - Toda afirmação sobre Kim & Thorne tem página e equação, e confere com `docs/equation_transcription.md` e com `outputs/report.md`.
   - Nenhum número do texto é digitado à mão: todos vêm de macros de `paper/generated/numbers.tex`.
   - Cada número traz sua origem onde isso importa ao leitor.
   - O texto não afirma mais do que os outputs mostram. Atenção a "Estabelece", "the residuals have the form of the neglected terms" e às frases sobre Klinkhammer e Visser.
2. **Linguagem:** neutra. Nada de "misprint", "typo" ou juízo sobre os autores.
3. **Gerador:**
   - Cada macro sai da coluna e da linha certas do CSV.
   - Os filtros escolhem as linhas que o texto descreve.
   - A regra da tag está correta.
4. **Código novo em `run.py`:**
   - Os rótulos dos CSVs conferem com o relatório.
   - As origens estão corretas.
   - As figuras `_en` usam os mesmos dados.
5. Classificar cada achado como **bloqueante** ou **nota**. Resposta final: o conteúdo de `docs/handoffs/2026-11-08-reviewer-para-claude-3.md`.

## Riscos

- A citação de Klinkhammer (p. 38) está pendente da conferência humana, e o texto diz isso.
- Os valores impressos da tabela F7 foram conferidos no scan por leitura do agente, não por leitura humana.
