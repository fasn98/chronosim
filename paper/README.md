# Artigo para a ReScience C (fase 2, Kim & Thorne 1991)

Rascunho de "[Re] Do vacuum fluctuations prevent the creation of closed timelike curves?".

| Arquivo | O que é |
|---|---|
| `NOTES.md` | Exigências da ReScience C, com a fonte de cada uma, e os riscos da submissão |
| `outline.md` | Esboço aprovado (Etapa 1) |
| `metadata.yaml` | Metadados no formato do template oficial |
| `content.tex` | Corpo do artigo, que o `article.tex` do template inclui |
| `bibliography.bib` | Referências, conforme `chronosim/sources.py` |
| `build_tables.py` | Gera `generated/numbers.tex` e `generated/tab_*.tex` a partir de `outputs/` |

**Nenhum número é digitado à mão:** todos vêm de `generated/numbers.tex`. Para regenerá-los:

```bash
.venv/bin/python -m chronosim.run      # outputs/
.venv/bin/python paper/build_tables.py # paper/generated/
```

**PDF:** ainda não foi gerado. Faltam duas coisas:

- **Arquivos do template** (https://github.com/ReScience/template): `article.tex`, `rescience.cls`, `yaml-to-latex.py`, `Makefile` e as pastas de fontes. Eles ainda não foram copiados para cá; o template é GPL-3+, e as fontes são OFL.
- **Ferramentas:** XeLaTeX, latexmk e biber, que não estão instalados na máquina em que o rascunho foi escrito.

**Figuras:** vêm de `outputs/`, com rótulos em português (pendente; ver a seção Limitations do artigo).
