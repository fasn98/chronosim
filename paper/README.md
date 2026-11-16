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
| `measure_timings.py` | Mede os tempos da seção Reproducibility (pytest, execução padrão, `build_tables.py`, `make -C paper`) e grava `generated/timings.json` e `generated/timings.tex`, com CPU, RAM e sistema |

**Nenhum número é digitado à mão:** todos vêm de `generated/numbers.tex`. Para regenerá-los:

```bash
.venv/bin/python -m chronosim.run      # outputs/
.venv/bin/python paper/build_tables.py # paper/generated/
```

**PDF:** `make -C paper` gera `paper/article.pdf`.

- O `Makefile` baixa o template oficial (https://github.com/ReScience/template) fixado no commit `0978d0e`, e copia para `paper/` só o que a compilação usa: `rescience.cls`, `header.tex`, `yaml-to-latex.py`, `article.py` e as pastas de fontes. Esses arquivos não entram no repositório (`paper/.gitignore`): o template é GPL-3+, e as fontes são OFL.
- `article.tex` é o do template, com dois acréscimos: `booktabs`, para as tabelas, e `\emergencystretch`, para evitar linhas que invadem a margem.
- **Ferramentas:** XeLaTeX, latexmk, biber, Python 3 com PyYAML e o locale `en_US.UTF-8`, que o `yaml-to-latex.py` do template exige. No Debian/Ubuntu: `apt-get install texlive-xetex texlive-latex-extra texlive-bibtex-extra biber latexmk texlive-fonts-recommended locales`.
- O PDF não é versionado. Ele acompanha a submissão e a release.

**Tempos (seção Reproducibility):** medidos na máquina que gera o PDF da release, antes da tag:

```bash
python paper/measure_timings.py   # roda pytest, chronosim.run, build_tables.py e make -C paper; grava generated/timings.json e timings.tex
make -C paper                     # de novo, para o PDF incluir os tempos
git add paper/generated/timings.json paper/generated/timings.tex outputs/
git commit -m "Timings measured on the release machine" -m "run: padrão (via paper/measure_timings.py)"
```

O script registra CPU, número de CPUs lógicas, RAM, sistema, kernel e a versão exata do Python que rodou. Se um dos comandos falhar, ele para e não grava nada.

**Figuras:** vêm de `outputs/`, nas versões em inglês (`*_en.png`). Elas são desenhadas a partir dos mesmos dados das versões em português, e `tests/test_figures_lang.py` confere isso.
