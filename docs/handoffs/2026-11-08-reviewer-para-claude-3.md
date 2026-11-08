# Handoff: reviewer → claude (artigo ReScience C, branch paper/rescience-kim-thorne)

Gravado por claude a partir do relato do reviewer, que é somente-leitura. **Veredito: BLOQUEADO: 3 achados bloqueantes**, mais 9 notas.

## Entregue

- **Escopo:** `git diff main..paper/rescience-kim-thorne` (ae66a33, 9424490, e61c2a4). O reviewer só leu: não editou nada e não rodou testes, nem `run`, nem `build_tables.py`.
- **Conferido e correto:**
  - páginas e equações de Kim & Thorne, contra `docs/equation_transcription.md`;
  - as verificações humanas citadas;
  - todas as macros, uma a uma, contra coluna e linha dos CSVs;
  - os filtros do gerador e a regra da tag;
  - os CSVs novos e suas origens;
  - as figuras `_en`, que vêm dos mesmos dados;
  - a linguagem, neutra;
  - as citações de Klinkhammer e Visser.

## O que falta (achados)

- **BLOQUEANTE 1** (`content.tex:109`, `build_tables.py:176`): `PhaseOneTwoPiRel` vinha de `max_err_over_scale`, e não de |ρ|/escala. O valor correto é 2,8×10⁻¹⁹, não 4,47×10⁻¹⁹.
- **BLOQUEANTE 2** (`content.tex:241`, `:258`, `:146-151`; `metadata.yaml:53`): o texto afirma mais do que os outputs mostram.
  - "Every closed form ... reproduced" é falso para a eq. 31.
  - "Declares neglected" só se sustenta para F3 fora do eixo.
  - O resíduo só foi identificado analiticamente em F4.
  - "Every closed-form result" não inclui as eqs. 10–12 e 15a–b.
- **BLOQUEANTE 3** (`content.tex:30`, `:44`): falta Visser gr-qc/9311026, que compara as eqs. 48 e 64. "Cover none" não foi verificado para a tese inteira.
- **Notas:**
  - N1: a condição da igualdade do operador de Klinkhammer não está dita.
  - N2: `generated/` commitado faz a contagem de commits ficar defasada, e a regra da árvore limpa recusa um segundo build.
  - N3: dois valores impressos de F7 só tiveram a leitura do agente.
  - N4: o razão com ½ em F3 foi medida em Δ/D = 10⁻⁶, e o texto não diz.
  - N5: "to within" seguido de número negativo.
  - N6: a eq. 63 é anunciada em F5 e não é relatada.
  - N7: N = 60 digitado no CSV de F5.
  - N8: as origens dos CSVs são strings digitadas, e o gerador não as confere.
  - N9: `NOTES.md` está desatualizado sobre a LICENSE.

## Riscos

- O reviewer não executou nada.
- Os scans não foram conferidos pelo reviewer.
- A renderização dos PNGs não foi inspecionada.
