# Checagem de literatura: o fator 2 da eq. 4a de Kim & Thorne (1991)

Pergunta: alguém já apontou que a eq. 4a de Kim & Thorne (Phys. Rev. D 43, 3929), sem o fator ½, dá resultados 2 vezes acima da convenção que reproduz o estado térmico (`outputs/report.md`)? Há errata?

**Resposta curta:** não encontrei errata nem ninguém que aponte o fator. Encontrei dois trabalhos que usam um operador equivalente ao com ½ e chegam, para a divergência perto de uma geodésica nula fechada, a uma expressão que é **metade** da eq. 67 de Kim & Thorne. Mesmo assim, os dois dizem estar de acordo com Kim & Thorne. Essa comparação depende das nossas transcrições de Klinkhammer e Visser, que tiveram leitura única.

## O que foi lido

### Klinkhammer (1992), tese de doutorado, Caltech

Baixada para `papers/Klinkhammer_g_1992.pdf`. A URL indicada (`thesis.library.caltech.edu`) não respondeu; o mesmo arquivo veio de `thesis.caltech.edu/6633/1/Klinkhammer_g_1992.pdf`. Tem 72 páginas, digitalizadas com OCR. As equações abaixo foram lidas nas imagens das páginas 30 a 32 da tese (páginas 36 a 38 do PDF), em leitura única.

- **Função de Hadamard:** G¹(X, X′) = ⟨φ(X)φ(X′) + φ(X′)φ(X)⟩ (p. 29), com singularidade Δ^(1/2)/(4π²σ) (eq. 3.2, p. 30). É a mesma normalização de Kim & Thorne (eqs. 17 e 20) e de Li & Gott.
- **Operador (eq. 3.3, p. 30):**
  D_ab = (1/6)(∇_a′∇_b + ∇_a∇_b′) − (1/12) g_ab ∇_c∇^c′ − (1/12)(∇_a∇_b + ∇_a′∇_b′) + (1/48) g_ab(∇_c∇^c + ∇_c′∇^c′) − R_ab + (1/4) g_ab R.
  Para G simétrico, em espaço plano e numa solução da equação de onda, isto é exatamente ½ × [ (2/3)∇_a∇_b′ − (1/3)∇_a∇_b − (1/6) g_ab ∇_c∇^c′ ]: **o operador com ½**. Teste numérico: `tests/test_literature_factor2.py`, aplicado a um par de imagens do estado térmico, com desvio máximo de 4,0e-20 contra `pointsplit.stress_tensor`.
- **Resultado (eqs. 3.7 e 3.8, pp. 31–32):** ⟨T_ab⟩ = −1/[12π²(δt)³] Σ_(n≥1) (Δ_n^(1/2)/ζ_n)(1 + 4e^(nh) + e^(2nh))/(e^(nh) − 1)³ k_a k_b.
- **Comparação:** com Δ_n^(1/2)/ζ_n da eq. 64 de Kim & Thorne e e^h = 1/ξ (eq. 33), a eq. 3.8 dá **exatamente metade** da eq. 67 de Kim & Thorne. É uma identidade algébrica, e o mesmo teste a verifica. O cálculo numérico com o operador com ½ concorda com a eq. 3.8 até o termo de ordem Δ/X que a eq. 67 também despreza.
- **O que a tese diz:** "their results for the scalar vacuum polarization, when specialized to the vicinity of the closed null geodesics of their spacetimes, are in perfect agreement with ours" (seção 4, p. 38 da tese, p. 44 do PDF). **A tese não comenta o fator 2.** A concordância declarada vale na forma (δt⁻³, estrutura de fluido nulo, dependência em h), não no coeficiente.
- **Versão publicada:** Visser cita a versão publicada como G. Klinkhammer, Phys. Rev. D 46, 3388 (1992). Não a li.

### Visser, arXiv:hep-th/9303023

Lido no texto do arXiv, sem conferir a versão publicada.

- **Operador (eq. 10):** os mesmos coeficientes da eq. 3.3 de Klinkhammer, que ele cita como fonte. É o operador com ½.
- **Resultado (eq. 34):** −Σ′ Δ^(1/2)/(24π²ζ) (1 + 4e^(nh) + e^(2nh))/(e^(nh) − 1)³ l_μl_ν/(δt)³, com a soma sobre as duas orientações. É igual à eq. 3.8 de Klinkhammer, como ele próprio afirma.
- **O que ele diz:** "Furthermore this result is consistent with that of Kim and Thorne [3] — their equation (67)." **Não comenta o fator 2.**

## Buscas que não encontraram nada

| Busca | Resultado |
|---|---|
| Crossref, registro de 10.1103/PhysRevD.43.3929 | sem `updated-by`, `update-to` ou relação de errata; 156 citações registradas |
| Crossref, busca bibliográfica por "erratum" com o título | nenhuma errata; os resultados são artigos não relacionados |
| Página do artigo na APS | nenhum aviso de errata ou correção encontrado |
| Busca na web: errata ou "factor of 2" com Kim & Thorne | nenhuma errata; nenhuma menção ao fator |
| Busca na web: Kim & Thorne com point-splitting, Hadamard, "factor of two", Visser, Klinkhammer | encontrou Visser hep-th/9303023, sem menção ao fator |
| Busca na web: operador 2/3, −1/3, −1/6 com Misner, Sushkov, Tanaka, Hiscock | nada específico sobre a eq. 4a |
| INSPIRE: 175 trabalhos que citam Kim & Thorne; 131 com arXiv, todos baixados como texto | 99 mencionam Kim & Thorne; nenhum tem "factor of 2", "factor of two", "twice" ou "half" junto de Kim & Thorne em contexto de T_μν; nenhum fala de erro, engano ou fator ausente junto de Kim & Thorne |
| Os mesmos 131: referências às eqs. 4, 51 ou 67 | só Visser (gr-qc/9311026), que compara a eq. 102 dele com as eqs. 48 e 64 de Kim & Thorne (determinante de Van Vleck), sem tocar no operador |

**Limites desta busca:** 44 dos 175 trabalhos citantes não têm versão no arXiv e não foram lidos. A busca no texto é por palavras e pode perder uma menção escrita de outra forma (por exemplo, numa equação). A versão publicada de Klinkhammer (PRD 46, 3388) não foi lida.

## Redação

Como não encontramos quem tenha apontado o fator, dentro dos limites acima, a redação do relatório não muda: este trabalho mostra a inconsistência de convenção de forma independente. Acrescenta-se que Klinkhammer (1992) e Visser (hep-th/9303023) chegam, com o operador com ½, à metade da eq. 67, e declaram concordância com Kim & Thorne sem comentar o fator.
