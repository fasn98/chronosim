# chronosim · relatório de validação

## Fase 1

Comparação entre a energia do vácuo **calculada** a partir da soma de imagens (Li & Gott 1998, eq. 8) com point-splitting (eq. 11) e a **fórmula fechada** (eq. 12), para um campo escalar conformemente acoplado no espaço de Misner. Unidades naturais (ħ = c = 1).

### Resultado

- Maior diferença no grid (15 pontos, 4 componentes cada), relativa à escala 1/(1440π²ξ⁴): **5.5e-16** [COMPUTED]
- Inclinação de log|ρ| contra log ξ (b = 1, ξ de 10⁻³ a 1): **-4.0** [COMPUTED]; a eq. 12 prevê exatamente −4

| ξ | b | ρ calculado | ρ eq. 12 | erro / escala | traço / escala | cauda de Kummer |
|---|---|---|---|---|---|---|
| 0.5 | 1 | 5.2604116411262466716 | 5.260411641126246671 | 5.5e-16 | 6.17e-19 | 2.09e-5 |
| 0.5 | 3 | 0.061607676766419908814 | 0.061607676766419908813 | 1.05e-18 | 7.14e-21 | 2.09e-5 |
| 0.5 | 2pi | 3.1703248505865360274e-22 | 0.0 | 4.47e-19 | 4.22e-21 | 2.09e-5 |
| 0.5 | 10 | -0.0028509938866864932549 | -0.0028509938866864932551 | 3.39e-19 | 3.12e-21 | 2.09e-5 |
| 0.5 | 30 | -0.0033708742831224759311 | -0.0033708742831224759313 | 2.68e-19 | 2.42e-21 | 2.09e-5 |
| 1 | 1 | 0.32877572757039041695 | 0.32877572757039041694 | 1.37e-16 | 1.54e-19 | 2.09e-5 |
| 1 | 3 | 0.0038504797979012443008 | 0.0038504797979012443008 | 2.63e-19 | 1.79e-21 | 2.09e-5 |
| 1 | 2pi | 4.9536325790414625494e-24 | 0.0 | 1.12e-19 | 1.05e-21 | 2.09e-5 |
| 1 | 10 | -0.00017818711791790582844 | -0.00017818711791790582845 | 8.47e-20 | 7.81e-22 | 2.09e-5 |
| 1 | 30 | -0.00021067964269515474571 | -0.00021067964269515474571 | 6.7e-20 | 6.05e-22 | 2.09e-5 |
| 2 | 1 | 0.020548482973149401059 | 0.020548482973149401059 | 3.43e-17 | 3.86e-20 | 2.09e-5 |
| 2 | 3 | 0.0002406549873688277688 | 0.0002406549873688277688 | 6.57e-20 | 4.46e-22 | 2.09e-5 |
| 2 | 2pi | 7.7400509047522908789e-26 | 0.0 | 2.79e-20 | 2.64e-22 | 2.09e-5 |
| 2 | 10 | -0.000011136694869869114278 | -0.000011136694869869114278 | 2.12e-20 | 1.95e-22 | 2.09e-5 |
| 2 | 30 | -0.000013167477668447171607 | -0.000013167477668447171607 | 1.67e-20 | 1.51e-22 | 2.09e-5 |

### O que a fase 1 estabelece e o que não estabelece

- **Estabelece:** a eq. 12 de Li & Gott decorre das eqs. 8 e 11, reproduzida de forma independente até a precisão numérica, incluindo a divergência ξ⁻⁴ no horizonte (o mecanismo invocado por Hawking, 1992) e o anulamento em b = 2π (o contra-exemplo de Li & Gott).
- **Não estabelece:** que viagens no tempo sejam possíveis ou impossíveis. O cálculo é semiclássico (campo quântico em espaço-tempo fixo, sem reação gravitacional) e Kay, Radzikowski & Wald (1997) provaram que a função de dois pontos não pode ser Hadamard no horizonte de Cauchy; a questão depende de uma teoria de gravidade quântica.

## Fase 2 (parcial): vácuo de Minkowski com imagens

**Comparação com Hiscock & Konkowski (1982): PENDENTE.** O artigo ainda não foi lido, e nenhuma fórmula dele foi usada. Tudo abaixo é só a soma numérica G = Σ_{n≠0} G_M(X, Λⁿ X′) com o point-splitting da eq. 11. As referências a Li & Gott vêm do arXiv, com a auditoria contra o PRL também pendente.

### Resultado

- Maior erro de diferenças finitas no grid (15 pontos, comparando δ = 1e-10, h = 1e-11 com δ = 1e-11, h = 1e-12), relativo a |ρ|: **1.16e-19** [COMPUTED]
- Inclinação de log|ρ| contra log ξ (b = 1, ξ de 10⁻³ a 1): −4 com desvio de **-3.83e-15** [COMPUTED]
- Lei de decaimento medida: exponencial, não n⁻⁴. No último termo somado, −ln(r)/b difere de 1 por no máximo **1.52e-22** [COMPUTED]

| ξ | b | ρ calculado | ξ⁴ρ | erro DF / ρ | (ρ−3p) / ρ | traço / ρ | fração da cauda | N |
|---|---|---|---|---|---|---|---|---|
| 0.5 | 1 | 5.2604116411262557133 | 0.32877572757039098 | 1.16e-19 | 4.58e-20 | 1.37e-22 | 4.29e-26 | 56 |
| 0.5 | 3 | 0.061637733786540367704 | 0.003852358361658773 | 7.68e-21 | 3.12e-21 | 2.25e-23 | 3.0e-28 | 21 |
| 0.5 | 2pi | 0.0015392811901640300987 | 9.6205074385251881e-5 | 7.79e-22 | 6.21e-22 | 9.89e-25 | 9.49e-31 | 11 |
| 0.5 | 10 | 0.000036814835464342067617 | 2.3009272165213792e-6 | 5.34e-22 | 5.35e-22 | 2.42e-26 | 3.97e-31 | 7 |
| 0.5 | 30 | 7.585003482259294486e-14 | 4.7406271764120591e-15 | 5.28e-22 | 5.33e-22 | 4.99e-35 | 8.19e-40 | 3 |
| 1 | 1 | 0.32877572757039098205 | 0.32877572757039098 | 2.9e-20 | 1.14e-20 | 3.43e-23 | 4.29e-26 | 56 |
| 1 | 3 | 0.0038523583616587729815 | 0.003852358361658773 | 1.92e-21 | 7.8e-22 | 5.64e-24 | 3.0e-28 | 21 |
| 1 | 2pi | 0.00009620507438525188117 | 9.6205074385251881e-5 | 1.95e-22 | 1.55e-22 | 2.47e-25 | 9.49e-31 | 11 |
| 1 | 10 | 2.3009272165213792261e-6 | 2.3009272165213792e-6 | 1.34e-22 | 1.34e-22 | 6.05e-27 | 3.97e-31 | 7 |
| 1 | 30 | 4.7406271764120590537e-15 | 4.7406271764120591e-15 | 1.32e-22 | 1.33e-22 | 1.25e-35 | 8.19e-40 | 3 |
| 2 | 1 | 0.020548482973149436378 | 0.32877572757039098 | 7.26e-21 | 2.86e-21 | 8.58e-24 | 4.29e-26 | 56 |
| 2 | 3 | 0.00024077239760367331134 | 0.003852358361658773 | 4.8e-22 | 1.95e-22 | 1.41e-24 | 3.0e-28 | 21 |
| 2 | 2pi | 6.0128171490782425731e-6 | 9.6205074385251881e-5 | 4.87e-23 | 3.88e-23 | 6.18e-26 | 9.49e-31 | 11 |
| 2 | 10 | 1.4380795103258620163e-7 | 2.3009272165213792e-6 | 3.34e-23 | 3.35e-23 | 1.51e-27 | 3.97e-31 | 7 |
| 2 | 30 | 2.9628919852575369086e-16 | 4.7406271764120591e-15 | 3.3e-23 | 3.33e-23 | 3.12e-36 | 8.19e-40 | 3 |

### Os dois vácuos em ξ = 1

D = ρ(Minkowski com imagens) − ρ(eq. 12). A coluna seguinte é uma **observação empírica**, sem derivação nem fonte: D·b⁴·exp(4π²/b) fica quase constante em b ≈ 1–3 e sobe para b maior. Em b = 0,8, D é só cerca de 3 vezes o erro de diferenças finitas, e a coluna não é confiável ali.

| b | ξ⁴ρ Minkowski [COMPUTED] | ξ⁴ρ eq. 12 [DERIVED] | D [COMPUTED] | erro DF | D·b⁴·exp(4π²/b) |
|---|---|---|---|---|---|
| 0.8 | 0.802979376528981 | 0.802979376528981 | 1.0919542e-19 | 3.75e-20 | 120.81978 |
| 1 | 0.328775727570391 | 0.32877572757039 | 5.651168e-16 | 9.54e-21 | 78.958182 |
| 1.5 | 0.0647739638131661 | 0.064773963755243 | 5.7923138e-11 | 7.62e-22 | 78.956835 |
| 2 | 0.0203506032383653 | 0.020350590036348 | 1.3202017e-8 | 1.19e-22 | 78.956837 |
| 3 | 0.00385235836165877 | 0.00385047979790124 | 1.8785638e-6 | 7.4e-24 | 78.958205 |
| 4 | 0.00108997909861866 | 0.00107401894047031 | 1.5960158e-5 | 8.29e-25 | 78.993596 |
| 5 | 0.0003624899022816 | 0.000315293102136562 | 4.71968e-5 | 1.3e-25 | 79.221736 |
| 2pi | 9.62050743852519e-5 | 0.0 | 9.6205074e-5 | 1.87e-26 | 80.291609 |
| 8 | 1.70460987474097e-5 | -0.000130766753022047 | 0.00014781285 | 2.44e-27 | 84.183991 |
| 10 | 2.30092721652138e-6 | -0.000178187117917906 | 0.00018048805 | 3.07e-28 | 93.535039 |
| 30 | 4.74062717641206e-15 | -0.000210679642695155 | 0.00021067964 | 6.26e-37 | 636.23305 |

### O que esta parte estabelece e o que não estabelece

- **Estabelece:** o valor numérico de ⟨T_ab⟩ para a soma de imagens de Minkowski com o operador da eq. 11, com convergência em segunda ordem, a estrutura ρ = 3p sem traço, escala ξ⁻⁴ e decaimento exponencial dos termos.
- **Não estabelece:** que este resultado reproduz Hiscock & Konkowski. Essa é a comparação pendente. Também não estabelece a forma de D(b), que é só um ajuste observado, nem nada sobre viagens no tempo.

## Constantes numéricas escolhidas por nós (HEURISTIC)

- `dps` = 80: working precision in decimal digits (mpmath)
- `delta` = 1.0e-10: point-splitting distance along z, units of length
- `h` = 1.0e-11: finite-difference step (must be < delta), units of length
- `n_terms` = 24: image pairs summed explicitly before the Kummer tail
- `mink_tol` = 1.0e-25: Minkowski image sum stops when |pair_00| <= mink_tol |sum_00|
