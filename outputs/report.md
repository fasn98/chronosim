# chronosim · fase 1 · relatório de validação

Comparação entre a energia do vácuo **calculada** a partir da soma de imagens (Li & Gott 1998, eq. 8) com point-splitting (eq. 11) e a **fórmula fechada** (eq. 12), para um campo escalar conformemente acoplado no espaço de Misner. Unidades naturais (ħ = c = 1).

## Resultado

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

## Constantes numéricas escolhidas por nós (HEURISTIC)

- `dps` = 80: working precision in decimal digits (mpmath)
- `delta` = 1.0e-10: point-splitting distance along z, units of length
- `h` = 1.0e-11: finite-difference step (must be < delta), units of length
- `n_terms` = 24: image pairs summed explicitly before the Kummer tail

## O que isto estabelece e o que não estabelece

- **Estabelece:** a eq. 12 de Li & Gott decorre das eqs. 8 e 11, reproduzida de forma independente até a precisão numérica, incluindo a divergência ξ⁻⁴ no horizonte (o mecanismo invocado por Hawking, 1992) e o anulamento em b = 2π (o contra-exemplo de Li & Gott).
- **Não estabelece:** que viagens no tempo sejam possíveis ou impossíveis. O cálculo é semiclássico (campo quântico em espaço-tempo fixo, sem reação gravitacional) e Kay, Radzikowski & Wald (1997) provaram que a função de dois pontos não pode ser Hadamard no horizonte de Cauchy; a questão depende de uma teoria de gravidade quântica.