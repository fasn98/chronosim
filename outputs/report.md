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

## Fase 2 (em andamento): Kim & Thorne (1991)

Plano em `docs/phase2_plan.md`, transcrições em `docs/equation_transcription.md`. Até aqui: teste árbitro do operador e item F1. F3, F4, F2, F5, F6 e F7 ainda não foram implementados.

### Teste árbitro: eq. 4a de Kim & Thorne contra a eq. 11 de Li & Gott

Operador das eqs. 4a e 4b de Kim & Thorne (p. 3933), como impressas, escrito sem usar o código da fase 1, aplicado à função de Rindler (Li & Gott, eq. 6) regularizada pela eq. 3 (p. 3933). A referência é o operador da eq. 11 de Li & Gott (arXiv gr-qc/9711074, p. 2), validado na fase 1. Escala: 1/(1440π²ξ⁴). Todos os valores são COMPUTED.

| ξ | ρ eq. 4a / escala | ρ eq. 11 / escala | maior desvio da razão em relação a 2 |
|---|---|---|---|
| 0.5 | -5.9999999999999999996 | -2.9999999999999999998 | -7.29e-19 |
| 1.3 | -5.9999999999999999999 | -3.0 | -1.08e-19 |
| 2 | -6.0 | -3.0 | -4.56e-20 |

**Resultado: dá o dobro.** Os dois artigos normalizam a função de Hadamard da mesma forma, 1/(4π²σ) com σ igual à metade do quadrado da separação (Kim & Thorne, eqs. 3, 17 e 20; Li & Gott, eq. 9). A eq. 4a não tem o fator ½ da eq. 11, e nenhuma definição em volta dela o compensa. Inconsistência entre os dois artigos, nas convenções como impressas.

### F1: geometria do espaço de recobrimento

Recobrimento construído por iteração do mapa de Poincaré das eqs. 25 e 28, sem usar as eqs. 29–38. Pontos a uma distância ε da garganta, com D = 1. Todos os valores são COMPUTED.

| Afirmação do artigo | Maior desvio medido |
|---|---|
| Origem do recobrimento em T = X = D/(1−ξ) (eq. 29) | 2.11e-81 |
| Cópias U_n = ξ⁻ⁿU₀, V_n = ξⁿV₀ (eq. 33) | 3.9e-70 |
| σ_N = ΔV·ΔU (eq. 34) | 1.48e-70 |
| ΔU (eq. 37) | 1.46e-70 |

**Discrepância na eq. 31 (p. 3940).** A cópia T_n da garganta, obtida das eqs. 25, 28 e 29, satisfaz a relação da eq. 31 com o sinal oposto do coeficiente impresso: a razão medida entre (T̃+c)/(X̃+c) e (1+ξ²ⁿ)/(1−ξ²ⁿ) é −1, com desvio máximo de 7.38e-80 (n = 1, 2, 3). As duas transcrições leem o sinal +. Nenhuma fórmula usada adiante depende da eq. 31.

Erro relativo da eq. 38 contra o intervalo exato. Ele cai linearmente com ε, como esperado da aproximação da eq. 35:

| β | N | ε = 1e-2 | ε = 1e-4 | ε = 1e-6 | ε = 1e-8 |
|---|---|---|---|---|---|
| 0.1 | 1 | 0.007795 | 7.855e-5 | 7.856e-7 | 7.856e-9 |
| 0.1 | 2 | 0.003883 | 3.898e-5 | 3.898e-7 | 3.898e-9 |
| 0.1 | 3 | 0.002577 | 2.583e-5 | 2.583e-7 | 2.583e-9 |
| 0.1 | 5 | 0.001537 | 1.539e-5 | 1.539e-7 | 1.539e-9 |
| 0.1 | 10 | 0.0007782 | 7.788e-6 | 7.788e-8 | 7.788e-10 |
| 0.5 | 1 | 0.006182 | 6.22e-5 | 6.22e-7 | 6.22e-9 |
| 0.5 | 2 | 0.00316 | 3.17e-5 | 3.17e-7 | 3.17e-9 |
| 0.5 | 3 | 0.002243 | 2.248e-5 | 2.248e-7 | 2.248e-9 |
| 0.5 | 5 | 0.001648 | 1.65e-5 | 1.65e-7 | 1.65e-9 |
| 0.5 | 10 | 0.001421 | 1.423e-5 | 1.423e-7 | 1.423e-9 |
| 0.9 | 1 | 0.00446 | 4.48e-5 | 4.48e-7 | 4.48e-9 |
| 0.9 | 2 | 0.002917 | 2.925e-5 | 2.925e-7 | 2.925e-9 |
| 0.9 | 3 | 0.00264 | 2.647e-5 | 2.647e-7 | 2.647e-9 |
| 0.9 | 5 | 0.002566 | 2.573e-5 | 2.573e-7 | 2.573e-9 |
| 0.9 | 10 | 0.002562 | 2.569e-5 | 2.569e-7 | 2.569e-9 |

### O que esta parte estabelece e o que não estabelece

- **Estabelece:** as eqs. 33, 34 e 37 são exatas para a geometria das eqs. 25 e 28, e a eq. 38 é uma aproximação cujo erro relativo é proporcional a ε/D, sobre o eixo de simetria.
- **Não estabelece:** nada ainda sobre ⟨T_μν⟩ de Kim & Thorne. A eq. 38 é só o ponto de partida dos itens seguintes, e o caso fora do eixo (eqs. 46 e 61) não foi testado. O árbitro mostra o desacordo entre as convenções impressas, mas não diz qual dos artigos o causa.
- **Disputa com Hawking:** Hawking (1992, p. 603 e seção VI) contesta a conclusão de Kim & Thorne. Ele argumenta que um eventual corte deveria ser em distância invariante, e não no tempo próprio, que depende do observador. Esta fase reproduz o cálculo de Kim & Thorne e não decide quem estava certo.

## Constantes numéricas escolhidas por nós (HEURISTIC)

- `dps` = 80: working precision in decimal digits (mpmath)
- `delta` = 1.0e-10: point-splitting distance along z, units of length
- `h` = 1.0e-11: finite-difference step (must be < delta), units of length
- `n_terms` = 24: image pairs summed explicitly before the Kummer tail
