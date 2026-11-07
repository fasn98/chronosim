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

Plano em `docs/phase2_plan.md`, transcrições em `docs/equation_transcription.md`. Até aqui: testes árbitros do operador e itens F1 e F3. F4, F2, F5, F6 e F7 ainda não foram implementados.

### Teste árbitro: eq. 4a de Kim & Thorne contra a eq. 11 de Li & Gott

Operador das eqs. 4a e 4b de Kim & Thorne (p. 3933), como impressas, escrito sem usar o código da fase 1, aplicado à função de Rindler (Li & Gott, eq. 6) regularizada pela eq. 3 (p. 3933). A referência é o operador da eq. 11 de Li & Gott (arXiv gr-qc/9711074, p. 2), validado na fase 1. Escala: 1/(1440π²ξ⁴). Todos os valores são COMPUTED.

| ξ | ρ eq. 4a / escala | ρ eq. 11 / escala | maior desvio da razão em relação a 2 |
|---|---|---|---|
| 0.5 | -5.9999999999999999996 | -2.9999999999999999998 | -7.29e-19 |
| 1.3 | -5.9999999999999999999 | -3.0 | -1.08e-19 |
| 2 | -6.0 | -3.0 | -4.56e-20 |

**Resultado: dá o dobro.** Os dois artigos normalizam a função de Hadamard da mesma forma, 1/(4π²σ) com σ igual à metade do quadrado da separação (Kim & Thorne, eqs. 3, 17 e 20; Li & Gott, eq. 9). A eq. 4a não tem o fator ½ da eq. 11, e nenhuma definição em volta dela o compensa. Inconsistência entre os dois artigos, nas convenções como impressas.

### Árbitro independente dos dois artigos: efeito Unruh

Referência: ρ = −π²T⁴/30 com T = 1/(2πξ), ou seja, −3/(1440π²ξ⁴). A relação é a atribuída a Candelas & Deutsch (1977), com volume e página confirmados no Crossref, mas **o texto do artigo não foi lido** (acesso fechado). A densidade térmica π²T⁴/30 foi conferida por integração numérica da integral de Bose.

| ξ | referência / escala | operador com ½ / referência | eq. 4a impressa / referência |
|---|---|---|---|
| 0.5 | -3.0 | 0.99999999999999999994 | 1.9999999999999999999 |
| 1.3 | -3.0 | 0.99999999999999999999 | 2.0 |
| 2 | -3.0 | 1.0 | 2.0 |

**Resultado:** a versão com ½ concorda com o valor de Unruh, e a eq. 4a impressa dá o dobro. É uma inconsistência de convenção entre a eq. 4a impressa e o valor de Unruh.

### F1: geometria do espaço de recobrimento

Recobrimento construído por iteração do mapa de Poincaré das eqs. 25 e 28, sem usar as eqs. 29–38. Pontos a uma distância ε da garganta, com D = 1. Todos os valores são COMPUTED.

| Afirmação do artigo | Maior desvio medido |
|---|---|
| Origem do recobrimento em T = X = D/(1−ξ) (eq. 29) | 2.11e-81 |
| Cópias U_n = ξ⁻ⁿU₀, V_n = ξⁿV₀ (eq. 33) | 3.9e-70 |
| σ_N = ΔV·ΔU (eq. 34) | 1.48e-70 |
| ΔU (eq. 37) | 1.46e-70 |

**Eq. 31 (p. 3940):** a eq. 31 como impressa é inconsistente com a Fig. 5 e com as eqs. 25, 28 e 29; a versão com o coeficiente de sinal oposto é consistente com as três. Nenhum resultado posterior depende desta equação.

- **Leitura:** coeficiente positivo (1+ξ²ⁿ)/(1−ξ²ⁿ), nas duas transcrições e na verificação humana.
- **Eqs. 25, 28 e 29:** a razão medida entre (T̃+c)/(X̃+c) e o coeficiente impresso é −1, com desvio máximo de 7.38e-80 (n = 1, 2, 3).
- **Fig. 5:** as cópias n ≥ 1 são desenhadas com inclinação negativa, e a eq. 31 como impressa dá inclinação positiva para 0 < ξ < 1.
- **Posição de T_0 e cruzamentos com o horizonte:** a construção reproduz X̃ = −D/(1−ξ) em n = 0, X̃ = −ξⁿD/(1−ξ) nos cruzamentos e os espaçamentos D, Dξ, Dξ² da Fig. 5. **Estes dois fatos não distinguem o sinal:** valem também para o coeficiente positivo, como um teste verifica. Só a inclinação distingue.

Erro relativo da eq. 38 contra o intervalo exato. Ele cai linearmente com ε. O texto apresenta a eq. 35 como aproximação ("very nearly the same", p. 3941), justificada pela proximidade da garganta e do horizonte, mas **não fala em ordem em ε**. A primeira ordem é resultado desta medição e da expressão exata do erro, não uma afirmação do artigo.

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

### F3: ⟨T^μν⟩ entre as bocas

Diferenças finitas aplicadas à eq. 66, com λ_N da eq. 61 e prefator lido literalmente, sem usar a eq. 67, com os dois operadores em paralelo: a eq. 4a impressa (`eq4_stress`, separação em y) e a eq. 4a com o fator ½ (`stress_tensor` da fase 1, separação em z). Δ é a distância a H_1, em unidades de λ_1. D = 1. Tudo é COMPUTED, salvo a coluna da eq. 67, que é DERIVED. A tabela completa está em `outputs/kim_thorne_f3.csv`.

| β | b/D | X/D | ρ/X | Δ/D | sem ½ / eq. 67, T^TT | com ½ / eq. 67, T^TT | maior desvio sem ½ (TT, TX, XX) | demais componentes / eq. 67 |
|---|---|---|---|---|---|---|---|---|
| 0.1 | 0.01 | 0.25 | 0 | 1e-2 | 0.999988803671631 | 0.499994401835816 | 2.24e-5 | 0.000446 |
| 0.1 | 0.01 | 0.25 | 0 | 1e-4 | 0.999999888036716 | 0.499999944018358 | 2.24e-7 | 4.46e-6 |
| 0.1 | 0.01 | 0.25 | 0 | 1e-6 | 0.999999998880367 | 0.499999999440184 | 2.24e-9 | 4.46e-8 |
| 0.1 | 0.1 | 0.25 | 0 | 1e-2 | 0.999988618229504 | 0.499994309114752 | 2.28e-5 | 0.000449 |
| 0.1 | 0.1 | 0.25 | 0 | 1e-4 | 0.999999886182295 | 0.499999943091148 | 2.28e-7 | 4.49e-6 |
| 0.1 | 0.1 | 0.25 | 0 | 1e-6 | 0.999999998861823 | 0.499999999430911 | 2.28e-9 | 4.49e-8 |
| 0.5 | 0.01 | 0.25 | 0 | 1e-2 | 0.999672641338953 | 0.499836320669476 | 0.000655 | 0.00244 |
| 0.5 | 0.01 | 0.25 | 0 | 1e-4 | 0.99999672641339 | 0.499998363206695 | 6.55e-6 | 2.44e-5 |
| 0.5 | 0.01 | 0.25 | 0 | 1e-6 | 0.999999967264134 | 0.499999983632067 | 6.55e-8 | 2.44e-7 |
| 0.5 | 0.1 | 0.25 | 0 | 1e-2 | 0.999668605561985 | 0.499834302780993 | 0.000663 | 0.00245 |
| 0.5 | 0.1 | 0.25 | 0 | 1e-4 | 0.99999668605562 | 0.49999834302781 | 6.63e-6 | 2.45e-5 |
| 0.5 | 0.1 | 0.25 | 0 | 1e-6 | 0.999999966860556 | 0.499999983430278 | 6.63e-8 | 2.45e-7 |
| 0.9 | 0.01 | 0.25 | 0 | 1e-2 | 0.99798971268927 | 0.498994856344635 | 0.00402 | 0.00641 |
| 0.9 | 0.01 | 0.25 | 0 | 1e-4 | 0.999979897126893 | 0.499989948563446 | 4.02e-5 | 6.41e-5 |
| 0.9 | 0.01 | 0.25 | 0 | 1e-6 | 0.999999798971269 | 0.499999899485634 | 4.02e-7 | 6.41e-7 |
| 0.9 | 0.1 | 0.25 | 0 | 1e-2 | 0.997979617743566 | 0.498989808871783 | 0.00404 | 0.00643 |
| 0.9 | 0.1 | 0.25 | 0 | 1e-4 | 0.999979796177436 | 0.499989898088718 | 4.04e-5 | 6.43e-5 |
| 0.9 | 0.1 | 0.25 | 0 | 1e-6 | 0.999999797961774 | 0.499999898980887 | 4.04e-7 | 6.43e-7 |
| 0.5 | 0.01 | 0.25 | 1e-1 | 1e-4 | 0.998588374737475 | 0.499294187368738 | 0.00495 | 0.0455 |
| 0.5 | 0.01 | 0.25 | 1e-1 | 1e-6 | 0.998591575493972 | 0.499295787746986 | 0.00495 | 0.0455 |
| 0.5 | 0.01 | 0.25 | 1e-2 | 1e-4 | 0.999982638819002 | 0.499991319409501 | 4.63e-5 | 0.00455 |
| 0.5 | 0.01 | 0.25 | 1e-2 | 1e-6 | 0.999985879030475 | 0.499992939515238 | 4.96e-5 | 0.00455 |
| 0.5 | 0.01 | 0.25 | 1e-3 | 1e-4 | 0.999996585533217 | 0.499998292766609 | 6.23e-6 | 0.000455 |
| 0.5 | 0.01 | 0.25 | 1e-3 | 1e-6 | 0.999999826384168 | 0.499999913192084 | 4.63e-7 | 0.000455 |

- **Maior desvio do operador sem ½ em relação à eq. 67:** sobre o eixo, 0.00404 (β = 0,9, Δ/D = 1e-2); fora do eixo, 0.00495 (ρ/X = 0,1).
- **Maior erro de diferenças finitas em T^TT:** 3.96e-20.
- **Inclinação de log|T^TT| contra log Δ, sobre o eixo (β = 0,5, b/D = 0,01, X/D = 0,25):** −3 com desvio de -4.39e-6 sem ½ e -4.39e-6 com ½.

**Resultado:** a versão sem ½ reproduz a eq. 67, e a versão com ½ dá metade dela. Os resultados impressos ficam, portanto, 2 vezes acima do valor ancorado pelo árbitro de Unruh. Um fator constante não altera a taxa de divergência, que é a conclusão qualitativa do artigo.

O desvio residual tem a forma dos termos que o artigo despreza. Sobre o eixo, ele cai linearmente com Δ em X/D = 0,25. Em X = D/2, T^TT e T^TX ficam no nível do erro numérico e só T^XX se desvia, quadraticamente, porque a derivada primeira do prefator da eq. 66 se anula ali. Fora do eixo, o desvio em kk cresce como (ρ/X)², e as demais componentes crescem como ρ/X, a correção fracionária que o texto diz ignorar (p. 3945).

### O que esta parte estabelece e o que não estabelece

- **Estabelece:** as eqs. 33, 34 e 37 são exatas para a geometria das eqs. 25 e 28, e a eq. 38 é uma aproximação cujo erro relativo é proporcional a ε/D, sobre o eixo de simetria.
- **Estabelece também:** a eq. 67 decorre das eqs. 61 e 66 pelo operador da eq. 4a impressa, com os desvios medidos acima, e a divergência é Δ⁻³.
- **Não estabelece:** que as eqs. 61 e 66 estejam certas: elas são entrada do cálculo, em óptica geométrica. A garganta (eq. 51) ainda não foi testada. A ancoragem pelo valor de Unruh foi feita no vácuo de Rindler, não no buraco de minhoca, e o texto de Candelas & Deutsch não foi lido.
- **Disputa com Hawking:** Hawking (1992, p. 603 e seção VI) contesta a conclusão de Kim & Thorne. Ele argumenta que um eventual corte deveria ser em distância invariante, e não no tempo próprio, que depende do observador. Esta fase reproduz o cálculo de Kim & Thorne e não decide quem estava certo.

## Constantes numéricas escolhidas por nós (HEURISTIC)

- `kt_series_tol` = 1.0e-30: the eqs. 66-67 series stops when (b xi / 2D)^N is below this
- `dps` = 80: working precision in decimal digits (mpmath)
- `delta` = 1.0e-10: point-splitting distance along z, units of length
- `h` = 1.0e-11: finite-difference step (must be < delta), units of length
- `n_terms` = 24: image pairs summed explicitly before the Kummer tail
- `kt_split_ratio` = 1.0e-10: point-splitting distance / distance to the horizon, Kim & Thorne runs
- `kt_step_ratio` = 1.0e-11: finite-difference step / distance to the horizon, Kim & Thorne runs
