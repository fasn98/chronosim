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

## Fase 2: Kim & Thorne (1991)

Plano em `docs/phase2_plan.md`, transcrições em `docs/equation_transcription.md`. Conteúdo: testes árbitros do operador e itens F1 a F7. Figuras: `fig3_kt_divergence.png` e `fig4_kt_cutoffs.png`. Checagem de literatura do fator 2: `docs/literature_check_factor2.md`.

### Árbitro principal: estado térmico no espaço de Minkowski

Não depende de nenhum artigo fechado. A função de Hadamard é a soma de imagens em tempo imaginário (período β = 1/T), derivada da condição KMS em `docs/thermal_images.md`, sem o termo n = 0, que é o vácuo de Minkowski. A referência é Stefan-Boltzmann para um escalar sem massa, π²T⁴/30, conferido pela integral de Bose. Valores COMPUTED.

| T | operador | ρ / (π²T⁴/30) | maior desvio de 3p/ρ em relação a 1 | traço / ρ | erro DF |
|---|---|---|---|---|---|
| 0.3 | eq. 4a impressa | 1.999999999999999999936 | -2.5e-20 | 8.36e-23 | 3.16e-20 |
| 0.3 | com ½ | 0.9999999999999999999682 | -2.5e-20 | 2.09e-23 | 3.14e-20 |
| 1 | eq. 4a impressa | 1.999999999999999999936 | -2.5e-20 | 8.36e-23 | 3.16e-20 |
| 1 | com ½ | 0.9999999999999999999682 | -2.5e-20 | 2.09e-23 | 3.14e-20 |
| 3 | eq. 4a impressa | 1.999999999999999999936 | -2.5e-20 | 8.36e-23 | 3.16e-20 |
| 3 | com ½ | 0.9999999999999999999682 | -2.5e-20 | 2.09e-23 | 3.14e-20 |

Os valores coincidem entre as temperaturas porque δ e h escalam com β.

**Resultado:** o operador com ½ reproduz Stefan-Boltzmann, e a eq. 4a impressa dá o dobro. O estado tem pressão ρ/3 e traço zero com os dois operadores. É uma inconsistência de convenção entre a eq. 4a impressa e o valor térmico.

A soma de imagens também foi comparada com a integral de modos com ocupação de Bose, em pontos separados (β = 1, 400 pares e o termo dominante da cauda):

| t | r | imagens | modos | diferença relativa |
|---|---|---|---|---|
| 0.1 | 0.2 | 0.159568016514404 | 0.159568016477603 | 2.31e-10 |
| 0.3 | 0.5 | 0.130419330870684 | 0.130419330597301 | 2.1e-9 |
| 0.7 | 1.5 | 0.0780197619218726 | 0.0780197599661581 | 2.51e-8 |

A diferença vem do corte da soma de imagens. Birrell & Davies (1982) fica como referência pendente: a página não pôde ser conferida, e a derivação está escrita por extenso em `docs/`.

### Comparação direta: eq. 4a de Kim & Thorne contra a eq. 11 de Li & Gott

Operador das eqs. 4a e 4b de Kim & Thorne (p. 3933), como impressas, escrito sem usar o código da fase 1, aplicado à função de Rindler (Li & Gott, eq. 6) regularizada pela eq. 3 (p. 3933). A referência é o operador da eq. 11 de Li & Gott (arXiv gr-qc/9711074, p. 2), validado na fase 1. Escala: 1/(1440π²ξ⁴). Todos os valores são COMPUTED.

| ξ | ρ eq. 4a / escala | ρ eq. 11 / escala | maior desvio da razão em relação a 2 |
|---|---|---|---|
| 0.5 | -5.9999999999999999996 | -2.9999999999999999998 | -7.29e-19 |
| 1.3 | -5.9999999999999999999 | -3.0 | -1.08e-19 |
| 2 | -6.0 | -3.0 | -4.56e-20 |

**Resultado: dá o dobro.** Os dois artigos normalizam a função de Hadamard da mesma forma, 1/(4π²σ) com σ igual à metade do quadrado da separação (Kim & Thorne, eqs. 3, 17 e 20; Li & Gott, eq. 9). A eq. 4a não tem o fator ½ da eq. 11, e nenhuma definição em volta dela o compensa. Inconsistência entre os dois artigos, nas convenções como impressas.

### Árbitro secundário: efeito Unruh

Referência: ρ = −π²T⁴/30 com T = 1/(2πξ), ou seja, −3/(1440π²ξ⁴). **Origem HEURISTIC:** a relação "vácuo de Rindler = −térmico" é citada de Candelas & Deutsch (1977), artigo não lido (volume e página confirmados só no Crossref; acesso fechado). O árbitro térmico continua sendo o principal.

| ξ | referência / escala | operador com ½ / referência | eq. 4a impressa / referência |
|---|---|---|---|
| 0.5 | -3.0 | 0.99999999999999999994 | 1.9999999999999999999 |
| 1.3 | -3.0 | 0.99999999999999999999 | 2.0 |
| 2 | -3.0 | 1.0 | 2.0 |

**Resultado:** a versão com ½ concorda com o valor de Unruh, e a eq. 4a impressa dá o dobro, como no árbitro térmico. É uma inconsistência de convenção entre a eq. 4a impressa e o valor de Unruh.

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

Diferenças finitas aplicadas à eq. 66, com λ_N da eq. 61 e prefator lido literalmente, sem usar a eq. 67, com os dois operadores em paralelo: a eq. 4a impressa (`eq4_stress`, separação em y) e a eq. 4a com o fator ½ (`stress_tensor` da fase 1, separação em z). Δ é a distância a H_1, em unidades de λ_1. D = 1. Tudo é COMPUTED, inclusive a coluna da eq. 67: ela é a série impressa, truncada em kt_series_tol. A tabela completa está em `outputs/kim_thorne_f3.csv`.

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

- **Maior desvio do operador sem ½ em relação à eq. 67:** sobre o eixo, 0,00404 (β = 0,9, Δ/D = 1e-2); fora do eixo, 0,00495 (ρ/X = 0,1, Δ/D = 1e-4).
- **Maior erro de diferenças finitas em T^TT:** 3.96e-20.
- **Inclinação de log|T^TT| contra log Δ, sobre o eixo (β = 0,5, b/D = 0,01, X/D = 0,25):** −3 com desvio de -4.39e-6 sem ½ e -4.39e-6 com ½.

**Resultado:** Os resultados impressos de Kim & Thorne (eq. 67) são reproduzidos pelo operador da eq. 4a como impresso, sem o fator ½. Portanto o fator é carregado de forma consistente até o resultado final, e os valores impressos ficam 2 vezes acima dos obtidos com a convenção que reproduz o árbitro térmico. A taxa de divergência (Δ⁻³) e a troca de sinal na hipersuperfície não dependem desse fator.

**Literatura** (`docs/literature_check_factor2.md`): não foi encontrada errata nem quem tenha apontado o fator. Klinkhammer (1992, tese, eqs. 3.3 e 3.8) e Visser (hep-th/9303023, eqs. 10 e 34) usam o operador com ½ e chegam à metade da eq. 67 (identidade verificada em `tests/test_literature_factor2.py`), mas declaram concordância com Kim & Thorne sem comentar o fator.

O desvio residual tem a forma dos termos que o artigo despreza. Sobre o eixo, ele cai linearmente com Δ em X/D = 0,25. Em X = D/2, T^TT e T^TX ficam no nível do erro numérico e só T^XX se desvia, quadraticamente, porque a derivada primeira do prefator da eq. 66 se anula ali. Fora do eixo, o desvio em kk cresce como (ρ/X)², e as demais componentes crescem como ρ/X, a correção fracionária que o texto diz ignorar (p. 3945).

### F4: ⟨T^μν⟩ na garganta

Diferenças finitas aplicadas à eq. 50, com λ_N da eq. 46, sem usar a eq. 51, com os dois operadores em paralelo, em pontos da eq. 42 a uma distância Δ antes de H_1, em unidades de λ_1. A comparação usa as 16 componentes do tensor, e o desvio é relativo à maior componente da eq. 51. Para a eq. 52, roda a forma impressa √(X²+Y²) e a forma √(Y²+Z²), o ρ da eq. 39. Tudo é COMPUTED, inclusive as colunas da eq. 51 (série impressa truncada em kt_series_tol). A tabela completa está em `outputs/kim_thorne_f4.csv`.

| β | b/D | θ | φ | Δ/b | sem ½ / eq. 51, T^TT | com ½ / eq. 51, T^TT | desvio com √(Y²+Z²) | desvio com √(X²+Y²) |
|---|---|---|---|---|---|---|---|---|
| 0.1 | 0.01 | 0.2 | 0 | 1e-4 | 1.0 | 0.5 | 1.865e-6 | 0.001047 |
| 0.1 | 0.01 | 0.7 | 0 | 1e-4 | 1.0 | 0.5 | 3.925e-6 | 0.06137 |
| 0.1 | 0.01 | 1.3 | 0 | 1e-4 | 1.0 | 0.5 | 9.61e-6 | 0.5779 |
| 0.5 | 0.01 | 0.2 | 0 | 1e-2 | 1.0 | 0.5 | 0.001233 | 0.001198 |
| 0.5 | 0.01 | 0.2 | 0 | 1e-4 | 1.0 | 0.5 | 1.232e-5 | 0.001198 |
| 0.5 | 0.01 | 0.2 | 0 | 1e-6 | 1.0 | 0.5 | 1.232e-7 | 0.001198 |
| 0.5 | 0.01 | 0.7 | 0 | 1e-2 | 1.0 | 0.5 | 0.002593 | 0.089 |
| 0.5 | 0.01 | 0.7 | 0 | 1e-4 | 1.0 | 0.5 | 2.593e-5 | 0.09157 |
| 0.5 | 0.01 | 0.7 | 0 | 1e-6 | 1.0 | 0.5 | 2.593e-7 | 0.0916 |
| 0.5 | 0.01 | 1.3 | 0 | 1e-2 | 1.0 | 0.5 | 0.004254 | 0.5712 |
| 0.5 | 0.01 | 1.3 | 0 | 1e-4 | 1.0 | 0.5 | 4.254e-5 | 0.5778 |
| 0.5 | 0.01 | 1.3 | 0 | 1e-6 | 1.0 | 0.5 | 4.254e-7 | 0.5779 |
| 0.9 | 0.01 | 0.2 | 0 | 1e-4 | 1.0 | 0.5 | 4.155e-5 | 0.0015 |
| 0.9 | 0.01 | 0.7 | 0 | 1e-4 | 1.0 | 0.5 | 6.069e-5 | 0.1332 |
| 0.9 | 0.01 | 1.3 | 0 | 1e-4 | 1.0 | 0.5 | 7.757e-5 | 0.5778 |
| 0.5 | 0.01 | 0.7 | 0.6 | 1e-2 | 1.0 | 0.5 | 0.00225 | 0.1232 |
| 0.5 | 0.01 | 0.7 | 0.6 | 1e-4 | 1.0 | 0.5 | 2.25e-5 | 0.1209 |
| 0.5 | 0.01 | 0.7 | 0.6 | 1e-6 | 1.0 | 0.5 | 2.25e-7 | 0.1209 |

- **Maior desvio do operador sem ½, com √(Y²+Z²), em todo o grid:** 0.007757 (Δ/b = 1e-2); em Δ/b = 1e-6, 7.757e-07.
- **Maior desvio com a leitura literal √(X²+Y²):** 0.5779; em Δ/b = 1e-6, 0.5779. Não cai com Δ.
- **Maior erro de diferenças finitas em T^TT:** 1.64e-16.
- **Inclinação de log|T^TT| contra log Δ (β = 0,5, b/D = 0,01, θ = 0,7):** −3 com desvio de -9.6e-13 sem ½ e -9.6e-13 com ½.

**Resultado:** Os resultados impressos de Kim & Thorne (eq. 51) são reproduzidos pelo operador da eq. 4a como impresso, sem o fator ½. Portanto o fator é carregado de forma consistente até o resultado final, e os valores impressos ficam 2 vezes acima dos obtidos com a convenção que reproduz o árbitro térmico. A taxa de divergência (Δ⁻³) e a troca de sinal na hipersuperfície não dependem desse fator.

**Eq. 52 (conferida):** a eq. 52 como impressa tem radicando X²+Y²; isso é inconsistente com a igualdade ρ̂ = ∂/∂ρ da própria equação, com a definição de ρ (eq. 39), com a descrição de ρ̂ como vetor unitário (p. 3943) e com a eq. 51, que só é reproduzida com √(Y²+Z²). O cálculo usa √(Y²+Z²). Com a forma impressa, a componente transversal de ρ̂ na garganta, para φ = 0, é cos(θ/2), e o desvio cresce com θ sem cair com Δ.

O desvio com √(Y²+Z²) cai linearmente com Δ/b. Ele é o termo com a segunda derivada de λ_N nas direções transversais, que a eq. 51 não contém: a expressão analítica desse termo reproduz o resíduo em YY, YZ e ZZ, e as demais componentes ficam no nível do erro numérico. Por isso T^TT coincide com a eq. 51 nos 15 dígitos da tabela.

### F2: desfocalização por óptica geométrica

A garganta é tratada como a lente fina divergente de distância focal b/2 que o texto descreve (p. 3942), e a amplitude de uma fonte pontual é 1/B, com B o elemento da matriz de transferência. Isso usa a imagem do artigo, mas não as eqs. 47, 48 e 64. A ligação com Δ_N^(1/2) é a eq. 24. COMPUTED.

| b/D | garganta, ψ/eq. 47 − 1 (N = 1, 2, 3, 5) | entre as bocas, X/D = 0,25, ψ/(eq. 64/ζ_N) − 1 (N = 1, 2, 3, 5) |
|---|---|---|
| 1e-2 | 0.0, -0.009901, -0.01968, -0.03895 | -0.02597, -0.03562, -0.04514, -0.06391 |
| 1e-4 | 0.0, -9.999e-5, -0.0002, -0.0003999 | -0.0002666, -0.0003666, -0.0004665, -0.0006664 |
| 1e-6 | 0.0, -1.0e-6, -2.0e-6, -4.0e-6 | -2.667e-6, -3.667e-6, -4.667e-6, -6.667e-6 |

As eqs. 47, 48 e 64 são o termo dominante em b/D: a correção é linear em b/D e cresce com N. A eq. 57 coincide com o lapso exato de V obtido do espaço de recobrimento sobre o horizonte, com desvio máximo de 4.22e-81.

### F5: hipersuperfícies polarizadas e horizonte de Cauchy

Raiz numérica de λ_N(x, x) = 0 (eqs. 46 e 61): coincide com a eq. 43 até 1.05e-80 e com a eq. 62 até 1.81e-77 (relativo à distância de H_N a T = X). O limite N → ∞, avaliado em N = 60, se aproxima das eqs. 44a e 63 como ξ^N:

| β | desvio de N = 60 contra a eq. 44a | desvio contra a eq. 63 |
|---|---|---|
| 0.1 | 0.00487 | 0.00325 |
| 0.5 | 9.71e-15 | 6.48e-15 |
| 0.9 | 8.68e-39 | 5.79e-39 |

### F6: quadrimomento

Quadratura numérica das eqs. 51 e 67 (validadas em F3 e F4) nos domínios escritos no artigo, em unidades de Planck (l_P = m_P = 1), com D = 1000, b = 10, X = 250 e β = 0,5. Por integrar as eqs. 51 e 67 impressas, estes resultados carregam a convenção sem ½; com ½, cada um cai à metade. COMPUTED.

- **Eq. 69:** desvio máximo de 1.99e-8 em β = 0,1, 0,5, 0,9, N = 1 e 3 e ρ/X = 0, 1e-3, 1e-2; independe de ρ, como o texto afirma (testado em ρ/X = 0, 1e-3, 1e-2; a independência também é algébrica).
- **Eq. 68** (N = 1 e 2, ΔT = l_P): desvio de -5.34e-11.
- **Eq. 54** (N = 1): a integral numérica se aproxima da eq. 54 quando ΔT/b diminui: ΔT/b = 0.1: desvio 0.6318; ΔT/b = 0.01: desvio 0.07867; ΔT/b = 0.001: desvio 0.009066. A eq. 54 é o termo dominante; o termo em tan θ (cortado em a/b = 1e-6) e o fator cos θ somam a correção.
- **Eq. 53** (N = 1, ΔT = l_P): o coeficiente de ln(b/a) da integral numérica difere do da eq. 53 por -1.25e-9. A razão total fica abaixo de 1 pelo termo finito, que a eq. 53 não inclui: a/b = 1e-2: 0.983216; a/b = 1e-4: 0.991602; a/b = 1e-6: 0.994401.

### F7: conferência de aritmética de estimativas de ordem de grandeza

**Isto confere só a aritmética entre os valores de entrada impressos e os números impressos, não se as estimativas estão corretas.** Entradas: l_P/D ~ 10⁻³⁵ com D = 1 m (p. 3935), de onde vem l_P ~ 10⁻³³ cm, e m_P ~ 10⁻⁵ g (p. 3943).

**Conferência humana das leituras:**

- **T⁰⁰_max:** o valor impresso 10⁴⁹ g/cm³ foi conferido na p. 3935.
- **Eq. 14:** a forma (l_P/D)(m_P/l_P³) foi conferida na p. 3935.
- **m_P:** o valor 10⁻⁵ g foi conferido na p. 3943.
- **l_P:** o valor 10⁻³³ cm não foi conferido na p. 3930. Ele é derivado da razão l_P/D ~ 10⁻³⁵ com D = 1 m, conferida na p. 3935.

| Grandeza | Fórmula | Valor impresso | Recalculado |
|---|---|---|---|
| l_P / D (D = 1 m) | l_P / D | 1.0e-35 (Kim & Thorne (1991), p. 3935, after eq. 14) | 1.0e-35 |
| (δL/L)_max, eq. 13, N − 1 = 0 | l_P / D | 1.0e-35 (Kim & Thorne (1991), p. 3935, after eq. 13) | 1.0e-35 |
| densidade de Planck m_P/l_P³ | m_P / l_P³ | 1.0e+94 (Hawking (1992), p. 603: 10^94 g/cc) | 1.0e+94 |
| T⁰⁰_max, eq. 14 | (l_P/D)(m_P/l_P³) | 1.0e+49 (Kim & Thorne (1991), p. 3935, after eq. 14) | 1.0e+59 |
| Δt do corte de Hawking (δg^VP = 1, eq. 9c) | l_P² / D, em l_P | 1.0e-35 (Kim & Thorne (1991), p. 3937: Delta t ~ 10^-35 l_P) | 1.0e-35 |

O texto após a eq. 14 (p. 3935) declara T⁰⁰_max como l_P/D ~ 10⁻³⁵ vezes a densidade de Planck m_P/l_P³ e dá o valor ~10⁴⁹ g/cm³. Com as entradas declaradas no próprio artigo (m_P ~ 10⁻⁵ g, l_P/D ~ 10⁻³⁵, D = 1 m), essa expressão dá ~10⁵⁹ g/cm³. O valor impresso é inconsistente com a expressão e as entradas que o acompanham. Trata-se de estimativa de ordem de grandeza; nenhuma conclusão do artigo depende dela.

As demais quantidades conferem. A densidade de Planck de 10⁹⁴ g/cm³ que Hawking imprime (1992, p. 603) coincide com m_P/l_P³ a partir das entradas de Kim & Thorne.

### O que esta parte estabelece e o que não estabelece

- **Estabelece:** as eqs. 33, 34 e 37 são exatas para a geometria das eqs. 25 e 28, e a eq. 38 é uma aproximação cujo erro relativo é proporcional a ε/D, sobre o eixo de simetria.
- **Estabelece também:** as eqs. 67 e 51 decorrem das eqs. 61 e 66 e das eqs. 46 e 50 pelo operador da eq. 4a impressa, com os desvios medidos acima, e a divergência é Δ⁻³. As eqs. 43, 44a, 62, 63 e 69 decorrem das expressões de partida; as eqs. 47, 48, 53, 54, 64 e 68 são os termos dominantes dos cálculos numéricos correspondentes.
- **Não estabelece:** que as eqs. 46, 50, 61 e 66 estejam certas: elas são entrada do cálculo, em óptica geométrica. A ancoragem pelos árbitros foi feita no estado térmico e no vácuo de Rindler, não no buraco de minhoca. O que os árbitros mostram é qual convenção reproduz esses valores.
- **Disputa com Hawking:** Hawking (1992, p. 603 e seção VI) contesta a conclusão de Kim & Thorne. Ele argumenta que um eventual corte deveria ser em distância invariante, e não no tempo próprio, que depende do observador. A `fig4_kt_cutoffs.png` mostra os dois cortes lado a lado. Esta fase reproduz o cálculo de Kim & Thorne e não decide quem estava certo.
- **F7** confere aritmética de estimativas de ordem de grandeza; não testa as estimativas.

## Constantes numéricas escolhidas por nós (HEURISTIC)

- `ekt_theta_grid` = 64: grid points on the visible arc of a mouth in the inner (aim) search
- `ekt_nu_grid` = 48: grid points of the contact-normal angle in the middle search
- `ekt_polish_dps` = 40: working precision when a bracketed root is polished with findroot
- `ekt_verify_tol` = 1.0e-25: tolerance of the verifier (distances, speeds, residuals) at polish precision
- `ekt_max_split` = 12: halvings of an s-interval whose root count changes before it is reported as unresolved
- `ekt_jump_ratio` = 1e-06: a bracketed sign change is a root only if the refined |f| falls below this fraction of the end values (rejects jumps)
- `ekt_bisect_rel_tol` = 1e-13: bisection stops at this fraction of the bracket width (floats); with mpmath, 10^-(dps-8)
- `ekt_float_depth` = 11: floats are used while k log10(2/B) is below this; mpmath above (itinerary windows ~ (B/2)^k)
- `ekt_future_traversals` = 6: traversals followed after the collision when the verifier builds the history
- `ekt_s_grid` = 240: points of the uniform part of the outer grid of collision positions
- `ekt_s_far` = 3: farthest collision position searched, in D, along the incoming line (before the traversal) or the exit line (after)
- `ekt_s_near_over_R` = '0.01': closest collision position to the line's mouth-2 entry searched before the traversal, as a fraction of R
- `ekt_s_near_after` = '1e-6': closest collision position to the mouth-1 exit searched after the traversal, in D
- `ekt_refine_max_split` = 30: further halvings of an unresolved nu- or s-interval in the adaptive refinement of the empty case
- `ekt_window_grid` = 160: extra collision positions in each refinement window around a known special point
- `ekt_cert_chunks` = 8: pieces of the s range of the empty case run as separate interval-exclusion tasks (parallelism only)
- `ekt_refine_factor` = 2: grid refinement factor of the stabilization search (empty case)
- `kt_series_tol` = 1.0e-30: the eqs. 66-67 series stops when (b xi / 2D)^N is below this
- `dps` = 80: working precision in decimal digits (mpmath)
- `delta` = 1.0e-10: point-splitting distance along z, units of length
- `h` = 1.0e-11: finite-difference step (must be < delta), units of length
- `n_terms` = 24: image pairs summed explicitly before the Kummer tail
- `kt_split_ratio` = 1.0e-10: point-splitting distance / distance to the horizon, Kim & Thorne runs
- `kt_step_ratio` = 1.0e-11: finite-difference step / distance to the horizon, Kim & Thorne runs
- `kt_quad_dps` = 15: working precision of the F6 integrals
- `kt_throat_thickness` = 1.0e-6: a / b, angular cutoff tan(theta) < b / a for eqs. 53-54
- `thermal_n_terms` = 24: imaginary-time image pairs summed explicitly before the Kummer tail
- `thermal_split_ratio` = 1.0e-10: point-splitting distance / beta in the thermal arbiter
- `thermal_step_ratio` = 1.0e-11: finite-difference step / beta in the thermal arbiter
- `thermal_mode_dps` = 30: working precision of the mode integral
- `thermal_mode_intervals` = 30: subintervals of the mode integral
- `thermal_mode_pairs` = 400: image pairs summed in the image-versus-mode comparison
- `thermal_mode_cutoff` = 60: mode integral stops at k = cutoff / beta, where n_k ~ e^-60
- `kt_f5_limit_n` = 60: N at which the N -> infinity limits of F5 (eqs. 44a, 63) are evaluated

## Fase 3 · Bolas de bilhar (Echeverria, Klinkhammer & Thorne, 1991): execução reduzida

Esta é a execução padrão, com um subconjunto pequeno de R1–R11. A reprodução completa está em `outputs/phase3_full.md`, escrita pela última execução de `python -m chronosim.run --full`; os números de lá não foram recalculados nesta execução. Todas as grandezas são **COMPUTED**.

- **R1, eq. 2.1b:** desvio máximo 1,38e-40 rad.
- **R2, eqs. 2.2–2.4:** troca de velocidades 144, troca espelhada 144, outra 0.
- **Bola em repouso (p. 1084):** verificador aceita em x = 0,3; 0,5; 0,7: sim; sim; sim.
- **R3, Fig. 6b** (B = 0,01, R = 1e-4): 1 solução; rapidez menos (1 − 2B): -0,000178; 3 intervalos não resolvidos.
- **R8, caso central (0,2; 0,6; 0,5; 0,01; 0,001):** 3 soluções geométricas (I, II, II), 3 raízes físicas das eqs. 4.2–4.3; maior resíduo das eqs. 4.2–4.3 ou desvio da A6: 1,48e-37; sign(β) = sign(d) em 3/3; 0 intervalos não resolvidos.
- **R10, borda superior, φ = 0,8:** β classe I = 0,00023052; eq. B32 = 0,00023069 (B = 10⁻²).
- **R11, classe II depois da travessia, φ = 0,85:** 1 soluções; β = -0,0007575.

**Não estabelece** nada além destes casos; o que a fase estabelece está em `outputs/phase3_full.md`.

