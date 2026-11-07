# Transcrição das equações de Kim & Thorne (1991)

O PDF (`papers/KIMprd91.pdf`) é um scan sem texto extraível, então cada equação usada é transcrita duas vezes, em passadas separadas, a partir de recortes ampliados.

- **Passada 1:** agente principal, recortes a 250 dpi.
- **Passada 2:** subagente sem acesso à passada 1 nem aos arquivos do projeto, recortes a 330 dpi.

As duas passadas vêm do mesmo tipo de modelo, então são separadas, mas não equivalem a dois revisores humanos. Erros de leitura comuns aos dois não seriam detectados aqui, e é por isso que as derivações (F3, F4) também servem de checagem.

**Escopo:** todas as equações usadas na fase 2 (árbitros e F1 a F7).

Notação: `Xt` e `Tt` são X̃ e T̃, `==` é o sinal de identidade impresso, e `phi^` é φ̂.

| Eq. | Pág. | Passada 1 | Passada 2 | Confere |
|---|---|---|---|---|
| 1 | 3932 | `G^(1)(x,x') = <Psi| phi^(x) phi^(x') + phi^(x) phi^(x') |Psi>` | `G^(1)(x, x') = <Psi | phi^(x) phi^(x') + phi^(x) phi^(x') | Psi>` | sim |
| 3 | 3933 | `G^(1)_reg = G^(1) - 1/(4 pi^2 sigma_0)` | `G^(1)_reg = G^(1) - 1/(4 pi^2 sigma_0)` | sim |
| 4a | 3933 | `T_(mu nu) = lim_(x->x') D_(mu nu) G^(1)_reg` | `T_(mu nu) = lim_(x -> x') D_(mu nu) G^(1)_reg` | sim |
| 4b | 3933 | `D_(mu nu) = (2/3) nabla_mu nabla_(nu') - (1/3) nabla_mu nabla_nu - (1/6) g_(mu nu) nabla_alpha nabla^(alpha')` | idem; marcou as frações como legíveis, porém pequenas no scan | sim |
| 17 | 3937 | `G^(1) = 1/(4 pi^2 sigma)` | `G^(1) = 1 / (4 * pi^2 * sigma)` | sim |
| 20 | 3938 | `G^(1) = sum_(N=-inf)^(+inf) (Delta_N^(1/2)/(4 pi^2)) (1/sigma_N + v_N ln|sigma_N| + w_N)` | idem | sim |
| 25a | 3939 | `T = tau, X = Y = Z = 0` | `T = tau, X = Y = Z = 0` | sim |
| 25b | 3939 | `T = D + gamma tau, X = D - beta gamma tau, Y = Z = 0` | idem | sim |
| 26 | 3939 | `X' = gamma(X-D) + beta gamma(T-D), Y' = Y, Z' = Z; T' = gamma(T-D) + beta gamma(X-D)` | idem | sim |
| 28 | 3939 | `(T=tau, X, Y, Z) on mouth 1 is the same as (T'=tau, X'=-X, Y'=Y, Z'=Z) on mouth 2` | idem | sim |
| 29 | 3940 | `Tt = T - D/(1-xi), Xt = X - D/(1-xi)` | idem | sim |
| 30 | 3940 | `xi == ((1-beta)/(1+beta))^(1/2)` | idem | sim |
| 31 | 3940 | `(Tt + xi^n D/(1-xi)) = ((1+xi^(2n))/(1-xi^(2n))) (Xt + xi^n D/(1-xi))` | idem | sim |
| 32 | 3940 | `U == Xt - Tt, V == (1/2)(Xt + Tt)` | idem | sim |
| 33 | 3940 | `U_n(x) = xi^(-n) U_0(x), V_n(x) = xi^n V_0(x)` | idem | sim |
| 34 | 3941 | `sigma_N = [V_N(x) - V_0(x')][U_N(x) - U_0(x')]` | idem | sim |
| 35 | 3941 | `V_N(x) - V_0(x') = zeta_N` | idem | sim |
| 36 | 3941 | `zeta_N = sum_(n=0)^(N-1) D xi^n = ((1-xi^N)/(1-xi)) D` | idem | sim |
| 37 | 3941 | `U_N(x) - U_0(x') = (X-T) xi^(-N) - (X'-T')` | idem | sim |
| 38 | 3941 | `sigma_N = zeta_N[(X-T) xi^(-N) - (X'-T')]` | idem | sim |
| 39 | 3941 | `rho == sqrt(Y^2+Z^2), rho' == sqrt(Y'^2+Z'^2)` | idem | sim |
| 42 | 3941 | `rho = b sin(theta), X + b = b cos(theta)` | idem | sim |
| 46 | 3942 | `lambda_N(x,x') == 2(b - sqrt(b^2-rho^2)) + X - T - (X'-T') xi^N` | idem | sim |
| 50 | 3942 | `G^(1)_reg = sum_(N=1)^inf (xi/(4 pi^2 D)) (b xi/(2D))^(N-1) (1/lambda_N(x,x') + 1/lambda_N(x',x))` | idem | sim |
| 51 | 3943 | `T^(mu nu) = (-xi/(6 pi^2 D)) sum_(N=1)^inf (b xi/(2D))^(N-1) [(1+4 xi^N+xi^(2N)) k^mu k^nu + (2+4 xi^N) tan(theta)(k^mu rhohat^nu + rhohat^mu k^nu) + 4 tan^2(theta) rhohat^mu rhohat^nu] / [b(1-cos(theta))(1+xi^N) - T(1-xi^N)]^3` | idem | sim |
| 9a | 3934 | `T^(mu nu) ~ (b/D)^(N or (N-1)) (l_P/D) (m_P/(Delta t)^3)` | idem | sim |
| 9b | 3934 | `R_(mu rho nu sigma) ~ G T_(mu nu) ~ (b/D)^(N or (N-1)) (l_P/D) (l_P/(Delta t)^3)` | idem | sim |
| 9c | 3934 | `delta L / L ≃ delta g^VP_(mu nu) ~ (b/D)^(N or (N-1)) (l_P/D) (l_P/Delta t)` | idem | sim |
| 10 | 3935 | `delta R^QG >~ l_P/(Delta t)^3, delta g^QG >~ l_P/Delta t` | idem | sim |
| 11 | 3935 | `delta g^QG / delta g^VP >~ (D/l_P) (D/b)^(N or (N-1))` | idem | sim |
| 12 | 3935 | `T^00 ~ -(b/D)^(N or (N-1)) (l_P/D) (m_P/(Delta t)^3)` | idem | sim |
| 13 | 3935 | `(delta L/L)_max ~ (delta g^VP)_max ~ (b/D)^(N or (N-1)) (l_P/D)` | idem | sim |
| 14 | 3935 | `T^00_max ~ (l_P/D)(m_P/(Delta t)^3)_(Delta t = l_P) = (l_P/D)(m_P/l_P^3)` | idem | sim |
| 15a | 3936 | `Delta M_1 ~ (1/G)(b^2/D)` | idem | sim |
| 15b | 3936 | `Delta M_- ~ (b/D) m_P` | idem | sim |
| 24 | 3939 | `Delta_N^(1/2) = zeta_N psi_N` | idem | sim |
| 43 | 3941 | `T_HN(theta) == ((1+xi^N)/(1-xi^N)) b (1 - cos theta), 0 <= theta < pi/2` | idem | sim |
| 44a | 3941 | `T_H- = lim T_HN = b (1 - cos theta), 0 <= theta < pi/2` | idem | sim |
| 47 | 3942 | `psi_N = (1/D) (b/(2D))^(N-1)` | idem | sim |
| 48 | 3942 | `Delta_N^(1/2) = (zeta_N/D) (b/(2D))^(N-1)` | idem | sim |
| 53 | 3943 | `P^mu_N == int int int T^(mu nu) n_nu b^2 sin theta = (-2b/(3 pi)) (l_P/Delta T)^2 ln(b/a) (b xi/(2D))^N (1+2 xi^N)/(1-xi^N)^3 k^mu` | idem | sim |
| 54 | 3943 | `P^mu_- == ... = (-m_P/(6 pi)) (l_P/Delta T) sum_N (b/(2D))^N (1+4 xi^N+xi^(2N))/(1-xi^N)^2 k^mu` | idem | sim |
| 57 | 3944 | `zeta_N = D - X + sum_(n=1)^(N-1) D xi^n + X xi^N = ((1-xi^N)/(1-xi)) D - (1-xi^N) X` | idem | sim |
| 62 | 3944 | `T_HN == X + rho^2/(2X(D-X)) (D/(1-xi^N) - X)` | idem | sim |
| 63 | 3944 | `T_H- == X + rho^2/(2X)` | idem | sim |
| 64 | 3945 | `Delta_N^(1/2) = (zeta_N D/(X(D-X))) (b/(2D))^N` | idem | sim |
| 68 | 3945 | `P^mu_- == int_0^inf int_(-inf)^(T_H- - Delta T) T^(mu X) dT 2 pi rho drho = (-m_P/(6 pi)) (l_P/Delta T) sum_N (b/(2D))^N (1+4 xi^N+xi^(2N))/(1-xi^N)^2 k^mu` | idem | sim |
| 69 | 3946 | `(dP^mu/dA)_N == int_(-inf)^(t_HN - l_P) T^(mu X) dT = (-D/(12 pi^2 X(D-X))) (b xi/(2D))^N (1+4 xi^N+xi^(2N))/(1-xi^N)^3 k^mu` | idem; a passada 2 nota que o limite superior é impresso com t minúsculo | sim |
| 52 | 3943 | `rhohat = d/drho = (1/sqrt(X^2+Y^2)) (Y d/dY + Z d/dZ); k = (d/dX + d/dT)` | idem; marcou o radicando como legível, a confirmar | sim |
| 61 | 3944 | `lambda_N(x,x') = (X-T) - (X'-T') xi^N + rho^2/(2X) + (rho'^2/(2(D-X'))) xi^N` | idem | sim |
| 66 | 3945 | `G^(1)_reg = sum_(N=1)^inf (D/(4 pi^2 X (D-X))) (b xi/(2D))^N (1/lambda_N(x,x') + 1/lambda_N(x',x))` | idem | sim |
| 67 | 3945 | `T^(mu nu) = (-D k^mu k^nu/(6 pi^2 X (D-X))) sum_(N=1)^inf (b xi/(2D))^N (1 + 4 xi^N + xi^(2N)) / {(X-T)(1-xi^N) + (rho^2/2X) + [rho^2/2(D-X)] xi^N}^3` | idem | sim |

**Resultado:** as 52 equações conferem. Os recortes das eqs. 12 (expoente), 13, 15b, 17, 24, 25a, 28, 29, 38, 43, 50, 57, 61, 62 e 67 erraram o alvo na primeira tentativa e foram refeitos nas duas passadas antes de contar como lidos.

**Eq. 52: conferida.** Ver a seção "Verificações humanas". O código mantém as duas formas do denominador para a comparação de F4, e o cálculo usa √(Y²+Z²).

**Derivação como checagem da transcrição:** as eqs. 61 e 66, aplicadas ao operador da eq. 4a, fecham na eq. 67, e as eqs. 46 e 50 fecham na eq. 51 (ver `outputs/report.md`, F3 e F4). Isso confere as leituras dessas seis equações e a da eq. 4b.

## Verificações humanas

**Terceira leitura, independente e humana (2026-11-06).** Fabio conferiu visualmente no PDF as eqs. 1, 3, 4a e 4b, e as quatro conferem com as transcrições: a eq. 4a não tem o fator ½; a eq. 4b tem os coeficientes 2/3, −1/3 e −1/6; a eq. 3 subtrai 1/(4π²σ₀), com σ₀ igual à metade do quadrado da separação; a eq. 1 imprime os dois termos na mesma ordem.

**Eq. 31 (p. 3940, 2026-11-06).** Fabio confirmou que a eq. 31 está impressa com coeficiente positivo (1+ξ²ⁿ)/(1−ξ²ⁿ).

**Eq. 52 (p. 3943), terceira leitura, independente e humana (2026-11-07).** Fabio confirmou que o denominador está impresso como √(X²+Y²), com numerador (Y ∂/∂Y + Z ∂/∂Z).

## Números impressos usados em F7

**Aguardando conferência humana.** As três leituras abaixo sustentam a comparação 10⁴⁹ contra 10⁵⁹ g/cm³ (reviewer, I6).

| Número | Pág. | Leituras |
|---|---|---|
| l_P ~ 10⁻³³ cm | 3930 | passada 1 (agente principal) |
| m_P ~ 10⁻⁵ g | 3943 | passada 1 (agente principal) |
| T⁰⁰_max ~ 10⁴⁹ g/cm³ | 3935 | passada 1 (com zoom) e passada 2 (subagente, citação literal) |

## Eq. 52

A eq. 52 como impressa tem radicando X²+Y²; isso é inconsistente com a igualdade ρ̂ = ∂/∂ρ da própria equação, com a definição de ρ (eq. 39), com a descrição de ρ̂ como vetor unitário (p. 3943) e com a eq. 51, que só é reproduzida com √(Y²+Z²). O cálculo usa √(Y²+Z²).

| Evidência | Origem |
|---|---|
| ρ̂ = ∂/∂ρ com ρ = √(Y²+Z²) (eq. 39): a regra da cadeia dá ∂/∂ρ = (Y ∂/∂Y + Z ∂/∂Z)/√(Y²+Z²) | a própria equação; teste numérico em `tests/test_kim_thorne_f4.py` |
| "unit transverse vector" (p. 3943); com √(X²+Y²), a componente transversal na garganta é cos(θ/2) | o texto |
| só √(Y²+Z²) reproduz a eq. 51; com √(X²+Y²), o desvio chega a 0,58 e não cai com Δ | F4 |

## Eq. 31

A eq. 31 como impressa é inconsistente com a Fig. 5 e com as eqs. 25, 28 e 29; a versão com o coeficiente de sinal oposto é consistente com as três. Nenhum resultado posterior depende desta equação.

| Evidência | Distingue o sinal? |
|---|---|
| Construção das eqs. 25, 28 e 29: a razão medida contra o coeficiente impresso é −1 | sim |
| Fig. 5: cópias n ≥ 1 desenhadas com inclinação negativa; a eq. 31 impressa dá inclinação positiva para 0 < ξ < 1 | sim |
| X̃ = −D/(1−ξ) em n = 0, posição de T_0 declarada no texto logo antes da eq. 31 | **não** |
| Cruzamentos com o horizonte em X̃ = −ξⁿD/(1−ξ), com os espaçamentos D, Dξ, Dξ² da Fig. 5 | **não** |

Os dois últimos fatos são verificados a partir da construção, mas valem para os dois sinais: em n = 0 o coeficiente diverge, e no horizonte T̃ = X̃ a equação só exige X̃ + c = 0. Um teste registra isso. Eles mostram que a versão de sinal oposto é consistente, sem excluir a impressa.

## Dimensão

`tests/test_kim_thorne_geometry.py` tem um teste de análise dimensional para cada fórmula implementada (eqs. 3, 17, 29, 30, 31, 32, 33, 36, 37, 38 e o operador da eq. 4): com c = G = ħ = 1, todo comprimento escala com λ, e a fórmula precisa escalar com a potência de λ da sua dimensão. As eqs. 61, 66 e 67 têm os seus em `tests/test_kim_thorne_f3.py`, e as eqs. 46, 50, 51 e 52 em `tests/test_kim_thorne_f4.py`.

## Definições em volta do fator ½ (árbitro)

Lidas com zoom para procurar um fator ½ equivalente. **Nenhum foi encontrado**, e a leitura humana acima confirma as eqs. 1, 3, 4a e 4b.

- **Eq. 1 (p. 3932):** define G^(1) como a soma de dois termos, impressos na mesma ordem, φ̂(x)φ̂(x′) + φ̂(x)φ̂(x′). A eq. 6 de Li & Gott (arXiv, p. 2) escreve φ(X)φ(X′) + φ(X′)φ(X), com as ordens trocadas. Os dois têm dois termos, sem fator na frente.
- **Eq. 3 (p. 3933):** "σ₀ is half the square of the infinitesimal vector separation". Dá 1/(4π²σ₀), igual a 1/(2π²·separação²), a eq. 9 de Li & Gott. O teste árbitro confirma numericamente que as duas normalizações coincidem.
- **Eqs. 4a e 4b (p. 3933):** T_μν = lim D_μν G_reg, com coeficientes 2/3, −1/3 e −1/6. A eq. 11 de Li & Gott (arXiv, p. 2) tem os mesmos coeficientes e um ½ na frente.
- **Eq. 20 (p. 3938):** coeficiente Δ_N^(1/2)/(4π²), o mesmo da eq. 17 para Δ = 1.

Os dois artigos atribuem essa fórmula a referências que não lemos: ref. 15 em Kim & Thorne, refs. 24 e 25 em Li & Gott.
