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

Verificação humana (Fabio, 2026-11-07), em resposta ao achado I6 do reviewer.

| Número | Pág. | Leituras | Status |
|---|---|---|---|
| T⁰⁰_max ~ 10⁴⁹ g/cm³, e a eq. 14 como (l_P/D)(m_P/(Δt)³) com Δt = l_P, igual a (l_P/D)(m_P/l_P³) | 3935 | passada 1 (com zoom), passada 2 (subagente, citação literal) e leitura humana: "is l_P/D ~ 10⁻³⁵ of the Planck density m_P/l_P³; i.e., it is T⁰⁰_max ~ 10⁴⁹ g/cm³"; o expoente 49 é legível sem ambiguidade | conferido |
| m_P ~ 10⁻⁵ g, no texto após a eq. 54 | 3943 | passada 1 e leitura humana | conferido |
| l_P ~ 10⁻³³ cm | 3930 | passada 1 | não conferido na p. 3930; implicado pela razão l_P/D ~ 10⁻³⁵ com D = 1 m declarada na p. 3935 |

O texto após a eq. 14 (p. 3935) declara T⁰⁰_max como l_P/D ~ 10⁻³⁵ vezes a densidade de Planck m_P/l_P³ e dá o valor ~10⁴⁹ g/cm³. Com as entradas declaradas no próprio artigo (m_P ~ 10⁻⁵ g, l_P/D ~ 10⁻³⁵, D = 1 m), essa expressão dá ~10⁵⁹ g/cm³. O valor impresso é inconsistente com a expressão e as entradas que o acompanham. Trata-se de estimativa de ordem de grandeza; nenhuma conclusão do artigo depende dela.

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


# Echeverria, Klinkhammer & Thorne (1991), fase 3

O PDF (`papers/ECHprd91.pdf`) também é um scan sem texto extraível. Seguimos a mesma disciplina da fase 2.

- **Passada 1:** agente principal, com recortes a 250 dpi.
- **Passada 2:** subagente sem acesso à passada 1 nem aos arquivos do projeto, com recortes a 330 dpi e zoom adicional em frações, sinais, linhas e expoentes. A passada 2 não marcou nenhum símbolo como ilegível.

Escopo: as equações usadas na reprodução R1–R11 e as do Apêndice A que servem de checagem.

Notação: `a` = α, `b` = β, `th` = θ, `ph` = φ, `eps` = ε, `lam1` = λ₁, `g1` = γ₁. No Apêndice A, as distâncias latinas a e a′ da Fig. 17 aparecem como `a_` e `a_'` (a passada 2 usou `aL` e `aL'`). As duas passadas diferem só na notação.

| Eq. | Pág. | Passada 1 | Passada 2 | Confere |
|---|---|---|---|---|
| 2.1a | 1082 | `v_out = v_in` | idem | sim |
| 2.1b | 1082 | `psi_out = psi_in - 2 phi` | idem | sim |
| 2.2a | 1082 | `v1' = v2   (speeds)` | idem | sim |
| 2.2b | 1082 | `v2' = v1` | idem | sim |
| 2.2c | 1082 | `v1' + v2' = v1 + v2   (vectors)` | idem | sim |
| 2.3a | 1082 | `v1' = v2, v2' = v1   (vectors)` | idem | sim |
| 2.3b | 1082 | `r1 - r2 = 2R (v2 - v1)/\|v2 - v1\|` | idem | sim |
| 2.4a | 1082 | `v1' = (v2)_reflected in line parallel to v1+v2 ; v2' = (v1)_reflected in line parallel to v1+v2` | idem | sim |
| 2.4b | 1082 | `r1 - r2 = 2 R s (v2 + v1)/\|v2 + v1\|, s = sign(v2 - v1)` | idem; na definição de s, v2 e v1 são módulos (itálico) | sim |
| text | 1083 | `v1' = (D - 2B)/T_d == 1 - 2B` | idem | sim |
| text | 1083 | `v2 = v1 / sin psi ; distance 1/cos psi ; time lapse (1/cos psi)(1/v2) = tan psi / v1` | idem | sim |
| 3.1 | 1083 | `tan psi = v1` | idem | sim |
| 3.2 | 1084 | `sin psi + tan psi = 2 v1` | idem | sim |
| text | 1084 | `T = -1/(1+cos psi) ; T = -cos psi/(1+cos psi) ; T = 0 ; T = cos psi/(1+cos psi) ; T = 1/(1+cos psi)` | idem | sim |
| 3.3 | 1085 | `R << B << 1` | idem | sim |
| 3.4 | 1085 | `v1 > 1, psi_A > B` | idem | sim |
| 3.5 | 1085 | `L >> 1 and L >> h` | idem | sim |
| 3.6 | 1086 | `v2 = v1' = v1 cos(h/L) = v1 [1 - (1/2)(h/L)^2]` | idem | sim |
| 3.7 | 1086 | `n = 2L / (v1 - 1 - (1/2)(h/L)^2 v1)` | idem | sim |
| 3.8 | 1086 | `n_min = 3 sqrt(3 v1/2) \|h\| / (v1 - 1)^(3/2)` | idem | sim |
| 3.9a | 1086 | `psi_(2k+2) - psi_(2k) = -2 phi_(2k+1) for 0 <= k <= n-1` | idem | sim |
| 3.9b | 1086 | `sin phi_(2k+1) - sin phi_(2k-1) = -psi_(2k)/B for 1 <= k <= n-1` | idem | sim |
| 3.10 | 1086 | `psi_0 = psi_A + h/L ; (text) psi_(2n) = psi_0 - pi` | idem | sim |
| 3.11 | 1087 | `v2 = v1' = v1 cos(psi_0 - psi_A) < v1 << 1` | idem | sim |
| 4.1a | 1088 | `psi_A > B` | idem | sim |
| 4.1b | 1088 | `v1 > 1` | idem | sim |
| text | 1089 | `h = -B sin theta, psi_A = theta + phi` | idem | sim |
| 4.2 | 1089 | `B sin a sin((a-b)/2) [sin th - sin(th + ph - (a-b)/2)] + sin b sin((a-b)/2) { sin(th+ph) - B [sin th + sin(th + ph + (a-b)/2)] } = (v1 + d) sin((a+b)/2) sin a sin b` | idem; a passada 2 marca o lado direito "(v1 + d)" como conferido com zoom | sim |
| 4.3 | 1089 | `B sin((a+b)/2 - th - ph)(sin a + sin b) + B sin th sin(a - b) - sin b sin(a - th - ph) = -d sin a sin b` | idem | sim |
| 4.4a | 1089 | `rho = sin[(a+b)/2] / sin[(a-b)/2]` | idem | sim |
| 4.4b | 1089 | `d = 2 s R / (1 + rho^2 + 2 rho cos a)^(1/2), s = sign(d)` | idem | sim |
| 4.5 | 1089 | `v2 = v1 rho` | idem | sim |
| 4.6a | 1089 | `(a,b) = (0,0), (pi,0), (0,pi), (pi,pi), (2 ph, 0)` | idem | sim |
| 4.6b | 1089 | `any solution with rho < 0` | idem | sim |
| 4.6c | 1089 | `any solution with sign(b) != sign(d) == s` | idem | sim |
| 4.7 | 1090 | `B/2 < ph < pi/2, B - ph < th < ph, v1min < v1 < v1max` | idem | sim |
| 4.8 | 1090 | `v1min, v1max = (cos th/cos ph)(1 - 2B cos ph) -/+ 2R/cos ph` | idem | sim |
| 4.9 | 1090 | `\|b\| << 1, \|eps\| << 1, eps == a - 2 ph` | idem | sim |
| 4.10 | 1091 | `v1 = v1min + lam1 (v1max - v1min)` | idem | sim |
| 4.11 | 1091 | `v2 = v1min + lam2 (v1max - v1min)` | idem | sim |
| 4.12a | 1091 | `b = 8 sin ph / (cos th (tan^2 ph - tan^2 th)) B R lam1` | idem | sim |
| 4.12b | 1091 | `eps = 8 cos ph / sin(th + ph) lam1 R` | idem | sim |
| 4.12c | 1091 | `lam2 = (1 + 2 cos ph/(tan^2 ph - tan^2 th) B) lam1` | idem | sim |
| 4.13a | 1091 | `b = -8 sin ph / (cos th (tan^2 ph - tan^2 th)) B R (1 - lam1)` | idem | sim |
| 4.13b | 1091 | `eps = -8 cos ph / sin(th + ph) (1 - lam1) R` | idem | sim |
| 4.13c | 1091 | `1 - lam2 = (1 + 2 cos ph/(tan^2 ph - tan^2 th) B)(1 - lam1)` | idem; o sinal "=" está borrado no scan, mas é "=" | sim |
| 4.14 | 1092 | `b = sin ph cos 2ph B, eps = cos 2ph cos^2 ph / sin ph B` | idem; a passada 2 registra que não há fator lam1, como em B34 | sim |
| 4.15 | 1092 | `b = eps = -4 sin ph / cos 2ph R lam1` | idem | sim |
| A1 | 1095 | `v2 sin((a-b)/2) = v1 sin((a+b)/2)` | idem | sim |
| A2 | 1095 | `\|d\|/v1 = 2R/\|v1 + v2\| = 2R/(v1^2 + v2^2 + 2 v1 v2 cos a)^(1/2)` | idem | sim |
| text | 1095 | `r1(T) = v1 T, r2(T) = v2 T - 2 s R (v1+v2)/\|v1+v2\| ; d2/v2 = d/v1 ; sign(d2) = sign(d) = s ; sign(b) = sign(d) = s ; a'/v2 + (c - d2)/v2 - 1 = 0` | idem | sim |
| A3 | 1096 | `v2 = (a' + c)/(1 + d/v1)` | idem | sim |
| A4 | 1096 | `c = sin(th+ph)/sin a (1 - 2B cos ph + 2B sin(psi/2) cos g1 / sin g2) - 2B sin(psi/2) cos g3 / sin g2` | idem; o último termo −2B fica fora dos parênteses grandes | sim |
| A5 | 1096 | `d = -a + sin g2/sin a (1 - 2B cos ph + 2B sin(psi/2) cos g1 / sin g2)` | idem | sim |
| A6 | 1096 | `a_ = B/sin b [sin(th - b) - sin th']` | idem | sim |
| A7 | 1096 | `a_' = a_ cos[(th + b + th')/2] / cos[(th - b + th')/2]` | idem | sim |
| A8 | 1096 | `psi = th - b - th'` | idem | sim |
| A9 | 1096 | `g1 = (th - 3 th' - b)/2` | idem | sim |
| A10 | 1096 | `g2 = th + ph - 2 th' - b` | idem | sim |
| A11 | 1096 | `g3 = (th - th' - b)/2 + ph` | idem | sim |
| A12 | 1096 | `th' = ph + th - (a + b)/2` | idem | sim |
| A13 | 1097 | `rho == v2/v1` | idem | sim |
| A14 | 1097 | `rho = sin[(a+b)/2]/sin[(a-b)/2]` | idem | sim |
| A15 | 1097 | `d = 2 s R/(1 + rho^2 + 2 rho cos a)^(1/2)` | idem | sim |
| A16 | 1097 | `a_' sin b = B[sin th - sin(b + th')] = B[sin th - sin(th + ph - (a-b)/2)]` | idem | sim |
| A17 | 1097 | `c sin a = sin(th+ph)(1 - 2B cos ph) - 2B sin(psi/2) cos(th + ph + g3)` | idem | sim |
| A18 | 1097 | `psi = -ph + (a-b)/2, g1 = -th - (3/2) ph + (3/4) a + (1/4) b, g2 = a - th - ph, g3 = ph/2 + (a-b)/4` | idem | sim |
| A19 | 1097 | `c sin a = sin(th+ph) - B [sin(th + ph + (a-b)/2) + sin th]` | idem | sim |
| B1 | 1097 | `M eps + N eps^2 + P1 b + Q1 b^2 + S1 b eps = 0` | idem | sim |
| B2 | 1097 | `M eps + N eps^2 + P2 b + Q2 b^2 + S2 b eps = 0` | idem | sim |
| B3 | 1097 | `M = (1/2) B sin 2ph cos th` | idem | sim |
| B4 | 1097 | `N = (1/2) B (cos 2ph cos th + (1/4) sin 2ph sin th)` | idem | sim |
| B5 | 1097 | `P1 = -sin(ph - th) - 2(4 lam1 + s - 2) R sin ph + B[(3/4) sin(2ph - th) - (1/4) sin(2ph + th) - sin th]` | idem | sim |
| B6 | 1098 | `Q1 = -cos ph cos th - (1/2) cot ph sin(ph + th) - 2(2 lam1 - 1) R cos ph + B[(1/2)(9 cos^2 ph - 1) cos th + (1/4) sin th cot ph (7 cos^2 ph - 3)]` | idem | sim |
| B7 | 1098 | `S1 = -cos 2ph cos th / cos ph + R[2(2 lam1 - 1)/cos ph + (4 - 8 lam1 - s) cos ph] + B[cos th ((3/2) cos^2 ph - 1) + (1/4) sin 2ph sin th]` | idem | sim |
| B8 | 1098 | `P2 = -sin(ph - th) + 2 s R sin ph + B(-2 sin th cos^2 ph + (1/2) sin 2ph cos th)` | idem | sim |
| B9 | 1098 | `Q2 = -s R cos ph + B((1/2) cos th - (3/8) sin 2ph sin th)` | idem | sim |
| B10 | 1098 | `S2 = -cos(ph - th) + s R cos ph + B(cos th cos^2 ph + (5/4) sin th sin 2ph)` | idem | sim |
| B11 | 1098 | `v1 = (1/cos ph)[(1 - 2B cos ph) cos th + 2R(2 lam1 - 1)]` | idem | sim |
| B12 | 1098 | `Q b + S eps + P = 0` | idem | sim |
| B13 | 1098 | `Q = -cos ph cos th - (1/2) cot ph sin(ph + th) - (4 lam1 - s - 2) R cos ph + B(-cos th + (9/2) cos^2 ph cos th + cos^2 ph cot ph sin th)` | idem | sim |
| B14 | 1098 | `S = cos th / cos ph - cos(ph + th) + R[2(1 - 2 lam1) cos 2ph / cos ph - 2 s cos ph] + B(-cos th + (1/2) cos th cos^2 ph - sin 2ph sin th)` | idem | sim |
| B15 | 1098 | `P = -4 R sin ph (2 lam1 + s - 1)` | idem | sim |
| B16 | 1098 | `b = sin 2ph cos th / (2 sin(ph - th)) B eps << eps` | idem | sim |
| B17 | 1098 | `eps = 8 cos ph / sin(th + ph) (lam1 - sigma) R` | idem | sim |
| B18 | 1098 | `sigma == (1/2)(1 - s) = 0 if s = +1, 1 if s = -1` | idem | sim |
| B19 | 1098 | `b = 4 sin 2ph cos ph cos th / (sin(ph - th) sin(ph + th)) (lam1 - sigma) B R = 8 sin ph / (cos th (tan^2 ph - tan^2 th)) (lam1 - sigma) B R` | idem | sim |
| B20 | 1098 | `lam2 - sigma = [1 + cos^3 ph cos^2 th / (2 sin(ph - th) sin(ph + th)) B](lam1 - sigma) = (1 + 2 cos ph/(tan^2 ph - tan^2 th) B)(lam1 - sigma)` | idem (as duas formas) | sim |
| B21 | 1098 | `mu == (ph - th)/B` | idem | sim |
| B22 | 1099 | `k1 b + k2 eps + k3 b eps = 0` | idem | sim |
| B23 | 1099 | `k1 = -mu B + 2 s R sin ph - B sin ph cos^2 ph` | idem | sim |
| B24 | 1099 | `k2 = B sin ph cos^2 ph` | idem | sim |
| B25 | 1099 | `k3 = -1` | idem | sim |
| B26 | 1099 | `Q = -2 cos^2 ph` | idem | sim |
| B27 | 1099 | `S = 2 sin^2 ph` | idem | sim |
| B28 | 1099 | `P = -8 R sin ph (lam1 - sigma)` | idem | sim |
| B29 | 1099 | `b^2 + p b + q = 0` | idem | sim |
| B30 | 1099 | `p = 2 tan^2 ph [2(lam1 - sigma) - s sin ph] R - (sin ph cos 2ph - mu tan^2 ph) B` | idem | sim |
| B31 | 1099 | `q = -4 B R sin^2 ph (lam1 - sigma)` | idem | sim |
| B32 | 1099 | `b = -p/2 + sqrt(p^2/4 - q)` | idem | sim |
| B33 | 1099 | `b = -q/p if p >> \|q\| ; b = -p if p << -\|q\|` | idem | sim |
| B34 | 1099 | `b = sin ph cos 2ph B, eps = cos 2ph cos^2 ph / sin ph B   if ph - pi/4 << -R lam1` | idem | sim |
| B35 | 1099 | `b = eps = -4 sin ph / cos 2ph R lam1   if ph - pi/4 >> +R lam1` | idem | sim |

**Resultado:** as 103 entradas conferem.

## Derivações usadas como checagem da transcrição

Cada uma é um teste em `tests/test_ekt_transcription.py`.

- **B13–B15:** iguais a B6 − B9, B7 − B10 e B5 − B8, porque B12 é (B1 − B2)/β. A parte em B de P se cancela, o que confere B5 e B8 termo a termo.
- **B11:** decorre das eqs. 4.8 e 4.10.
- **B26–B28:** são o limite θ → φ, B, R → 0 de B13–B15.
- **B19:** as duas formas impressas são iguais.
- **3.8:** é o mínimo da eq. 3.7 em L.
- **B34 e B35:** são os dois limites de B32 longe de π/4.
- **4.12c:** é a segunda forma da B20.

## Eq. B20: duas formas impressas

As duas formas impressas da eq. B20 diferem por um fator 4 no termo em B. A primeira, cos³φ cos²θ / (2 sin(φ−θ) sin(φ+θ)), vale 1/4 da segunda, 2 cos φ / (tan²φ − tan²θ), porque sin(φ−θ) sin(φ+θ) = sin²φ − sin²θ. A eq. 4.12c usa a segunda. Qual delas a geometria segue é medido em R9 (`outputs/report.md`).

## Eq. B30: o termo em R

Eliminando ε entre as eqs. B22 e B12, com B23–B28, obtive para o termo em R do coeficiente p de B29 a forma 2 tan²φ [2(λ₁ − σ)/sin φ − s sin φ] R. O impresso é 2 tan²φ [2(λ₁ − σ) − s sin φ] R, e as passadas 1 e 2 leem a eq. B30 da mesma forma. O termo em B e a eq. B31 conferem com a minha eliminação.

A geometria independente decide a favor da forma impressa. Na borda superior, com B = 10⁻², R = 10⁻⁵ e λ₁ = ½ (R10), a eq. B32 com o p impresso fica dentro de 0,07–0,3% do β geométrico para φ entre 0,75 e 0,82. Com o p da minha eliminação, o desvio é de 0,9–1,2%. A diferença fica, portanto, na minha eliminação, e não é apontada como inconsistência do artigo. O código usa a forma impressa.

## Eq. 4.6: raízes que a lista não cobre

As raízes não espúrias das eqs. 4.2–4.3 que a geometria independente não encontra têm, todas, alguma distância do Apêndice A negativa (a = PQ, a′ = PQ′ ou c; eqs. A6, A7 e A19). Nenhuma bola percorre distância negativa ao longo do caminho. A eq. 4.6 lista as raízes introduzidas pelas manipulações do Apêndice A (p. 1097), e estas têm outra origem. Ver R8 em `outputs/report.md`.

## Dimensão

As equações estão em unidades D = T_d = 1. Restaurando D e T_d, cada termo precisa ter a mesma dimensão.

- **Eq. 4.2:** o fator "(v₁ + d)" do lado direito soma velocidade e distância, como impresso (a passada 2 conferiu com zoom). Com T_d restaurado, ele é (v₁T_d + d), e a derivação no Apêndice A mostra de onde vem: a eq. A3, v₂ = (a′ + c)/(1 + d/v₁), em que o 1 é T_d (p. 1095, "the time delay ΔT = 1"). Não é inconsistência.
- **Demais equações:** são homogêneas com B, R e d em unidades de D, v em D/T_d e ângulos adimensionais.
- **Teste de escala:** o teste do modelo ("só v₁T_d/D importa"), decidido para o mapa, ainda não foi feito.

## Trechos do texto que parecem trocados

As duas passadas os leem da mesma forma. Nenhum é usado no cálculo.

- **p. 1087:** "decrease the trial-and-error φ₁ toward 2ψ₀", enquanto a p. 1086 dá φ₁ = ½ψ₀[1 + O(B)].
- **p. 1095:** "implies that sin β = s". O contexto, e a relação logo acima, falam de sign(β).
- **p. 1098, após a eq. B20:** "which we quoted and discussed in Sec. IV B [Eqs. (4.14) and (4.15)]". A forma das eqs. B17–B20 é a das eqs. 4.12–4.13.
- **p. 1099, após a eq. B35:** "in (B34) β and φ are independent of λ₁". Pelo contexto, seria ε.

## Citações das lacunas (plano, seção "O que o artigo declara como não demonstrado")

A transcrição dupla foi feita. A passada 1 é a que está no plano (`docs/phase3_plan.md`). A passada 2, literal, está na tabela abaixo.

- **Conferem:** as 20 citações, com a ressalva do item seguinte.
- **Q11, Q12 e Q18:** a passada 2 transcreveu a frase inteira, e o plano cita só um trecho. A frase completa da Q11 continua: ", since the maximum radius then can lie at the wormhole mouth, and the wormhole rather than a collision can be responsible for deflecting the ball back inward toward smaller radii." A Q18a continua: "(i.e., have self-consistent solutions to the equations of motion)."
- **Q14:** é a penúltima frase da legenda da Fig. 12, não a última.
- **Q12:** fica **pendente de conferência humana** até o usuário conferir.

| Citação | Pág. | Passada 2 (literal) |
|---|---|---|
| Q1 | 1079 (right col.) | However, we have searched hard for initial trajectories with zero multiplicity and have found none. On the other hand, our search has not covered all initial trajectories (see especially Sec. V), so we cannot guarantee the nonexistence of zero-multiplicity trajectories. |
| Q2 | 1079 (right col.) | However, we have not found a manageable way to analyze the case of coplanarity to within a distance ~ B. |
| Q3 | 1080 (right col.) | There presumably are also multiple-collision solutions, and we speculate about some possible, rather strange ones in the paragraph containing Eq. (3.11). |
| Q4 | 1080 (right col.) | We suspect, but have not proved, that the (finite-measure) initial trajectories with speeds v_1 << D/T_d ≡ 1 and with impact parameters h >> D ≡ 1 also have unit multiplicity; see Sec. III. |
| Q5 | 1081 (left col., end of Sec. I C) | Our spot checks in these extreme regions have not turned up any initial trajectories for which numerical solutions do not exist. |
| Q6a | 1081 (right col., Sec. II A) | Just how small R must be for this depends on one's model for the internal structure of the ball. |
| Q6b | 1081 (right col.) | (Our choice of B/2 rather than B/4 or 9B/10 is quite arbitrary.) |
| Q8 | 1087 (left col.) | (Readers who seek higher rigor than we do might worry that our analysis has examined only the leading-order effects in the small parameters B, R/B, 1/L, and h/L and has not proved rigorously that higher-order corrections are negligible. We are not worried.) |
| Q7 | 1083 (left col., Sec. III intro) | We also argue, but do not prove firmly, that there is only a single solution for any ball with (i) an initial speed that is sufficiently slow but not zero, and (ii) an initial path of motion that, if extended forever, remains far from the wormhole. |
| Q9 | 1087 (right col.) | We have searched cursorily for such peculiar solutions, without success, and we suspect they do not exist. However, we have no proof. |
| Q10 | 1088 (left col.) | We have not been able to rule out such a solution rigorously, but it seems exceedingly unlikely that one could exist.   [context: "(i) As seen by the ball there is an infinite sequence of self-collisions that goes on and on forever."] |
| Q11 | 1088 (left col.) | Note that this proof fails if the ball is initially at rest inside radius r_min, since the maximum radius then can lie at the wormhole mouth, and the wormhole rather than a collision can be responsible for deflecting the ball back inward toward smaller radii. |
| Q12a | 1088 (left col. bottom -> right col. top) | These equations have solutions in all regimes we have examined (the multiplicity is nonzero), but their high nonlinearity has prevented us from proving definitively that there always is a solution. |
| Q12b | 1088 (right col.) | In Sec. IV C and Appendix B we shall examine the intermediate case R/B < 1/2 but B << 1. In this case we shall show that for a wide range of dangerous initial trajectories there is always at least one self-consistent solution, and we shall argue that this is probably so for all initial trajectories, i.e., the multiplicity is probably always nonzero. |
| Q13a | 1089 (right col.) | It is far from obvious, just looking at the equations, whether there exist values of the wormhole and ball radii B, R and initial-trajectory parameters v_1, θ, φ that produce zero solutions. |
| Q13b | 1089 (right col.) | Numerical exploration, and the analytic considerations of the next section, have not turned up any such zero-multiplicity trajectories. |
| Q14 | 1090 (Fig. 12 caption) | Analytic solutions cannot be derived by the techniques of Appendix B for the white regions of the triangle (left corner and lower left edge), but numerical solutions have been found in spot checks throughout that white region.   [NOTE: this is the sentence about white regions; it is the second-to-last sentence of the caption. The actual last sentence is: "Figure 14 below shows, as an example, a solution (part analytic, part numerical) all along the upper left edge of the dangerous triangle, including the left corner."] |
| Q15a | 1091 (left col.) | At all points near the left corner and lower left edge where we have tried, we found numerical solutions. |
| Q15b | 1091 (left col.) | Thus, it seems likely that solutions exist everywhere in the dangerous triangle. |
| Q16 | 1091 (right col.) | We suspect, but have not proved, that the missing class-II solution actually exists, but the ball first encounters its collision shortly after passing through the wormhole, rather than before, and therefore this solution is beyond the domain of validity of our analysis. |
| Q17 | 1092 (right col., end of Sec. IV) | We find no evidence, when B << 1 and R/B < 1/2, for initial trajectories with zero multiplicity. |
| Q18a | 1092 (right col., Sec. V) | Thereby we shall accumulate evidence which suggests, but does not really prove, that all noncoplanar initial trajectories have multiplicity larger than zero (i.e., have self-consistent solutions to the equations of motion). |
| Q18b | 1092 (right col., Sec. V) | We have not found a good way to analyze dangerous initial trajectories near the boundary of the region (5.1). |
| Q19a | 1094 (left col.) | Thus, we can say nothing about the existence of solutions over the full range. |
| Q19b | 1094 (left col.) | However, these conclusions make us suspect that even when z_1 ~ BR, all initial trajectories will have at least one self-consistent solution.   [preceding sentence: "To recapitulate, these conclusions hold only for a neighborhood of coplanarity (typically z_1 << BR) that is much smaller than the full range of dangerous initial trajectories (typically z_1 ≲ BR)."] |
| Q20 | 1094 (right col., Sec. VI) | On the other hand, we have seen no evidence for a stronger type of ill posedness: generic initial trajectories with no self-consistent solutions. |
