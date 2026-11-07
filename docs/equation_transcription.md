# Transcrição das equações de Kim & Thorne (1991)

O PDF (`papers/KIMprd91.pdf`) é um scan sem texto extraível, então cada equação usada é transcrita duas vezes, em passadas separadas, a partir de recortes ampliados.

- **Passada 1:** agente principal, recortes a 250 dpi.
- **Passada 2:** subagente sem acesso à passada 1 nem aos arquivos do projeto, recortes a 330 dpi.

As duas passadas vêm do mesmo tipo de modelo, então são separadas, mas não equivalem a dois revisores humanos. Erros de leitura comuns aos dois não seriam detectados aqui, e é por isso que as derivações (F3, F4) também servem de checagem.

**Escopo atual:** equações do teste árbitro e de F1. As demais entram aqui antes do item que as usa.

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

**Resultado:** as 20 equações conferem. Os recortes das eqs. 17, 25a, 28, 29 e 38 erraram o alvo na primeira tentativa e foram refeitos nas duas passadas antes de contar como lidos.

## Dimensão

`tests/test_kim_thorne_geometry.py` tem um teste de análise dimensional para cada fórmula implementada (eqs. 3, 17, 29, 30, 31, 32, 33, 36, 37, 38 e o operador da eq. 4): com c = G = ħ = 1, todo comprimento escala com λ, e a fórmula precisa escalar com a potência de λ da sua dimensão.

## Definições em volta do fator ½ (árbitro)

Lidas com zoom para procurar um fator ½ equivalente. **Nenhum foi encontrado.**

- **Eq. 1 (p. 3932):** define G^(1) como a soma de dois termos, impressos na mesma ordem, φ̂(x)φ̂(x′) + φ̂(x)φ̂(x′). A eq. 6 de Li & Gott (arXiv, p. 2) escreve φ(X)φ(X′) + φ(X′)φ(X), com as ordens trocadas. Os dois têm dois termos, sem fator na frente.
- **Eq. 3 (p. 3933):** "σ₀ is half the square of the infinitesimal vector separation". Dá 1/(4π²σ₀), igual a 1/(2π²·separação²), a eq. 9 de Li & Gott. O teste árbitro confirma numericamente que as duas normalizações coincidem.
- **Eqs. 4a e 4b (p. 3933):** T_μν = lim D_μν G_reg, com coeficientes 2/3, −1/3 e −1/6. A eq. 11 de Li & Gott (arXiv, p. 2) tem os mesmos coeficientes e um ½ na frente.
- **Eq. 20 (p. 3938):** coeficiente Δ_N^(1/2)/(4π²), o mesmo da eq. 17 para Δ = 1.

Os dois artigos atribuem essa fórmula a referências que não lemos: ref. 15 em Kim & Thorne, refs. 24 e 25 em Li & Gott.
