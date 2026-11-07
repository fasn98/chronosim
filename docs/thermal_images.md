# Estado térmico: imagens em tempo imaginário a partir da condição KMS

Derivação usada pelo árbitro principal do fator ½ (`chronosim/thermal.py`). Ela é feita por extenso aqui porque a página de Birrell & Davies (1982) não pôde ser conferida; essa referência está marcada como pendente em `sources.py`.

Campo escalar livre sem massa, com ħ = c = 1, ρ̂ = e^(−βH)/Z e φ(t, **x**) = e^(iHt) φ(0, **x**) e^(−iHt).

## 1. Condição KMS

Seja W(t, **x**) = ⟨φ(t, **x**) φ(0, **0**)⟩. Como φ(t + iβ) = e^(−βH) φ(t) e^(βH), a ciclicidade do traço dá

⟨φ(0) φ(t + iβ)⟩ = (1/Z) Tr[e^(−βH) φ(0) e^(−βH) φ(t) e^(βH)] = (1/Z) Tr[e^(−βH) φ(t) φ(0)] = ⟨φ(t) φ(0)⟩.

## 2. Ocupação dos modos

Com φ = ∫ d³k / [(2π)³ √(2k)] (a_k e^(−ikt + i**k·x**) + h.c.) e [a_k, a†_k′] = (2π)³ δ³(**k** − **k**′), a invariância por translação dá

W(t, **x**) = ∫ d³k / [(2π)³ 2k] [(1 + n_k) e^(−ikt + i**k·x**) + n_k e^(ikt − i**k·x**)].

A condição do passo 1 exige (1 + n_k) e^(−βk) = n_k, ou seja, n_k = 1/(e^(βk) − 1), a ocupação de Bose.

## 3. Soma de imagens

Expandindo 1 + n_k = Σ_(m≥0) e^(−mβk) e n_k = Σ_(m≥1) e^(−mβk), cada fator e^(−mβk) desloca t por ∓imβ. Com o vácuo W₀(t, r) = 1/[4π²(r² − (t − iε)²)],

W(t, r) = Σ_(m∈ℤ) 1/[4π²(r² − (t − imβ)²)].

## 4. Função de Hadamard

G^(1) = W(x, x′) + W(x′, x). Trocar os pontos troca t por −t, o que só reordena a soma em m. Portanto, para m ≠ 0,

G^(1)_reg = Σ_(m≠0) 1/[2π²(r² − (t + imβ)²)].

O termo m = 0 é G_M = 1/[2π²(r² − t²)], a eq. 9 de Li & Gott e, com σ = (r² − t²)/2, a eq. 17 de Kim & Thorne. Regularizar remove esse termo.

Os termos m e −m são complexos conjugados, então cada par é real. Na coincidência, cada par é função só de mβ e decai como m⁻⁴ depois de duas derivadas. A cauda de Kummer usa isso, e o diagnóstico `decay_spread` confere.

## 5. Conferências numéricas

- **Forma de modos:** do passo 2, integrando os ângulos, G^(1)_reg = (1/2π²) ∫₀^∞ dk 2n_k cos(kt) sin(kr)/r. A soma de imagens coincide com essa integral em pontos separados (`outputs/report.md`).
- **Coincidência:** a soma de imagens dá G^(1)_reg = T²/6, igual a (1/π²) ∫₀^∞ k n_k dk.
- **Referência:** ρ = ∫ d³k/(2π)³ k n_k = π²T⁴/30 (Stefan-Boltzmann, um grau de liberdade), conferida por integração numérica.

## O que esta derivação não decide

Ela fixa G^(1) e a referência ρ, sem usar nenhum dos dois artigos. Com ela, a comparação mede qual operador leva de G^(1) a ρ. Ela não diz por que a eq. 4a foi impressa sem o fator ½.
