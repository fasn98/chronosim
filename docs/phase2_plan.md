# Fase 2 · plano: Kim & Thorne (1991)

**Estado: aguardando aprovação. Nada foi implementado.**

Fonte: Kim, S.-W. & Thorne, K. S. (1991), *Phys. Rev. D* 43, 3929–3947, DOI 10.1103/PhysRevD.43.3929. Volume, data e páginas conferidos no PDF (`papers/KIMprd91.pdf`, 19 páginas, lido inteiro). O PDF é um scan sem texto extraível: as equações foram lidas como imagem, e as eqs. 14, 51, 53, 54, 67, 68 e 69 foram conferidas de novo com zoom.

Unidades do artigo: c = G = ħ = 1, ou seja, l_P = m_P = 1 (seção II.E). Atenção à notação: o ξ deste artigo é o inverso do desvio Doppler (eq. 30), não a coordenada de Rindler da fase 1.

## O sistema

Um buraco de minhoca com garganta infinitesimal, em espaço plano. A boca 1 está em repouso, e a boca 2 se move em direção a ela com velocidade β (eqs. 25–28). O raio da garganta é b, e D é a distância entre as bocas quando o horizonte de Cauchy se forma, com b ≪ D. O campo é escalar, sem massa, conformemente acoplado. O resultado vale em aproximação de óptica geométrica, e só a parte dominante da divergência.

## Fórmulas escolhidas

| # | O que o artigo afirma | Equações | Cálculo independente | Teste de concordância |
|---|---|---|---|---|
| F1 | Geometria do espaço de recobrimento e intervalo geodésico σ_N | 25–30, 33, 34–38 | Aplicar a identificação das eqs. 25 e 28 às linhas de mundo e calcular σ_N = ΔU·ΔV de forma exata no recobrimento, sem a aproximação da eq. 35 | Eq. 33 exata até o erro de arredondamento. Eq. 38 contra o σ_N exato: o erro relativo deve cair com a proximidade do horizonte, e a lei é medida |
| F2 | Desfocalização por óptica geométrica: Δ_N^(1/2) | 47, 48, 64 | Traçado de raios com matrizes ABCD: lentes divergentes de distância focal b/2 separadas por D | Razão contra a eq. 48 e a eq. 64 tende a 1 quando b/D → 0, com a ordem da correção medida |
| F3 | ⟨T^μν⟩ entre as bocas | Entrada: 61, 66, 4. Alvo: 67 | Aplicar `pointsplit.stress_tensor` à eq. 66 por diferenças finitas, sem usar a eq. 67 | Razão contra a eq. 67 tende a 1 perto de H_N. O desvio é medido em função de (X−T)/X e ρ/X, que o artigo despreza |
| F4 | ⟨T^μν⟩ na garganta | Entrada: 46, 50, 4. Alvo: 51 | O mesmo de F3, agora com as componentes em k e ρ̂ | Os três coeficientes da eq. 51 (termos em kk, kρ̂, ρ̂ρ̂), em um grid de θ < π/2 |
| F5 | Localização das hipersuperfícies polarizadas | 43, 44a, 62, 63 | Raiz numérica de λ_N = 0 (eqs. 46 e 61) com os pontos colapsados | Raiz contra as eqs. 43 e 62, e limite N → ∞ contra as eqs. 44a e 63 |
| F6 | Quadrimomento transportado | 53, 54, 68, 69 | Integração numérica (`mp.quad`) de F3 e F4, nos domínios escritos no artigo | Eqs. 54, 68 e 69 na razão direta. Eq. 53 só no coeficiente de ln(b/a), porque o termo finito fica fora da fórmula |
| F7 | Escalas físicas: onde cada conjectura corta a divergência | 9c, 12–14, 15b, e a seção II.F | Avaliar as fórmulas com D = 1 m, e resolver δg^VP = 1 para Δt | Aritmética contra os números impressos: l_P/D ~ 10⁻³⁵, Δt ~ 10⁻³⁵ l_P na condição de Hawking e T⁰⁰_max (ver ponto de atenção 2) |

**F7 é mais fraca que as outras.** As eqs. 9c e 12–14 são estimativas de ordem de grandeza (`~`), então só há aritmética a conferir. Isso não testa se a estimativa está correta. Os valores saem marcados como ordem de grandeza.

## Fora do escopo, e por quê

- **Eqs. 10 e 11 (flutuações de gravidade quântica):** a eq. 10 é uma hipótese vinda de outra referência, com desigualdades, e não dá para reproduzi-la. Ela entra só como entrada citada em F7.
- **Eqs. 6–8:** estimativas genéricas, cobertas de modo mais forte por F3 e F4.
- **Eqs. 17–24 (forma de Hadamard):** F2 reproduz o resultado da eq. 24 para este espaço-tempo, mas a validade geral do ponto de partida (eq. 20) não é testada. O próprio artigo aponta essa hipótese como possível falha, na seção VII.
- **Eq. 44b e caústicas:** não entram em nenhum alvo escolhido.
- **Fórmulas de Misner citadas na seção V.A.1:** o artigo diz que usa uma versão "ligeiramente alterada" do cálculo de Hiscock & Konkowski. Usamos só o que está impresso aqui, sem tomar isso como a fórmula deles.

## Pontos de atenção encontrados na leitura

Estes pontos vêm de uma conferência feita **à mão**, ainda sem execução nenhuma.

### 1. Fator 2 no operador (decisão necessária)

- A eq. 4a de Kim & Thorne é T_μν = lim D_μν G_reg, **sem** fator ½.
- A eq. 11 de Li & Gott, como transcrita na fase 1, tem o fator ½.
- As duas funções de Hadamard têm a mesma normalização: 1/(4π²σ) aqui (eq. 17), com σ igual à metade do quadrado da distância, e 1/(2π²·distância²) em Li & Gott.
- À mão, as eqs. 51, 53, 54, 67, 68 e 69 saem do operador **sem** o ½. Elas são coerentes com a eq. 4a, mas ficam 2 vezes maiores que o mesmo cálculo com o operador da fase 1.
- A fase 1 reproduz a eq. 12 de Li & Gott com o ½. Isso sugere que a diferença está na eq. 4a, mas não decide a questão: a eq. 1 define G^(1) com os dois ordenamentos, e não fizemos essa dedução de forma independente aqui.

**Proposta:** comparar com o operador da eq. 4a como impresso, porque é o cálculo deles, e relatar à parte a razão contra `stress_tensor` da fase 1. O código da fase 1 não muda. Para as conclusões do artigo, que são de ordem de grandeza, um fator 2 não faz diferença, mas o relatório diz isso às claras.

### 2. Valor numérico depois da eq. 14

O texto diz que T⁰⁰_max é l_P/D ~ 10⁻³⁵ da densidade de Planck, "i.e., ... ~ 10⁴⁹ g/cm³". Com os valores do próprio artigo (l_P ~ 10⁻³³ cm, m_P ~ 10⁻⁵ g), m_P/l_P³ ~ 10⁹⁴ g/cm³, o que daria ~ 10⁵⁹ g/cm³. A diferença é de 10¹⁰, e pode ser erro tipográfico. **Proposta:** relatar os dois valores sem corrigir nenhum.

### 3. Prefator da eq. 66

A eq. 66 usa X na forma simétrica em x e x′ só no λ, e o prefator D/[4π²X(D−X)] depende de X. Derivar esse prefator produz termos que o artigo despreza. **Proposta:** ler literalmente, com X do primeiro argumento, e medir o efeito pelo desvio em F3.

### 4. Citação de Hawking (1992)

A contestação de Hawking que li é a descrita na seção II.F, que cita um seminário de 1990 (ref. 13), anterior ao artigo de 1992. **Hawking (1992) em si não foi lido.** A seção "o que estabelece" pode citar a conclusão pela via de Kim & Thorne e dizer isso. Para citar direto o artigo de 1992, preciso do PDF.

### 5. Referência interna no texto

Logo após a eq. 52, o texto diz "the stress-energy tensor (52)", mas o tensor é a eq. 51. É um detalhe, sem efeito sobre as contas.

## Desenho da implementação

- **`chronosim/kim_thorne.py`:** fórmulas transcritas, cada uma com a equação citada (`DERIVED`). Avaliações de ordem de grandeza ficam marcadas no `source`.
- **`chronosim/wormhole.py`:** recobrimento, matrizes ABCD, raízes e integrais (`COMPUTED`). Usa `pointsplit.stress_tensor` sem alterá-lo. Os índices são convertidos de forma explícita, porque ele devolve T_ab e o artigo escreve T^μν.
- **`sources.py`:** adiciona `KIM_THORNE_1991`.
- **`run.py`:** tabelas, CSV e relatório. **Figuras:** T⁰⁰ contra Δt em log-log, com a camada de dipolo do artigo (fig. 3a), e δg^VP contra Δt com os dois cortes, em Δt ~ l_P (Kim & Thorne) e DΔt ~ l_P² (Hawking).
- **`heuristic()` previstas:** parâmetros de avaliação (β, b/D, X/D, ρ/X, θ), corte N da série, passos δ e h adequados à escala de Δt, tolerância de `mp.quad` e cutoff angular de F6. Valores físicos do artigo (D = 1 m, l_P, m_P) entram como citados.
- **Dependências novas:** nenhuma.
- **Testes:** tolerâncias definidas depois de medir o erro, com o valor medido ao lado, como nas fases anteriores. Os critérios acima são a previsão do que deve valer. Se algum falhar, relato a diferença e paro, sem ajustar nada para concordar.

## O que a fase vai estabelecer e o que não vai

- **Estabelece, se os testes passarem:** as fórmulas fechadas das eqs. 51, 54, 67, 68 e 69 decorrem das expressões de partida do próprio artigo (eqs. 46, 50, 61, 66), e a geometria e a desfocalização se reproduzem de forma independente.
- **Não estabelece:** que a função de Hadamard (eq. 20) esteja certa, porque ela é entrada do cálculo. Também não estabelece nenhuma das conjecturas sobre o corte da divergência.
- **Disputa com Hawking:** Hawking contestou a conclusão de Kim & Thorne. O corte em Δt ~ l_P depende do referencial, e ele propõe DΔt ~ l_P², em que δg^VP chega a 1. Esta fase reproduz o cálculo deles, e **não decide quem estava certo**. A figura mostra os dois cortes lado a lado, sem escolher um.
- **Escopo:** um único campo, em espaço-tempo fixo (cálculo semiclássico), para um buraco de minhoca idealizado com b ≪ D.

## Decisões para aprovação

1. **Fator 2:** usar a eq. 4a como impressa e relatar a razão contra a fase 1, como proposto no ponto 1?
2. **Valor 10⁴⁹:** relatar sem corrigir, como proposto no ponto 2?
3. **Escopo:** implementar de F1 a F7, ou cortar algum item?
4. **Hawking (1992):** citar pela via de Kim & Thorne, ou você vai obter o PDF de 1992?
