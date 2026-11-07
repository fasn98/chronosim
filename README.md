# chronosim

Re-implementações independentes, com proveniência verificável, dos resultados publicados por trás da **Conjectura de Proteção da Cronologia** de Stephen Hawking (1992).

Segue a mesma disciplina do [optosim](https://github.com/fasn98/optosim): cada número carrega a sua origem, nada é inventado, e os limites são ditos na cara.

## A pergunta física

A relatividade geral permite espaços-tempo com **curvas fechadas tipo tempo** (caminhos que voltam ao próprio passado). Hawking conjecturou que a física quântica impede que eles se formem: perto do "horizonte da cronologia", a energia do vácuo cresceria sem limite e destruiria a máquina do tempo antes que ela funcionasse.

O espaço de **Misner** é o exemplo mais simples para testar isso. É o espaço de Minkowski em que os pontos são identificados por um boost de rapidez *b*. Na região R (x > |t|), em coordenadas de Rindler (ξ, η), a identificação η ~ η + nb cria curvas fechadas tipo tempo, e o horizonte da cronologia fica em ξ → 0.

## Fase 1 (esta versão)

Para um campo escalar conformemente acoplado, Li & Gott (1998) escrevem a função de Hadamard como uma **soma de imagens** (eq. 8), regularizam subtraindo o vácuo de Minkowski (eq. 10) e obtêm por point-splitting (eq. 11) uma **fórmula fechada** (eq. 12):

⟨T_ab⟩ = [(2π/b)⁴ − 1] / (1440 π² ξ⁴) · (4ξ² dη_a dη_b + g_ab)

O `chronosim` calcula a soma de imagens e o point-splitting **numericamente**, sem usar a eq. 12, e compara com ela. A fórmula contém os dois lados do debate:

- **o mecanismo de Hawking:** |ρ| ∝ ξ⁻⁴ diverge no horizonte;
- **o contra-exemplo de Li & Gott:** em b = 2π, ⟨T_ab⟩ = 0 em todo o espaço.

### Resultados (`python -m chronosim.run`)

| Teste | Resultado [origem] |
|---|---|
| Soma de imagens + point-splitting vs. eq. 12, 15 pontos (ξ ∈ {0,5; 1; 2}, b ∈ {1; 3; 2π; 10; 30}), 4 componentes | maior diferença **5,5e-16** relativa à escala 1/(1440π²ξ⁴) [COMPUTED] |
| Parte de Rindler (n = 0) vs. limite b → ∞ da eq. 12 | concordância melhor que 1e-15 [COMPUTED] |
| Traço (o campo é conforme, deve ser zero) e componentes fora da diagonal | < 1e-18 da escala [COMPUTED] |
| Inclinação de log\|ρ\| contra log ξ, ξ de 10⁻³ a 1 | **−4** até 15 dígitos [COMPUTED] |
| b = 2π | ρ/escala ~ 1e-19, zero dentro do erro numérico [COMPUTED] |

As figuras ficam em `outputs/`: a divergência no horizonte (`fig1_horizon_divergence.png`) e o anulamento em b = 2π (`fig2_self_consistent_b.png`).

## Proveniência por construção

Todo número de saída é um `Quantity` com uma origem:

| Origem | Significado |
|---|---|
| `CITED` | fórmula ou constante transcrita de uma fonte, com a equação indicada |
| `DERIVED` | avaliado exatamente (forma fechada, precisão arbitrária) a partir de fórmulas citadas |
| `COMPUTED` | resultado de um procedimento numérico do pacote (aproximação por construção) |
| `HEURISTIC` | valor escolhido por nós, sem fonte |

A origem só pode **piorar** numa computação. Não existe caminho de código que promova um valor `COMPUTED` a `DERIVED`. Nenhuma grandeza é `MEASURED`: este é um projeto puramente teórico e não há medição envolvida.

**Constantes `HEURISTIC`, reportadas em toda execução:** precisão de trabalho (80 dígitos), distância de point-splitting δ = 1e-10, passo de diferenças finitas h = 1e-11 e número de pares de imagens somados explicitamente (24).

**Guarda anti-fabricação:** um teste percorre a árvore sintática de todos os módulos e **falha o build** se algum deles puder alcançar um gerador de números aleatórios.

## Método numérico

- Precisão arbitrária com `mpmath`, para que o cancelamento entre a imagem n = 0 e o vácuo de Minkowski (ambos ~1/δ²) não perca dígitos.
- Derivadas por diferenças centrais. O limite X′ → X é feito separando os pontos em +δ e −δ ao longo de z e tirando a média, o que cancela os termos ímpares.
- γ é calculado como 2·asinh(√(u/2)), com u = cosh γ − 1 obtido diretamente. É exato e evita o cancelamento do acosh perto de 1.
- A soma de imagens usa 24 pares explícitos e a cauda pela aceleração de Kummer, supondo que os termos decaem como n⁻⁴. Essa hipótese é **verificada**, não assumida: o diagnóstico `decay_spread` mede quanto n⁴·termo(n) varia nos últimos termos, e um teste falha se a variação passar de 1e-14. A cauda corresponde a cerca de 2e-5 do total.

## O que isto estabelece e o que não estabelece

- ✅ **Estabelece:** a eq. 12 de Li & Gott decorre das eqs. 8 e 11 e se reproduz de forma independente até a precisão numérica. Isso inclui a divergência ξ⁻⁴ no horizonte e o anulamento em b = 2π.
- ❌ **Não estabelece:** que viagens no tempo sejam possíveis ou impossíveis. O cálculo é **semiclássico**: campo quântico num espaço-tempo fixo, sem reação gravitacional. Kay, Radzikowski & Wald (1997) provaram que, num horizonte de Cauchy gerado compactamente, a função de dois pontos não pode ser Hadamard. Ou seja, a própria teoria semiclássica deixa de valer ali, e a resposta depende de uma teoria de gravidade quântica.
- **Escopo:** um único campo (escalar sem massa, acoplamento conforme), em espaço plano, só na região R.

## Estado das fontes

| Fonte | Estado |
|---|---|
| Li & Gott (1998) | Transcrito da versão arXiv gr-qc/9711074. A comparação com a versão publicada (PRL 80, 2980) está **PENDENTE** por falta de acesso. |
| Kim & Thorne (1991) | Lido na versão publicada, em scan. Transcrições duplas em `docs/equation_transcription.md`. |
| Hawking (1992) | Lido; título, volume e página conferidos no próprio PDF. |
| Klinkhammer (1992), tese | Lida (eqs. 3.3, 3.7, 3.8), ver `docs/literature_check_factor2.md`. |
| Visser, hep-th/9303023 | Lido no arXiv (eqs. 10 e 34). Versão publicada não conferida. |
| Birrell & Davies (1982) | **PENDENTE:** a página sobre funções de Green térmicas não pôde ser conferida. A derivação usada está em `docs/thermal_images.md`. |
| Candelas & Deutsch (1977) | Volume e página confirmados só no Crossref. **PENDENTE:** o texto não foi lido (acesso fechado), e nenhuma equação dele foi conferida. |
| Hiscock & Konkowski (1982) | **PENDENTE:** o artigo não está disponível em versão aberta. A fase 4 fica bloqueada até o PDF ser obtido, e nenhuma fórmula dele é reconstruída de memória. |

PDFs de artigos ficam em `papers/`, que não entra no repositório.

## Instalação e uso

```bash
git clone https://github.com/fasn98/chronosim.git
cd chronosim
python3 -m venv .venv && source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -r requirements.txt
python -m pytest            # ~370 testes, ~2,5 min
python -m chronosim.run     # validação completa, figuras e relatório em outputs/ (~5 min)
```

## Roteiro

1. ✅ **Misner, vácuo de Li-Gott:** soma de imagens contra a forma fechada.
2. **Máquina do tempo de buraco de minhoca (Kim & Thorne, 1991):** o tensor ⟨T_ab⟩ perto do horizonte de Cauchy e a estimativa de onde a gravidade quântica passa a dominar. 🟡 Cálculo completo (F1 a F7, árbitros, figuras); falta a revisão do squad antes do merge. Resultados em `outputs/report.md`; plano em `docs/phase2_plan.md`.
3. **Bolas de bilhar (Echeverria, Klinkhammer & Thorne, 1991):** mecânica clássica em buracos de minhoca com curvas fechadas tipo tempo e busca de soluções autoconsistentes.
4. 🔒 **Misner, vácuo de Minkowski com imagens (Hiscock & Konkowski, 1982):** bloqueada até o PDF do artigo ser obtido.

## Referências

- Li, L.-X. & Gott, J. R. (1998). Self-consistent vacuum for Misner space and the chronology protection conjecture. *Phys. Rev. Lett.* 80, 2980. DOI 10.1103/PhysRevLett.80.2980 · arXiv:gr-qc/9711074
- Kim, S.-W. & Thorne, K. S. (1991). Do vacuum fluctuations prevent the creation of closed timelike curves? *Phys. Rev. D* 43, 3929. DOI 10.1103/PhysRevD.43.3929
- Echeverria, F., Klinkhammer, G. & Thorne, K. S. (1991). Billiard balls in wormhole spacetimes with closed timelike curves: Classical theory. *Phys. Rev. D* 44, 1077. DOI 10.1103/PhysRevD.44.1077
- Hiscock, W. A. & Konkowski, D. A. (1982). Quantum vacuum energy in Taub-NUT (Newman-Unti-Tamburino)-type cosmologies. *Phys. Rev. D* 26, 1225. DOI 10.1103/PhysRevD.26.1225 · **não lido; ver Estado das fontes**
- Hawking, S. W. (1992). Chronology protection conjecture. *Phys. Rev. D* 46, 603. DOI 10.1103/PhysRevD.46.603
- Kay, B. S., Radzikowski, M. J. & Wald, R. M. (1997). Quantum field theory on spacetimes with a compactly generated Cauchy horizon. *Commun. Math. Phys.* 183, 533. arXiv:gr-qc/9603012
- Misner, C. W. (1967). Taub-NUT space as a counterexample to almost anything. *Lectures in Applied Mathematics* 8, 160.

**Autor:** Fabio Augusto Sales Nogueira · ORCID [0009-0001-9303-8707](https://orcid.org/0009-0001-9303-8707)
