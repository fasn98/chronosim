# Fase 3 · plano: bolas de bilhar em buracos de minhoca (Echeverria, Klinkhammer & Thorne, 1991)

**Estado: aprovado, com as decisões abaixo. Ordem: reprodução R1–R11 contra as eqs. 4.2–4.3, com parada ao final; depois o mapa e a certificação. Reprodução R1–R11 calculada (`outputs/report.md`); parada para revisão do usuário. Mapa e certificação não iniciados.**

## Decisões do usuário (2026-11-07)

1. **Seção V (não coplanar):** fica como etapa opcional, depois do caso coplanar.
2. **Certificação intervalar:** sim, em duas etapas.
   - **Etapa 1:** teste de Krawczyk (`mpmath.iv`, arredondamento para fora) nas soluções das eqs. 4.2–4.3, num subconjunto de condições iniciais que inclui os casos "paradoxais" do artigo.
     - Cada solução certificada é uma prova de existência **para aquela condição inicial**.
     - Só o que o teste provar recebe o rótulo CERTIFIED.
     - A hipótese do teorema (função, caixa, critério de inclusão) fica documentada em `docs/`.
   - **Etapa 2:** só depois de validar a etapa 1, e com parada antes. Avaliar a viabilidade de certificar regiões de condições iniciais, por Krawczyk paramétrico ou por subdivisão adaptativa, e o teste de exclusão intervalar onde não houver solução.
   - **Origem CERTIFIED:** não existe hoje. O encaixe no sistema de origens será proposto, e `provenance.py` só muda depois de aprovado.
   - **No relatório:** uma prova sobre regiões não demonstra a afirmação para todas as condições iniciais, a menos que as regiões cubram o espaço inteiro.
3. **Grade:** aprovada como proposta. B, R/B, N_max e L_max são `heuristic()`. Teste de sensibilidade: as contagens não podem mudar de N_max = 12 para 16, nem de L_max = 50 D para 80 D. Se mudarem, relata-se onde.
4. **T_d:** substituído pelo teste de escala. O relatório deixa claro que a observação "só v₁·T_d/D importa" é deste projeto, não do artigo, e que é o teste que a confirma.
5. **Citações de lacunas:** nenhuma das 20 é usada antes da transcrição dupla. A citação da p. 1088 ("their high nonlinearity has prevented us from proving definitively that there always is a solution") fica **pendente de conferência humana** até o usuário conferir.

Fonte: Echeverria, F., Klinkhammer, G. & Thorne, K. S. (1991), *Phys. Rev. D* 44, 1077–1099, DOI 10.1103/PhysRevD.44.1077. Volume, número (4), data (15 de agosto de 1991) e páginas foram conferidos no PDF (`papers/ECHprd91.pdf`, 23 páginas, lido inteiro). O PDF é um scan sem texto extraível: tudo foi lido como imagem, a 110 dpi, e as eqs. 3.7, 3.8, 4.7 e 4.8 foram relidas a 250 dpi. A página k do PDF é a p. 1076 + k da revista. Abaixo, "p." é sempre a página da revista.

**Status das citações deste plano:** as citações literais abaixo vêm de **uma única leitura**. Elas entram na transcrição dupla da etapa B, com a segunda passada feita por um subagente sem acesso à primeira, antes de qualquer uso em código ou relatório. Os números de equação foram conferidos. As fórmulas não estão reescritas aqui, de propósito: elas só entram em `docs/equation_transcription.md` depois das duas passadas.

## O modelo

- **Espaço-tempo** (p. 1079, Fig. 1b): o "espaço-tempo de máquina do tempo eterna". É Minkowski plano, de onde se removem duas bolas em repouso, separadas por D, com defasagem temporal T_d entre elas. As superfícies das duas bolas são identificadas, então as bocas têm raio b e o buraco de minhoca tem comprimento desprezível. O caso de paradoxo dos gêmeos (Fig. 1a) é tratado como se fosse este, porque os dados iniciais ficam muito depois do retorno da boca 2 (p. 1079).
- **Identificação das bocas** (p. 1079, Fig. 2): pontos diametralmente opostos são identificados, ou seja, pontos ligados por reflexão no plano a meio caminho entre as bocas. A bola entra pela boca 2 e sai pela boca 1, voltando no tempo por T_d (p. 1081, p. 1083).
- **Unidades** (p. 1079): distâncias em D, tempos em T_d. B ≡ b/D é o raio da boca, R ≡ r/D é o raio da bola, e v é medida em D/T_d.
  - **Observação minha, não do artigo:** com essas unidades, numa mecânica não relativística, T_d não aparece como parâmetro separado. Mudar T_d com v₁ fixa equivale a mudar v₁. O eixo "defasagem temporal" do mapeamento é então o mesmo que o eixo v₁·T_d/D. Proponho verificar isso com um teste de escala, em vez de varrer T_d.
- **Bola** (p. 1078, p. 1081): é uma partícula clássica, não relativística, de esfera rígida, sem graus de liberdade internos. As colisões consigo mesma são perfeitamente elásticas. A bola é um objeto de teste, e o buraco de minhoca não recua (p. 1081).
- **Modelo da travessia** (p. 1081, itens i–iii):
  1. Não se consideram colisões em que o centro da bola está de um lado da garganta e a superfície que colide está do outro.
  2. Mesmo com R até B/2, o centro atravessa como o de uma bola infinitesimal.
  3. A bola recupera a forma esférica instantaneamente.

  O artigo diz: "(Our choice of B/2 rather than B/4 or 9B/10 is quite arbitrary.)" (p. 1081).
- **Regras de travessia** (eqs. 2.1a, 2.1b, p. 1082, Fig. 4): a velocidade em módulo é conservada, e a direção muda conforme a posição angular φ na boca.
- **Colisão única** (eqs. 2.2a–2.4b, p. 1082, Fig. 5): com uma só autocolisão, as leis de conservação somadas à regra de travessia admitem exatamente duas soluções. Uma é a "troca de velocidades" (2.3a, 2.3b). A outra é a "troca espelhada" (2.4a, 2.4b).

### Aproximações declaradas

| Aproximação | Onde vale | Página |
|---|---|---|
| Não relativístico: "travels with a speed small compared to light so special-relativistic effects can be ignored" | todo o artigo | p. 1078 |
| Paradoxo dos gêmeos tratado como máquina do tempo eterna | todo o artigo | p. 1079 |
| Bola como objeto de teste, sem recuo do buraco de minhoca | todo o artigo | p. 1081 |
| Uma só autocolisão | seção II.B em diante | p. 1080, p. 1082 |
| R ≪ B ≪ D ≡ 1 | seções III e IV.A | p. 1080, eq. 3.3 |
| Bola de tamanho infinitesimal "for ease of presentation" | resto da seção III | p. 1083 |
| Ordem dominante em B, R/B, 1/L, h/L | seção III.B | p. 1087 |
| B < ½ (bocas não se sobrepõem) e R/B < ½ (marés desprezíveis) | seção IV.B, apêndice A | p. 1088 |
| B ≪ 1, R/B < ½, uma travessia, colisão antes da travessia | seção IV.C, apêndice B | p. 1090 |
| Expansão até ordem quadrática em β e ε = α − 2φ | apêndice B | p. 1097 |
| Quase coplanar, z₁ ≪ BR (desacoplamento) | seção V | p. 1094, eq. 5.12 |

## O que o artigo afirma

| # | Afirmação | Onde |
|---|---|---|
| A1 | Regras de travessia coplanar | eqs. 2.1a, 2.1b, p. 1082 |
| A2 | Uma colisão única é troca de velocidades ou troca espelhada | eqs. 2.2–2.4, p. 1082 |
| A3 | Toda trajetória coplanar apontada para entre as bocas tem multiplicidade ≥ 2. A segunda solução tem uma travessia e troca de velocidades, com v₁′ = (D − 2B)/T_d ≡ 1 − 2B, qualquer que seja v₁ | p. 1083, Fig. 6 |
| A4 | Trajetória perpendicular à linha dos centros, pelo ponto médio: solução com uma travessia e troca espelhada, tan ψ = v₁ (B ≪ 1) | eq. 3.1, p. 1083, Fig. 7b |
| A5 | No limite v₁ → 0 dessa solução, "this solution is really a continuous infinity of solutions: the time T of the collision is completely arbitrary" | p. 1084 |
| A6 | Solução com duas travessias: sin ψ + tan ψ = 2v₁, com os instantes T = −1/(1+cos ψ), −cos ψ/(1+cos ψ), 0, cos ψ/(1+cos ψ), 1/(1+cos ψ) | eq. 3.2, p. 1084, Fig. 7c |
| A7 | Existe solução com três travessias (Fig. 7d). O ângulo fica "as an exercise for the reader", sem valor impresso | p. 1084 |
| A8 | Para R ≪ B ≪ 1, v₁ > 1 e ψ_A > B, há uma solução para cada n > n_min de travessias, ou seja, multiplicidade infinita | eqs. 3.3–3.10, pp. 1085–1087, Figs. 8 e 9 |
| A9 | Bola inicialmente em repouso, fora do elipsoide r_min: só a solução trivial. A prova tem uma ressalva (ver lacunas) | p. 1087–1088, Fig. 10 |
| A10 | Bola em repouso dentro de r_min, sobre a linha dos centros: existe solução não trivial | p. 1088 |
| A11 | Para R ≪ B ≪ 1, toda trajetória coplanar perigosa satisfaz ψ_A > B e v₁ > 1, e portanto tem multiplicidade infinita | eqs. 4.1a, 4.1b, p. 1088 |
| A12 | Para B < ½ e R/B < ½, uma colisão de troca espelhada antes de uma travessia: as soluções autoconsistentes correspondem às soluções não espúrias de duas equações em (α, β). As espúrias são listadas | eqs. 4.2–4.6c, p. 1089; derivação nas eqs. A1–A19, pp. 1094–1097 |
| A13 | Triângulo perigoso, para B ≪ 1 e R/B < ½ | eqs. 4.7, 4.8, Fig. 12, p. 1090 |
| A14 | Na região central do triângulo há duas soluções analíticas, de classe I e de classe II, lineares em λ₁ | eqs. 4.10–4.13c, Fig. 13, p. 1091; eqs. B16–B20 |
| A15 | Na borda superior esquerda, β(φ) muda bruscamente em φ = π/4 e em φ = √B. Há solução numérica até o canto. Os parâmetros são R = 10⁻⁵, B = 10⁻², λ₁ = ½ | eqs. 4.14, 4.15, Fig. 14, p. 1092; eqs. B29–B35 |
| A16 | Na região da borda superior, com s = −1, há uma faixa de φ sem solução de classe II não espúria, perto de φ ~ π/4 | p. 1099, após a eq. B31; p. 1091 |
| A17 | Trajetórias quase coplanares: perigosas só se z₁ < B sin[½ arctan(2R/L₂ cos θ)]. Com z₁ ≪ BR, cada solução coplanar tem uma correspondente não coplanar, e vice-versa | eqs. 5.1–5.12, pp. 1092–1094 |
| A18 | Conclusão: o problema de Cauchy é mal posto por excesso de soluções. Não se viu evidência de multiplicidade zero | p. 1094, seção VI |

## O que o artigo declara como não demonstrado ou em aberto

Só entram aqui os trechos em que o texto declara a lacuna. Citações literais, de leitura única (ver status acima).

1. p. 1079: "However, we have searched hard for initial trajectories with zero multiplicity and have found none. On the other hand, our search has not covered all initial trajectories (see especially Sec. V), so we cannot guarantee the nonexistence of zero-multiplicity trajectories."
2. p. 1079: "However, we have not found a manageable way to analyze the case of coplanarity to within a distance ~ B."
3. p. 1080: "There presumably are also multiple-collision solutions, and we speculate about some possible, rather strange ones in the paragraph containing Eq. (3.11)."
4. p. 1080: "We suspect, but have not proved, that the (finite-measure) initial trajectories with speeds v₁ ≪ D/T_d ≡ 1 and with impact parameters h ≫ D ≡ 1 also have unit multiplicity; see Sec. III."
5. p. 1081: "Our spot checks in these extreme regions have not turned up any initial trajectories for which numerical solutions do not exist."
6. p. 1081: "Just how small R must be for this depends on one's model for the internal structure of the ball."
7. p. 1083: "We also argue, but do not prove firmly, that there is only a single solution for any ball with (i) an initial speed that is sufficiently slow but not zero, and (ii) an initial path of motion that, if extended forever, remains far from the wormhole."
8. p. 1087: "(Readers who seek higher rigor than we do might worry that our analysis has examined only the leading-order effects in the small parameters B, R/B, 1/L, and h/L and has not proved rigorously that higher-order corrections are negligible. We are not worried.)"
9. p. 1087, sobre soluções com várias colisões perto do buraco de minhoca: "We have searched cursorily for such peculiar solutions, without success, and we suspect they do not exist. However, we have no proof."
10. p. 1088, sobre uma sequência infinita de autocolisões: "We have not been able to rule out such a solution rigorously, but it seems exceedingly unlikely that one could exist."
11. p. 1088: "Note that this proof fails if the ball is initially at rest inside radius r_min".
12. **[pendente de conferência humana]** p. 1088, sobre as eqs. 4.2–4.3: "but their high nonlinearity has prevented us from proving definitively that there always is a solution." E logo depois: "we shall argue that this is probably so for all initial trajectories, i.e., the multiplicity is probably always nonzero."
13. p. 1089: "It is far from obvious, just looking at the equations, whether there exist values of the wormhole and ball radii B, R and initial-trajectory parameters v₁, θ, φ that produce zero solutions. Numerical exploration, and the analytic considerations of the next section, have not turned up any such zero-multiplicity trajectories."
14. p. 1090, legenda da Fig. 12: "Analytic solutions cannot be derived by the techniques of Appendix B for the white regions of the triangle (left corner and lower left edge), but numerical solutions have been found in spot checks throughout that white region."
15. p. 1091: "At all points near the left corner and lower left edge where we have tried, we found numerical solutions. Thus, it seems likely that solutions exist everywhere in the dangerous triangle."
16. p. 1091: "We suspect, but have not proved, that the missing class-II solution actually exists, but the ball first encounters its collision shortly after passing through the wormhole, rather than before, and therefore this solution is beyond the domain of validity of our analysis."
17. p. 1092: "We find no evidence, when B ≪ 1 and R/B < ½, for initial trajectories with zero multiplicity."
18. p. 1092: "Thereby we shall accumulate evidence which suggests, but does not really prove, that all noncoplanar initial trajectories have multiplicity larger than zero". E também: "We have not found a good way to analyze dangerous initial trajectories near the boundary of the region (5.1)."
19. p. 1094: "Thus, we can say nothing about the existence of solutions over the full range." E também: "However, these conclusions make us suspect that even when z₁ ~ BR, all initial trajectories will have at least one self-consistent solution."
20. p. 1094, seção VI: "On the other hand, we have seen no evidence for a stronger type of ill posedness: generic initial trajectories with no self-consistent solutions."

## Proposta 1 · reprodução independente dos resultados do artigo

O cálculo independente é um **solucionador geométrico exato**, com B e R finitos, sem as aproximações do artigo. Ele usa só o modelo:
- retas entre eventos;
- a identificação por reflexão da Fig. 2;
- o atraso T_d;
- a colisão elástica de esferas iguais, na forma geral (normal de contato e transferência de momento ao longo dela).

Ele **não** usa as eqs. 2.3/2.4 nem 4.2/4.3. Estas ficam como alvo.

| # | Alvo | Cálculo independente | Teste de concordância |
|---|---|---|---|
| R1 | Eq. 2.1b (A1) | Travessia pela identificação da Fig. 2 | Igualdade exata em uma grade de φ e ψ_in. Isto é uma derivação usada como checagem da transcrição |
| R2 | Eqs. 2.3, 2.4 (A2) | Colisão elástica geral com v₁′ = v₂ e v₂′ = v₁ impostos | As únicas soluções numéricas são as duas trocas. Também se confere r₁ − r₂ das eqs. 2.3b e 2.4b |
| R3 | v₁′ = 1 − 2B (A3) | Solução da Fig. 6b com R = 0, e depois com R finito (Fig. 6c) | Exata para R = 0. Desvio medido em função de R |
| R4 | Eq. 3.1 e o contínuo em v₁ → 0 (A4, A5) | Solução da Fig. 7b com B finito | Razão contra tan ψ = v₁ tende a 1 com B → 0, com a ordem medida. O contínuo é testado com v₁ = 0 |
| R5 | Eq. 3.2 e os cinco instantes (A6) | Solução da Fig. 7c com B finito | O mesmo de R4 |
| R6 | Fig. 7d (A7) | Solução com três travessias | Sem valor impresso. Sai como COMPUTED, sem alvo |
| R7 | Eqs. 3.6–3.10 (A8) | Soluções com n travessias, n até N_max, para (v₁, ψ_A, h) fixos | L(n) contra a eq. 3.7. Menor n encontrado contra a eq. 3.8. O mapa de lente divergente (3.9a, 3.9b) é conferido travessia por travessia |
| R8 | Eqs. 4.2, 4.3, 4.6 (A12) | O solucionador geométrico, comparado com as raízes das eqs. 4.2 e 4.3 transcritas | As mesmas soluções (α, β) nos dois caminhos. Este é o árbitro principal da transcrição, como o árbitro térmico na fase 2 |
| R9 | Eqs. 4.12, 4.13 e Fig. 13 (A14) | Solucionador geométrico na região central do triângulo | Desvio contra as formas analíticas, em função de B e R. A ordem do erro é medida |
| R10 | Eqs. 4.14, 4.15, B32 e Fig. 14 (A15) | O mesmo, com R = 10⁻⁵, B = 10⁻², λ₁ = ½ | Curva β(φ) contra as eqs. 4.14, 4.15 e B32, em cada faixa |
| R11 | Classe II ausente perto de φ ~ π/4 (A16) | Busca de classe II com a colisão **depois** da travessia, que o solucionador geral permite | Pergunta aberta do próprio artigo (lacuna 16). Ver Q3 |

**Fora do escopo, a menos que você peça:**
- **A9 e A10** (prova com elipsoides): é uma prova, e não há o que reproduzir numericamente. O mapa dá evidência no caso v₁ = 0 (Q5), sem tocar na prova.
- **Seção V (não coplanar):** exige um solucionador 3D e seis parâmetros. Proponho deixar para uma etapa opcional depois do mapa coplanar, com no máximo a checagem do critério 5.1.

## Proposta 2 · mapeamento sistemático

**Espaço de condições iniciais (coplanar):**
- a trajetória inicial é uma reta no plano, dada por (h, ψ_A) da Fig. 9a, ou por (θ, φ) da Fig. 11, com h = −B sin θ e ψ_A = θ + φ (p. 1089);
- a velocidade v₁ completa o terço.

**Parâmetros:** B ∈ (0, ½) e R/B ∈ (0, ½), que é o domínio de validade declarado na p. 1088.

**Esboço da grade** (a densidade final será fixada pelo critério de estabilização abaixo):
- **B:** {0,01; 0,03; 0,1; 0,2; 0,3; 0,45}.
- **R/B:** {0,01; 0,1; 0,25; 0,45}.
- **(h, ψ_A):** grade densa sobre a faixa perigosa e uma margem não perigosa, por exemplo 64 × 64 por par (B, R).
- **v₁:** saída contínua (ver método), amostrada em uma grade log de 10⁻³ a 10².

**Classificação de cada condição inicial:**
- **0:** "nenhuma solução encontrada sob [critério]";
- **1;**
- **k finito;**
- **≥ k crescente com N_max:** uma solução nova para cada n até N_max, o que é compatível com multiplicidade infinita mas não a demonstra;
- **contínuo:** família a um parâmetro com resíduo nulo ao longo de um intervalo.

O critério de cada classe será registrado como `heuristic()`.

**O que o número conta:** soluções com **uma colisão** e no máximo N_max travessias, com a colisão dentro de uma janela |L| ≤ L_max. A trajetória livre conta quando é consistente, ou seja, quando não é perigosa. Soluções com várias colisões ficam fora do modelo do artigo e fora da contagem, então a multiplicidade do mapa é um **limite inferior** dentro desse modelo. Essa ressalva vai em cada figura e tabela.

## Proposta 3 · método numérico: achar todas as soluções, de forma verificável

### Formulação por tiro com 1D aninhado

Ela generaliza o método de tentativa e erro da Fig. 8 (p. 1084–1085) para B e R finitos e para colisão geral (troca de velocidades ou espelhada). Fixada a reta inicial:

1. **Parâmetro externo s:** posição do centro da bola jovem no evento de colisão, ao longo da reta inicial.
2. **Parâmetro do meio ν:** ângulo da normal de contato. O centro da bola velha fica em Q = P − 2R n(ν).
3. **Parâmetro interno χ:** direção de saída da bola jovem a partir de P. Para cada itinerário simbólico (a sequência de bocas atravessadas, até N_max travessias), propaga-se a reta com as regras de travessia. A condição é que a trajetória de retorno passe por Q: distância transversal com sinal, achada por **bracketing** em χ.
4. **Velocidade:** a rapidez em módulo se conserva nas travessias, então ela sai explicitamente da condição de tempo, comprimento do laço / rapidez = n·T_d.
5. **Lei de colisão:** de v₁′ e v₂, a lei dá uma condição escalar C(s, ν) = 0, achada por **bracketing** em ν, e dá v₁ como **saída**.
6. **Multiplicidade:** para cada ramo, v₁(s) é uma função de s. A multiplicidade em v₁* é o número de cruzamentos v₁(s) = v₁*, também por bracketing. Uma única varredura por reta inicial dá então a multiplicidade **para todas as velocidades de uma vez**.

Todas as buscas são 1D, por mudança de sinal em grade determinística, refinadas por bisseção ou Brent e polidas com `mp.findroot` em precisão estendida. O solucionador também inclui itinerários em que a bola atravessa antes da colisão, para R11 e Q3.

### Verificação de cada solução, independente do solucionador

Um verificador reconstrói a história completa, com todas as encarnações e todos os segmentos, e confere:
- retas com rapidez constante;
- a regra de travessia em cada boca;
- a lei elástica na colisão, com tolerância registrada;
- que não há **nenhuma outra** sobreposição bola-bola em instante nenhum, nem sobreposição de um segmento com o interior de uma boca fora das travessias;
- a exclusão (i) da p. 1081.

Um candidato que falha é descartado e registrado com o motivo.

### Garantia contra solução perdida

Ela é numérica, não rigorosa.
- **Pontos fracos do bracketing:** ele perde raízes duplas (tangências) e pares de raízes mais próximos que o passo da grade.
- **Critério de estabilização:** a grade de cada nível é refinada (×2) até que a contagem fique igual em duas refinações seguidas, e a posição de cada raiz não mude além de uma tolerância. O fator, o número de refinações e as tolerâncias são `heuristic()` e aparecem em toda execução.
- **Detecção de tangência:** registra-se o mínimo de |resíduo| entre raízes. Um mínimo pequeno sem mudança de sinal sinaliza possível raiz dupla e dispara refinamento local.
- **Opcional, se você quiser algo mais forte num subconjunto:** certificação por aritmética intervalar (`mpmath.iv`, teste de Krawczyk) nas eqs. 4.2–4.3 transcritas. Isso **exclui** raízes numa caixa e **certifica** raízes únicas. Mesmo assim vale só para aquelas equações e para aquele recorte, não para o problema geral.

### Precisão e proveniência

- A varredura usa `float`, e o polimento e a verificação usam mpmath. Ambas as precisões são registradas.
- Todo resultado do mapa é **COMPUTED**, no máximo. Formas fechadas do artigo, avaliadas como referência, podem ser DERIVED, mas nenhum número do mapa.
- Nada usa gerador aleatório: as grades são determinísticas, e a guarda anti-fabricação continua valendo.

### Transcrição

As equações que entram em código passam por duas leituras, a segunda por subagente sem acesso à primeira, e são registradas em `docs/equation_transcription.md`:
- 2.1–2.4;
- 3.1, 3.2 e os instantes da p. 1084;
- 3.6–3.10;
- 4.1, 4.2–4.8 e 4.10–4.15;
- A13–A19 como checagem da 4.2–4.3;
- B11, B16–B20 e B29–B35;
- 5.1 só se a seção V entrar.

**Análise dimensional:** restaurando D e T_d em cada equação, testa-se a homogeneidade sob D → λD e T_d → μT_d.

## Perguntas em que o mapeamento pode dar evidência

**Evidência numérica não é prova.** Um mapa em grade finita, com N_max e L_max finitos, mostra o que foi encontrado dentro de um critério. Ele não mostra que algo não existe fora dele, nem demonstra uma afirmação universal.

- **Q1.** Fora de R ≪ B ≪ 1 (por exemplo B = 0,3, R/B = 0,45), existe trajetória coplanar perigosa sem solução de uma colisão? É a lacuna 12–13. O resultado seria relatado como "nenhuma solução encontrada sob [critério]", nunca como "paradoxo".
- **Q2.** No canto esquerdo e na borda inferior esquerda do triângulo perigoso, onde o artigo fez só *spot checks* (lacunas 14–15), há solução em toda a grade?
- **Q3.** Perto de φ ~ π/4, na borda superior, a classe II ausente aparece quando se permite a colisão depois da travessia (lacuna 16)?
- **Q4.** Trajetórias lentas e longe do buraco de minhoca (v₁ ≪ 1, h ≫ 1) têm uma única solução de uma colisão, dentro do critério (lacunas 4 e 7)?
- **Q5.** Com v₁ = 0, fora e dentro de r_min, o mapa encontra só a solução trivial, ou também a de A10? Isso dá evidência numérica, não a prova da p. 1088.
- **Q6.** Para v₁ > 1 e ψ_A > B, aparece uma solução para cada n do menor encontrado até N_max? Isso continua valendo com B e R/B finitos, fora do regime da eq. 3.3? É a lacuna 8, sobre correções de ordem superior.
- **Q7.** Como a contagem (até N_max) depende de B e R/B? O artigo não tem esse dado, e o resultado seria só descritivo.

## O que esta fase estabelece e o que não estabelece

- **Estabelece, se os testes passarem:** reprodução independente, com B e R finitos, das soluções que o artigo constrói (A1–A8, A12–A16), e um mapa de contagem de soluções de uma colisão, com critério de busca declarado.
- **A fase é clássica.** Ela testa a alternativa de autoconsistência de Novikov: sempre existir ao menos uma solução autoconsistente. Ela **não** testa o mecanismo quântico de Hawking (a divergência de ⟨T_μν⟩ das fases 1 e 2), e não diz nada sobre ele. O artigo remete a mecânica quântica ao "paper II" (p. 1081), que não foi lido.
- **Não estabelece:**
  - nenhuma afirmação universal, nem "toda trajetória tem solução", nem "existe trajetória sem solução";
  - nada sobre soluções com várias colisões;
  - nada fora do modelo de bola da p. 1081;
  - nada sobre o caso relativístico, sobre o paradoxo dos gêmeos com bocas em movimento, nem sobre trajetórias não coplanares (a menos que a etapa opcional entre).

## Decisões tomadas na etapa A (respondidas acima)

1. **Seção V (não coplanar):** fica fora, opcional depois, ou entra com escopo mínimo (checagem do critério 5.1)?
2. **Certificação intervalar** (Krawczyk) num subconjunto: quer, ou o critério de estabilização basta?
3. **Grade de parâmetros:** os valores de B e R/B acima servem? E os limites N_max (proponho 12) e L_max (proponho 50 D)?
4. **T_d como eixo:** concorda em trocá-lo pelo teste de escala, já que T_d e v₁ não são independentes neste modelo?
