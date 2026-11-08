# chronosim · fase 3, reprodução completa (`python -m chronosim.run --full`)

## Fase 3 · Bolas de bilhar (Echeverria, Klinkhammer & Thorne, 1991): reprodução R1–R11

Todas as grandezas desta seção são **COMPUTED**, exceto a exclusão intervalar do caso vazio, quando ela aparece marcada como CERTIFIED (R8, caso vazio).

- **Modelo independente** (`chronosim/billiards.py`): o modelo exato que o artigo descreve em palavras e figuras.
  - Bocas de raio B identificadas por reflexão no plano x = ½ (Fig. 2).
  - Defasagem T_d: quem entra na boca 2 sai da boca 1 em T − 1.
  - Travessia em linha reta vista da garganta.
  - Colisão elástica de esferas iguais, na forma geral.
  - Não usa nenhuma equação do artigo.
- **Equações do artigo** (`chronosim/ekt.py`): transcritas duas vezes; a segunda passada foi cega. Registro em `docs/equation_transcription.md`.
- **Unidades:** D = T_d = 1.
- **Busca:** bracketing unidimensional aninhado, em grades determinísticas, com polimento em mpmath (40 dígitos).
  - Cada solução passa por um verificador que reconstrói a história inteira e confere as leis do modelo, incluindo a ausência de qualquer outra sobreposição entre encarnações.
  - "Intervalos não resolvidos" são tudo o que o bracketing não conseguiu seguir: trechos de s ou de ν em que o número de raízes mudou e não pôde ser ligado após 12 bisseções, bisseções (em s ou em ν) que perderam o ramo, e colchetes da busca da mira descartados (f indefinida dentro da bisseção, ou troca de sinal rejeitada como salto). Eles são contados, nunca descartados em silêncio, com uma exceção: intervalos da grade da mira com uma ponta indefinida (borda do domínio do laço) não são colchetes e não são registrados.
  - Cruzamentos rejeitados (v₁(s) passando pelo alvo num polo, ou impulso normal ≤ 0) e candidatos recusados pelo verificador entram como "candidatos rejeitados", com o motivo.

**O que estas contagens são.** São soluções com **uma colisão** e o itinerário indicado, encontradas sob as grades `ekt_*` (HEURISTIC). Não são a multiplicidade do problema, e não estabelecem nenhuma afirmação universal.

- **R1, eq. 2.1b.** O modelo de identificação reproduz ψ_out = ψ_in − 2φ com desvio máximo de 1,38e-40 rad numa grade de 25 × 25 em (φ, ψ_in). É uma derivação da regra, usada como checagem da transcrição.
- **R2, eqs. 2.2–2.4.** Foram testadas 144 colisões elásticas (pares v₁, v₂), impondo a troca de módulos das eqs. 2.2a–b.
  - Toda normal de contato que satisfaz a condição é uma das duas da eq. 2.3 ou 2.4: troca de velocidades 144 vezes, troca espelhada 144 vezes, outra 0 vezes.
  - A direção da linha dos centros desvia das eqs. 2.3b/2.4b no máximo 2,3e-41.

### R3, R4–R6: a trajetória simétrica das Figs. 6 e 7

A trajetória é perpendicular à linha dos centros e passa pelo ponto médio, com v₁ = 0,5.

- **R3, Fig. 6b, troca de velocidades** (B = 0,01, R = 1e-4): 1 solução. A rapidez entre colisões menos (1 − 2B) dá -0,000178. Houve 3 intervalos não resolvidos.
- **R3, Fig. 6b, troca de velocidades** (B = 0,01, R = 1e-5): 1 solução. A rapidez entre colisões menos (1 − 2B) dá -1,78e-5. Houve 7 intervalos não resolvidos.
- **R3, Fig. 6b, troca de velocidades** (B = 0,01, R = 1e-6): 1 solução. A rapidez entre colisões menos (1 − 2B) dá -1,78e-6. Houve 46 intervalos não resolvidos.
  - Desvio / R = -1,781, -1,782, -1,782: o desvio é proporcional a R, e v₁′ = 1 − 2B (p. 1083) é o limite R → 0.
- **Limite v₁ → 0 (p. 1084).** Com a bola em repouso sobre a linha dos centros em x = 0,3; 0,5; 0,7, o verificador aceita a solução: sim; sim; sim.
  - A rapidez do laço é 0,979998; menos (1 − 2B − 2R), dá 1,13e-41. O texto dá v₂ = 1, no limite B ≪ 1.
  - Nada na história depende do instante da colisão, que é o contínuo de soluções de que fala o artigo.

| B | soluções com 1, 2 e 3 travessias | tan ψ / v₁ − 1 (eq. 3.1) | (sin ψ + tan ψ)/(2v₁) − 1 (eq. 3.2) | instantes da p. 1084, maior desvio | intervalos não resolvidos (1/2/3) |
|---|---|---|---|---|---|
| 0,04 | y = -0,268948581; y = -0,278565901; y = -0,286884495 | 0,07579 | 0,04383 | 0,00257 | 48/44/42 |
| 0,02 | y = -0,259204821; y = -0,27177034; y = -0,279178381 | 0,03682 | 0,0211 | 0,0012 | 47/57/47 |
| 0,01 | y = -0,254536707; y = -0,268574866; y = -0,275548933 | 0,01815 | 0,01035 | 0,000583 | 53/57/49 |
| 0,005 | y = -0,25225214; y = -0,267024647; y = -0,273786081 | 0,009009 | 0,005129 | 0,000287 | 46/56/50 |

ψ é medido na geometria: tan ψ = |y|/(1/2), onde y é a altura da colisão. Ajustando log|desvio| contra log B, a ordem medida é 1,02 para a eq. 3.1 e 1,03 para a eq. 3.2. As duas são, portanto, o limite B → 0 da geometria. Para a Fig. 7d (três travessias), o artigo não imprime valor ("exercise for the reader"), e a coluna dá o y medido.

### R7: muitas travessias (seção III B)

**v₁ = 3, ψ_A = 0,6, h = 0,5**, B = 0,05, R = 10⁻⁴, colisão a até 50 D. O n_min da eq. 3.8 é 1,125.

| n | L medido | n da eq. 3.7 no L medido | intervalos não resolvidos | fração do espaço de busca, conservadora (fina) |
|---|---|---|---|---|
| 1 | 1,906699 | 2,0104 | 0 | 0,0% (0,0%) |
| 2 | 2,956811 | 3,0216 | 0 | 0,0% (0,0%) |
| 3 | 4,011271 | 4,0586 | 0 | 0,0% (0,0%) |
| 4 | 5,064029 | 5,1013 | 0 | 0,0% (0,0%) |
| 5 | 6,115843 | 6,1467 | 0 | 0,0% (0,0%) |
| 6 | 7,167124 | 7,1934 | 11207 | 8,66% (0,0309%) |

**v₁ = 2, ψ_A = 0,6, h = 1**, B = 0,05, R = 10⁻⁴, colisão a até 50 D. O n_min da eq. 3.8 é 5,196.

| n | L medido | n da eq. 3.7 no L medido | intervalos não resolvidos | fração do espaço de busca, conservadora (fina) |
|---|---|---|---|---|
| 1 | nenhuma solução encontrada | — | 2 | 0,000128% (0,000128%) |
| 2 | nenhuma solução encontrada | — | 18 | 0,00167% (0,00167%) |
| 3 | nenhuma solução encontrada | — | 24 | 0,00214% (0,00214%) |
| 4 | nenhuma solução encontrada | — | 34 | 0,271% (0,271%) |
| 5 | 0,9403691 | -14,374 | 132 | 1,05% (0,213%) |
| 6 | nenhuma solução encontrada | — | 14010 | 2,47% (0,175%) |

- **Primeiro conjunto:** soluções encontradas para n = 1, 2, 3, 4, 5, 6. A eq. 3.7, avaliada no L medido, dá n mais 1,01 a 1,19 (n ≫ 1, L ≫ 1 é o regime dela). O menor n encontrado é 1; a eq. 3.8, também assintótica, dá 1,125.
- **Segundo conjunto:** soluções encontradas só para n = 5, com 14220 intervalos não resolvidos no total, sem refinamento; a fração do espaço de busca que cobrem está na tabela, pelas medidas definidas em R8 (caso vazio). A busca é **inconclusiva** aqui: a ausência das outras soluções não é afirmada.

### R8: árbitro das eqs. 4.2–4.3 contra a geometria independente

Para cada solução geométrica verificada (colisão antes de uma travessia), α, β e d são medidos na geometria e postos nas eqs. 4.2 e 4.3 como impressas. Em separado, as eqs. 4.2–4.3 são resolvidas diretamente, por findroot a partir de uma grade de pontos de partida, e suas raízes são confrontadas com as geométricas.

| (θ, φ, λ₁, B, R) | movimento livre consistente? | soluções geométricas | raízes físicas das eqs. 4.2–4.3 | coincidem | raízes com distância negativa | maior resíduo das eqs. 4.2–4.3 | maior desvio de d (4.4b) e de v₂ (4.5) | sign(β) = sign(d) (4.6c) | maior desvio de a (A6) contra PQ medido |
|---|---|---|---|---|---|---|---|---|---|
| (0,2; 0,6; 0,5; 0,01; 0,001) | não | 3 (I, II, II) | 3 | 3/3 | 0 | 1,22e-41 | 4,25e-40 | sim | 1,48e-37 |
| (0,2; 0,6; 0,25; 0,2; 0,05) | não | 1 (I) | 1 | 1/1 | 0 | 6,91e-42 | 1,15e-41 | sim | 1,11e-40 |
| (0,0; 0,4; 0,75; 0,3; 0,135) | não | 0 | 0 | 0/0 | 0 | — | — | — | — |
| (0,05; 0,1; 0,5; 0,01; 0,001) | não | 1 (I) | 1 | 1/1 | 0 | 1,6e-42 | 3,44e-41 | sim | 4,9e-39 |
| (1,0; 1,05; 0,5; 0,01; 0,001) | não | 2 (I, II) | 2 | 2/2 | 1 | 5,07e-42 | 4,59e-41 | sim | 5,51e-38 |
| (-0,78; 0,9; 0,5; 0,1; 0,02) | não | 1 (I) | 1 | 1/1 | 0 | 7,53e-42 | 7,17e-43 | sim | 6,89e-40 |

As raízes das eqs. 4.2–4.3 são procuradas só no domínio 0 < α < π, |β| < π/2.

- **Árbitro.** Nos 5 casos com solução geométrica, as raízes físicas das eqs. 4.2–4.3 e as soluções geométricas são as mesmas, uma a uma.
  - O maior resíduo ou desvio entre resíduo das eqs. 4.2–4.3, d (4.4b), v₂ (4.5) e a (A6) é 1,48e-37.
  - A regra de sinal sign(β) = sign(d) da eq. 4.6c vale em todas as 8 soluções. A eq. A6 é comparada com a distância PQ medida em cada uma.
  - Isso confere a transcrição das eqs. 4.2, 4.3, 4.4a, 4.4b, 4.5 e A6, e o modelo da travessia, nestes casos.
- **Raízes que a eq. 4.6 não lista,** com alguma das distâncias a, a′ ou c do Apêndice A negativa: (θ, φ) = (1,0; 1,05): 1.
  - A geometria não a realiza, porque nenhuma bola percorre distância negativa ao longo do caminho.
  - A lista da eq. 4.6 cobre as raízes introduzidas pelas manipulações do Apêndice A (p. 1097). Esta tem outra origem.
- **Caso (0,0; 0,4; 0,75; 0,3; 0,135):** a trajetória é perigosa, porque o movimento livre se sobrepõe a si mesmo (distância mínima 0,135 < 2R). Com uma colisão antes da travessia, 1 travessia entre as visitas e as grades `ekt_*` padrão: **nenhuma solução encontrada na região resolvida; 424 intervalos não resolvidos, cobrindo 10,4% do espaço de busca pela medida conservadora (0,000179% pela medida fina).**
  - **Espaço de busca:** Ω = [s_min; s_max] × [0; 2π) em (posição da colisão s, ângulo ν da normal de contato); a mira é resolvida dentro de cada ponto de Ω.
  - **Medida fina:** soma das células (s × largura em ν) dos intervalos não resolvidos. Intervalos em ν registrados fora da grade de s (dentro de uma bisseção em s) não têm célula e não entram nela.
  - **Medida conservadora:** fração de s coberta pela união dos intervalos da grade que cada item não resolvido pode estragar, com ν inteiro. Um intervalo em ν não resolvido num ponto da grade pode esconder um ramo de raízes e estragar a ligação de v₁(s) nos dois intervalos vizinhos; um item fora da grade conta com o intervalo que o contém.
  - A medida conservadora é cota superior da região em que uma solução poderia ter escapado **só sob o critério de bracketing e na grade declarada**. Ela não cobre uma troca de sinal dupla entre dois nós. Os colchetes da busca da mira descartados (f indefinida na bisseção, ou salto) são registrados e contados como itens no ponto (s, ν) em que ocorreram. Não são registrados os intervalos da grade da mira com uma ponta indefinida (a borda do domínio do laço, along ≤ 0): uma raiz entre o último nó definido e essa borda não aparece em contagem nenhuma.
  - **Refinamento adaptativo** (até `ekt_refine_max_split` = 30 halvings a mais em cada intervalo): de 263 intervalos não resolvidos (mira, colchete perdido (f indefinida na bisseção): 3; ν, contagem mudou: 256; s, contagem mudou: 4), cobrindo 30,1% pela medida conservadora, sobram 424 (mira, colchete perdido (f indefinida na bisseção): 394; ν, contagem mudou: 22; s, contagem mudou: 8); 145 deles estão fora da grade de s e não entram na medida fina. Motivos dos que sobram: 391: aim_bracket_lost: an aim bracket dropped during the refinement; 21: a matched branch changes sign: a middle root the main search did not follow; 8: the middle-root search at a refined end of this leaf logged unresolved nu-intervals; 3: aim_bracket_lost; 1: aim-root count changes and C takes both signs at the ends.
  - Critério de resolvido: num trecho em que a contagem de raízes muda, a grandeza que decide a solução (a troca normal C, em ν; v₁ menos o alvo, em s) tem um único sinal em todas as raízes das duas pontas. É o critério de bracketing da busca principal, não uma prova.
  - Intervalos com bisseção perdida não são refinados e continuam contados. 0 soluções verificadas nesta busca; 1 candidatos rejeitados.
  - **Certificação intervalar** (`certify.py`, mpmath.iv a 53 bits, arredondamento para fora; enunciado e hipóteses em `docs/phase3_certification.md`).
    - **Itens não resolvidos depois do refinamento:** 424, sendo 30 do mesmo tipo dos que já eram contados e 394 colchetes da busca da mira descartados (registrados a partir desta versão). Resultado: 424 excluídos (30 dos 30; 394 dos 394), 0 indecididos; nenhuma raiz certificada, porque o teste de Krawczyk não foi implementado nesta versão.
    - Caixa de cada item: o intervalo em s que ele pode estragar (o da medida conservadora, arredondado para fora) × ν e mira inteiros. 27 caixas distintas, 165343 subcaixas examinadas.
    - **Ω inteiro, com uma travessia:** 656220 subcaixas em 8 pedaços de s, todas excluídas (miss ≠ 0: 142943; C ≠ 0: 84835; impulso ≤ 0: 68222; V ≠ 0: 24706; bola velha ainda não saiu da boca 1 (along ≤ 0): 4490; entrada na boca 2 impossível (d·n_X ≥ 0): 2918).
    - **Provado, em aritmética intervalar, que não existe solução de uma colisão para esta condição inicial no modelo de [uma colisão; a bola jovem colide na linha inicial (Fig. 11) antes da boca 2, em s ∈ [-3,0; -0,00135]; exatamente 1 travessia entre as duas visitas (entra na boca 2 e sai da boca 1 uma unidade de tempo antes); entrada transversal na boca 2; a bola velha já fora da boca 1 no instante da colisão; impulso normal > 0; modelo geométrico de billiards.py (EKT 1991, Sec. II A, p. 1081); mpmath.iv a 53 bits, arredondamento para fora]; soluções com várias colisões não foram examinadas.** Origem: CERTIFIED, rótulo “exclusão certificada” (rótulo proposto, à espera de aprovação); soluções na caixa: 0.
    - **Não coberto pela prova:** duas ou mais travessias entre as visitas, colisão depois de uma travessia, e posições de colisão fora desse intervalo de s. Para esses vale o que a busca por bracketing diz acima e nas variações abaixo (COMPUTED).
    - Nenhuma caixa ficou indecidida, então nenhuma raiz precisou do teste de Krawczyk; ele não foi implementado nesta versão.
  - Posição da colisão s (distância com sinal ao longo da linha, a partir da entrada na boca 2; s < 0 antes dela) em [-3; -0,00135] antes da travessia, e em [1,0e-6; 3] depois dela. A grade de s tem `ekt_s_grid` pontos uniformes mais aproximações geométricas às pontas; ν e a mira usam `ekt_nu_grid` e `ekt_theta_grid`.
  - Mesma busca, antes do refinamento: 0 raízes físicas das eqs. 4.2–4.3, 1 candidatos geométricos rejeitados (v1(s) crosses the target across a pole, not a root), 263 intervalos não resolvidos (mira, colchete perdido (f indefinida na bisseção): 3; ν, contagem mudou: 256; s, contagem mudou: 4), cobrindo 30,1% do espaço de busca pela medida conservadora (0,000558% pela medida fina).
  - Variações, sem refinamento:
    - Colisão antes da travessia, 2 travessias entre as visitas: 0 soluções, 2 candidatos rejeitados, 735 intervalos não resolvidos (ν, bisseção perdeu o ramo: 7; ν, contagem mudou: 696; s, contagem mudou: 32), cobrindo 57,5% do espaço de busca pela medida conservadora (0,00305% pela medida fina).
    - Colisão antes da travessia, 3 travessias entre as visitas: 0 soluções, 2 candidatos rejeitados, 682 intervalos não resolvidos (ν, bisseção perdeu o ramo: 6; ν, contagem mudou: 648; s, contagem mudou: 28), cobrindo 58,6% do espaço de busca pela medida conservadora (0,00312% pela medida fina).
    - Colisão depois de uma travessia, 1 travessia entre as visitas: 0 soluções, 1 candidatos rejeitados, 230 intervalos não resolvidos (ν, contagem mudou: 228; s, contagem mudou: 2), cobrindo 9,21% do espaço de busca pela medida conservadora (0,000296% pela medida fina).
    - Colisão depois de uma travessia, 2 travessias entre as visitas: 0 soluções, 0 candidatos rejeitados, 186 intervalos não resolvidos (ν, bisseção perdeu o ramo: 7; ν, contagem mudou: 176; s, contagem mudou: 3), cobrindo 9,21% do espaço de busca pela medida conservadora (0,000398% pela medida fina).
    - Colisão antes da travessia, 1 travessia entre as visitas, grades dobradas (ν, mira e s): 0 soluções, 1 candidatos rejeitados, 347 intervalos não resolvidos (mira, colchete perdido (f indefinida na bisseção): 1; ν, contagem mudou: 344; s, contagem mudou: 2), cobrindo 29,4% do espaço de busca pela medida conservadora (0,000242% pela medida fina).
    - Colisão antes da travessia, 2 travessias entre as visitas, grades dobradas (ν, mira e s): 0 soluções, 2 candidatos rejeitados, 796 intervalos não resolvidos (ν, contagem mudou: 760; s, contagem mudou: 36), cobrindo 66,4% do espaço de busca pela medida conservadora (0,00198% pela medida fina).
    - Colisão antes da travessia, 3 travessias entre as visitas, grades dobradas (ν, mira e s): 0 soluções, 2 candidatos rejeitados, 1646 intervalos não resolvidos (mira, colchete rejeitado como salto: 2; ν, bisseção perdeu o ramo: 6; ν, contagem mudou: 1570; s, contagem mudou: 68), cobrindo 73,1% do espaço de busca pela medida conservadora (0,00353% pela medida fina).
    - Colisão antes da travessia, 1 travessia entre as visitas, λ₁ = 0,1: 0 soluções, 2 candidatos rejeitados, 263 intervalos não resolvidos (mira, colchete perdido (f indefinida na bisseção): 3; ν, contagem mudou: 256; s, contagem mudou: 4), cobrindo 30,1% do espaço de busca pela medida conservadora (0,000558% pela medida fina).
    - Colisão antes da travessia, 1 travessia entre as visitas, λ₁ = 0,25: 0 soluções, 1 candidatos rejeitados, 263 intervalos não resolvidos (mira, colchete perdido (f indefinida na bisseção): 3; ν, contagem mudou: 256; s, contagem mudou: 4), cobrindo 30,1% do espaço de busca pela medida conservadora (0,000558% pela medida fina).
    - Colisão antes da travessia, 1 travessia entre as visitas, λ₁ = 0,5: 0 soluções, 1 candidatos rejeitados, 263 intervalos não resolvidos (mira, colchete perdido (f indefinida na bisseção): 3; ν, contagem mudou: 256; s, contagem mudou: 4), cobrindo 30,1% do espaço de busca pela medida conservadora (0,000558% pela medida fina).
    - Colisão antes da travessia, 1 travessia entre as visitas, λ₁ = 0,9: 0 soluções, 1 candidatos rejeitados, 263 intervalos não resolvidos (mira, colchete perdido (f indefinida na bisseção): 3; ν, contagem mudou: 256; s, contagem mudou: 4), cobrindo 30,1% do espaço de busca pela medida conservadora (0,000558% pela medida fina).
  - Isto **não** é declarado paradoxo nem multiplicidade zero. Soluções com várias colisões, que o artigo menciona (p. 1080), estão fora do modelo buscado; a grade é finita; e os intervalos não resolvidos acima são lugares onde uma solução poderia ter escapado. É o primeiro alvo do mapa (pergunta Q1).

### R9: região central do triângulo, eqs. 4.12–4.13 e as duas formas da eq. B20

Parâmetros: (θ, φ) = (0,2; 0,6), λ₁ = ½, R = 0,1 B.

| B | classe | β / (4.12a ou 4.13a) − 1 | ε / (4.12b ou 4.13b) − 1 | (λ₂ − σ) / (4.12c ou 4.13c) − 1 | (λ₂ − σ) / (B20, 1ª forma) − 1 |
|---|---|---|---|---|---|
| 0,04 | I | 0,2077 | 0,22 | 0,01366 | 0,1269 |
| 0,04 | II | 0,3144 | 0,2483 | 0,03446 | 0,15 |
| 0,02 | I | 0,09972 | 0,1038 | 0,003738 | 0,06084 |
| 0,02 | II | 0,1351 | 0,109 | 0,007865 | 0,06521 |
| 0,01 | I | 0,04875 | 0,05033 | 0,0009735 | 0,02972 |
| 0,01 | II | 0,06327 | 0,05145 | 0,001903 | 0,03068 |
| 0,005 | I | 0,02409 | 0,02478 | 0,0002482 | 0,01468 |
| 0,005 | II | 0,03067 | 0,02504 | 0,0004693 | 0,0149 |

- **β e ε (classe I):** os desvios das eqs. 4.12a–b caem com ordem 1,04 e 1,05 em B. As formas analíticas são o limite B → 0 da geometria.
- **λ₂ − σ:** com a segunda forma da eq. B20 (igual à 4.12c), o desvio cai com ordem 1,93, então o termo de ordem B confere. Com a primeira forma impressa, cai com ordem 1,04, então o coeficiente de ordem B dela não confere com a geometria.
- **As duas formas impressas da eq. B20 diferem por um fator 4 no termo em B** (ver `docs/equation_transcription.md`). A geometria segue a segunda forma.

### R10: borda superior θ = φ (Fig. 14), B = 10⁻², R = 10⁻⁵, λ₁ = ½

| φ | β classe I (geometria) | eq. B32 | eq. 4.14 | eq. 4.15 | classe II (geometria) |
|---|---|---|---|---|---|
| 0,05 | 0,001483 | 0,00049825 | 0,0004973 | -1,005e-6 | nenhuma |
| 0,1 | 0,0014612 | 0,00098029 | 0,0009784 | -2,037e-6 | nenhuma |
| 0,2 | 0,002047 | 0,0018335 | 0,00183 | -4,314e-6 | nenhuma |
| 0,3 | 0,0025618 | 0,0024448 | 0,002439 | -7,161e-6 | nenhuma |
| 0,5 | 0,0026422 | 0,0026049 | 0,00259 | -1,775e-5 | nenhuma |
| 0,7 | 0,0011649 | 0,0011614 | 0,001095 | -7,58e-5 | nenhuma |
| 0,75 | 0,00062446 | 0,00062527 | 0,0004822 | -0,0001927 | nenhuma |
| 0,78 | 0,0003506 | 0,00035152 | 7,593e-5 | -0,001303 | nenhuma |
| 0,8 | 0,00023052 | 0,00023069 | -0,0002095 | 0,0004913 | nenhuma |
| 0,82 | 0,00015967 | 0,00015932 | -0,0005056 | 0,0002115 | nenhuma |
| 0,9 | 6,6744e-5 | 6,6233e-5 | -0,00178 | 6,895e-5 | nenhuma |
| 1,1 | 3,0287e-5 | 3,0067e-5 | -0,005245 | 3,029e-5 | nenhuma |
| 1,3 | 2,2497e-5 | 2,2403e-5 | -0,008257 | 2,249e-5 | nenhuma |
| 1,5 | 2,0099e-5 | 2,009e-5 | -0,009875 | 2,015e-5 | nenhuma |

- **Classe I.** Para φ entre 0,7 e 0,9, a geometria acompanha a eq. B32 com desvio relativo máximo de 0,0077. Longe de π/4 ela se aproxima das eqs. 4.14 e 4.15, cada uma do seu lado. Tudo isso num único B = 10⁻²: nenhum limite B → 0 foi medido para B32, 4.14 ou 4.15.
  - Perto de φ → 0, β volta a crescer, e a eq. B32 deixa de acompanhar. O artigo atribui isso à quebra da expansão do Apêndice B perto do canto (p. 1092), e é onde a Fig. 14 mostra a parte numérica.
  - A forma da curva confere, a olho, com a Fig. 14: mínimo perto de φ = 0,1, máximo de cerca de 2,6 × 10⁻³ perto de φ = 0,5, e cerca de 2 × 10⁻⁵ em φ → π/2. A comparação com a figura é visual; não há números impressos.
- **Classe II.** Com a colisão antes da travessia, a geometria não encontra solução de classe II em nenhum φ da borda.
  - φ = 0,3: as eqs. 4.2–4.3 têm 0 raízes físicas de classe II e 0 de classe II com distância negativa (β = —, a = —). A eq. B29 com s = −1 dá β = 0,002433; 7,178e-6.
  - φ = 0,9: as eqs. 4.2–4.3 têm 0 raízes físicas de classe II e 2 de classe II com distância negativa (β = -0,001698; -7,349e-5, a = -0,0011; -5,71e-5). A eq. B29 com s = −1 dá β = -7,216e-5; -0,001701.
  - φ = 1,1: as eqs. 4.2–4.3 têm 0 raízes físicas de classe II e 2 de classe II com distância negativa (β = -0,005181; -3,085e-5, a = -0,00165; -1,55e-5). A eq. B29 com s = −1 dá β = -3,051e-5; -0,005206.
  - Onde a eq. B29 prevê classe II (φ > π/4), as raízes correspondentes das eqs. 4.2–4.3 têm a = PQ < 0: a colisão ficaria além da boca 2 ao longo do caminho, e a geometria não as realiza com a colisão antes da travessia.
  - O artigo suspeita (p. 1091, lacuna 16) que a classe II que falta exista com a colisão depois da travessia. R11 procura essa solução.

### R11: classe II com a colisão depois de uma travessia (lacuna 16, p. 1091)

| φ (borda superior) | soluções | β | posição após a saída da boca 1 | intervalos não resolvidos |
|---|---|---|---|---|
| 0,7 | 0 | — | — | 45 |
| 0,75 | 0 | — | — | 45 |
| 0,78 | 0 | — | — | 57 |
| 0,79 | 0 | — | — | 47 |
| 0,82 | 0 | — | — | 41 |
| 0,85 | 1 | -0,0007575 | 0,0008133 | 43 |

- O primeiro encontro com a colisão fica sobre a linha que sai livremente da boca 1, e há mais uma travessia entre as visitas.
- Soluções verificadas aparecem em φ = 0,85. Em φ = 0,7, 0,75, 0,78, 0,79, 0,82, nenhuma solução foi encontrada sob este critério.
- Só foi buscado um tipo de "colisão depois da travessia": uma travessia antes e uma entre as visitas. O artigo não especifica o que tem em mente.

### O que a reprodução estabelece e o que não estabelece

- ✅ **Estabelece:**
  - As eqs. 3.1, 3.2 (Figs. 6 e 7) e 4.12–4.13 são o limite B → 0 de um modelo geométrico exato e independente: os desvios caem com ordem ≈ 1 em B (R4–R6, R9).
  - As eqs. 4.14, 4.15 e B32 (Fig. 14) foram comparadas só em B = 10⁻²: a classe I acompanha a B32 perto de φ = π/4 (R10). Isso não é um teste de limite.
  - Na seção III B (Fig. 9), há uma solução para cada n de 1 a 6 no primeiro conjunto de R7, com a eq. 3.7 conferindo só assintoticamente. O segundo conjunto é inconclusivo.
  - As eqs. 4.2–4.3 são satisfeitas pela geometria no arredondamento de 40 dígitos, nos casos de R8 com solução (ver a tabela de R8).
- ❌ **Não estabelece:**
  - Nenhuma afirmação universal. As contagens valem para as grades, os itinerários e as janelas declarados.
  - Nada sobre soluções com várias colisões.
  - **A fase é clássica.** Ela testa a alternativa de autoconsistência de Novikov, não o mecanismo quântico de Hawking.
  - O mapa sistemático e a certificação intervalar ainda não foram feitos.
  - O teste de escala que confirmaria a observação deste projeto, de que só v₁T_d/D importa (a observação não é do artigo), ainda não foi feito.

