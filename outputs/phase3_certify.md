# chronosim · fase 3, certificação intervalar (`python -m chronosim.run --certify`)

Grandezas CERTIFIED desta fase; tudo o mais é COMPUTED. Enunciados em `docs/phase3_certification.md` e `docs/phase3_borda.md`.

## Caso vazio do R8 (0,0; 0,4; 0,75; 0,3; 0,135), intervalo principal (máquina pela mira)

  - **Certificação intervalar** (`certify.py`, mpmath.iv a 53 bits, arredondamento para fora; enunciado e hipóteses em `docs/phase3_certification.md`).
    - **Itens não resolvidos depois do refinamento:** 424, sendo 30 do mesmo tipo dos que já eram contados e 394 colchetes da busca da mira descartados (registrados a partir desta versão). Resultado: 424 excluídos (30 dos 30; 394 dos 394), 0 indecididos; o teste de Krawczyk não foi aplicado a esses itens.
    - Caixa de cada item: o intervalo em s que ele pode estragar (o da medida conservadora, arredondado para fora) × ν e mira inteiros. 27 caixas distintas, 165343 subcaixas examinadas.
    - **Ω inteiro, com uma travessia:** 656220 subcaixas em 8 pedaços de s, todas excluídas (miss ≠ 0: 142943; C ≠ 0: 84835; impulso ≤ 0: 68222; V ≠ 0: 24706; bola velha ainda não saiu da boca 1 (along ≤ 0): 4490; entrada na boca 2 impossível (d·n_X ≥ 0): 2918).
    - **Provado, em aritmética intervalar, que não existe solução de uma colisão para esta condição inicial no modelo de [uma colisão; a bola jovem colide na linha inicial (Fig. 11) antes da boca 2, em s ∈ [-3,0; -0,00135]; exatamente 1 travessia entre as duas visitas (entra na boca 2 e sai da boca 1 uma unidade de tempo antes); entrada transversal na boca 2; a bola velha já fora da boca 1 no instante da colisão; impulso normal > 0; modelo geométrico de billiards.py (EKT 1991, Sec. II A, p. 1081); mpmath.iv a 53 bits, arredondamento para fora]; soluções com várias colisões não foram examinadas.** Origem: CERTIFIED, rótulo “exclusão certificada”; soluções na caixa: 0.
    - **Controle positivo** (a mesma máquina, com o mesmo código, a mesma precisão e as mesmas hipóteses, só com os parâmetros do caso; o espaço de busca de cada caso é dividido em pedaços de s; critério: nenhuma caixa que contenha uma solução conhecida, ou uma cópia periódica dela, pode ser excluída):

      | caso | pedaços de s examinados até o fim | com solução: examinados até o fim | volume examinado | caixas examinadas | excluídas | indecididas (largura mínima) | não examinadas (orçamento) | soluções conhecidas | em caixa excluída | Krawczyk: raiz única | Krawczyk nas caixas sem solução conhecida: sem raiz / indecidido / raiz sem solução conhecida |
      |---|---|---|---|---|---|---|---|---|---|---|---|
      | R8 (0,2; 0,6; 0,5; 0,01; 0,001) | 13 de 19 | 2 de 2 | 99,997 % | 3423709 | 1704329 | 7496 | 78 | 3 | 0 | 3 de 3 | 7493 / 0 / 0 |
      | R8 (0,2; 0,6; 0,25; 0,2; 0,05) | 11 de 13 | 1 de 1 | 99,9995 % | 1501059 | 748985 | 1538 | 26 | 1 | 0 | 1 de 1 | 1537 / 0 / 0 |
      | R8 (0,05; 0,1; 0,5; 0,01; 0,001) | 9 de 19 | 1 de 1 | 98,3 % | 14056719 | 6953329 | 74958 | 164 | 1 | 0 | 1 de 1 | 74956 / 0 / 0 |
      | R8 (1,0; 1,05; 0,5; 0,01; 0,001) | 10 de 19 | 1 de 1 | 99,97 % | 3414084 | 1704065 | 2913 | 147 | 2 | 0 | 2 de 2 | 2911 / 0 / 0 |
      | R8 (-0,78; 0,9; 0,5; 0,1; 0,02) | 14 de 14 | 1 de 1 | 100 % | 544920 | 272317 | 150 | 0 | 1 | 0 | 1 de 1 | 149 / 0 / 0 |
      | Fig. 7(b), B = 0,04, R = 1e-7, v₁ = 0,5 | 28 de 28 | 2 de 2 | 100 % | 134600 | 66781 | 533 | 0 | 2 | 0 | 2 de 2 | 529 / 0 / 0 |
      | Fig. 7(b), B = 0,02, R = 1e-7, v₁ = 0,5 | 28 de 28 | 2 de 2 | 100 % | 125104 | 62084 | 482 | 0 | 2 | 0 | 2 de 2 | 478 / 0 / 0 |
      | Fig. 7(b), B = 0,01, R = 1e-7, v₁ = 0,5 | 28 de 28 | 2 de 2 | 100 % | 114802 | 56962 | 453 | 0 | 2 | 0 | 2 de 2 | 449 / 0 / 0 |
      | Fig. 7(b), B = 0,005, R = 1e-7, v₁ = 0,5 | 28 de 28 | 2 de 2 | 100 % | 114604 | 56866 | 450 | 0 | 2 | 0 | 2 de 2 | 446 / 0 / 0 |

      Onde estão as soluções conhecidas (cada solução contada uma vez):

      | caso | solução (s; ν; mira) | pedaço de s | cópias periódicas na caixa | em caixa excluída | caixa indecidida (maior lado em ν, mira) | Krawczyk |
      |---|---|---|---|---|---|---|
      | R8 (0,2; 0,6; 0,5; 0,01; 0,001) | -0,41364649; 0,19755483; -2,5392041 | [-0,75; -0,375] | 1 | 0 | 5,9e-9 | raiz única em cada cópia, raiz certificada do modelo geométrico (billiards.py) (CERTIFIED) |
      | R8 (0,2; 0,6; 0,5; 0,01; 0,001) | -0,40816375; 3,3440408; -2,5439834 | [-0,75; -0,375] | 1 | 0 | 5,9e-9 | raiz única em cada cópia, raiz certificada do modelo geométrico (billiards.py) (CERTIFIED) |
      | R8 (0,2; 0,6; 0,5; 0,01; 0,001) | -0,020857865; 3,5689055; -2,6983111 | [-0,02344; -0,01172] | 1 | 0 | 5,9e-9 | raiz única em cada cópia, raiz certificada do modelo geométrico (billiards.py) (CERTIFIED) |
      | R8 (0,2; 0,6; 0,25; 0,2; 0,05) | -0,37394705; 0,036956525; -2,4352158 | [-0,375; -0,1875] | 1 | 0 | 5,9e-9 | raiz única em cada cópia, raiz certificada do modelo geométrico (billiards.py) (CERTIFIED) |
      | R8 (0,05; 0,1; 0,5; 0,01; 0,001) | -0,52569286; 6,2733988; -2,9829218 | [-0,75; -0,375] | 2 | 0 | 5,9e-9 | raiz única em cada cópia, raiz certificada do modelo geométrico (billiards.py) (CERTIFIED) |
      | R8 (1,0; 1,05; 0,5; 0,01; 0,001) | -0,058025721; 0,99880141; -2,0904963 | [-0,09375; -0,04688] | 1 | 0 | 5,9e-9 | raiz única em cada cópia, raiz certificada do modelo geométrico (billiards.py) (CERTIFIED) |
      | R8 (1,0; 1,05; 0,5; 0,01; 0,001) | -0,056618444; 4,142793; -2,0926885 | [-0,09375; -0,04688] | 1 | 0 | 5,9e-9 | raiz única em cada cópia, raiz certificada do modelo geométrico (billiards.py) (CERTIFIED) |
      | R8 (-0,78; 0,9; 0,5; 0,1; 0,02) | -0,91676232; 5,2823579; -2,034807 | [-1,5; -0,75] | 1 | 0 | 5,9e-9 | raiz única em cada cópia, raiz certificada do modelo geométrico (billiards.py) (CERTIFIED) |
      | Fig. 7(b), B = 0,04, R = 1e-7, v₁ = 0,5 | -0,26894858; 6,2831852; -2,6480888 | [-0,375; -0,1875] | 2 | 0 | 5,9e-9 | raiz única em cada cópia, raiz certificada do modelo geométrico (billiards.py) (CERTIFIED) |
      | Fig. 7(b), B = 0,04, R = 1e-7, v₁ = 0,5 | -4,7751341e-8; 5,785363; 3,1415927 | [-8,941e-8; -4,47e-8] | 2 | 0 | 5,9e-9 | raiz única em cada cópia, raiz certificada do modelo geométrico (billiards.py) (CERTIFIED) |
      | Fig. 7(b), B = 0,02, R = 1e-7, v₁ = 0,5 | -0,25920482; 6,2831852; -2,663326 | [-0,375; -0,1875] | 2 | 0 | 5,9e-9 | raiz única em cada cópia, raiz certificada do modelo geométrico (billiards.py) (CERTIFIED) |
      | Fig. 7(b), B = 0,02, R = 1e-7, v₁ = 0,5 | -4,619345e-8; 5,8030103; 3,1415927 | [-8,941e-8; -4,47e-8] | 2 | 0 | 5,9e-9 | raiz única em cada cópia, raiz certificada do modelo geométrico (billiards.py) (CERTIFIED) |
      | Fig. 7(b), B = 0,01, R = 1e-7, v₁ = 0,5 | -0,25453671; 6,2831852; -2,6707126 | [-0,375; -0,1875] | 2 | 0 | 5,9e-9 | raiz única em cada cópia, raiz certificada do modelo geométrico (billiards.py) (CERTIFIED) |
      | Fig. 7(b), B = 0,01, R = 1e-7, v₁ = 0,5 | -4,5447042e-8; 5,8114078; 3,1415927 | [-8,941e-8; -4,47e-8] | 2 | 0 | 5,9e-9 | raiz única em cada cópia, raiz certificada do modelo geométrico (billiards.py) (CERTIFIED) |
      | Fig. 7(b), B = 0,005, R = 1e-7, v₁ = 0,5 | -0,25225214; 6,2831852; -2,674348 | [-0,375; -0,1875] | 2 | 0 | 5,9e-9 | raiz única em cada cópia, raiz certificada do modelo geométrico (billiards.py) (CERTIFIED) |
      | Fig. 7(b), B = 0,005, R = 1e-7, v₁ = 0,5 | -4,5081657e-8; 5,8155055; 3,1415927 | [-8,941e-8; -4,47e-8] | 2 | 0 | 5,9e-9 | raiz única em cada cópia, raiz certificada do modelo geométrico (billiards.py) (CERTIFIED) |

      - **O controle passou:** nenhuma das 16 soluções conhecidas (25 pontos, contando as cópias periódicas) caiu em caixa excluída; todas ficaram em caixas indecididas de largura mínima, e o Krawczyk provou raiz única em 16 delas, em todas as cópias.
      - **Só os pedaços com solução foram examinados até o fim** em: R8 (0,2; 0,6; 0,5; 0,01; 0,001); R8 (0,2; 0,6; 0,25; 0,2; 0,05); R8 (0,05; 0,1; 0,5; 0,01; 0,001); R8 (1,0; 1,05; 0,5; 0,01; 0,001). Nos outros pedaços desses casos, o orçamento acabou antes, e o que sobrou está na coluna “não examinadas”; ali nada é afirmado, nem exclusão nem raiz.
      - **Cópias periódicas:** a caixa cobre ν ∈ [-0,01; 6,29319] e mira ∈ [-3,15; 3,15], além de [0; 2π) e [−π; π). Uma solução com ν perto de 0 ou de 2π, ou com mira perto de ±π, aparece duas vezes na caixa (quatro, se as duas coisas valem); cada cópia é vigiada e testada, mas a solução é contada uma vez.
      - **Rótulo:** cada solução certificada recebe a origem CERTIFIED com o rótulo “raiz certificada do modelo geométrico (billiards.py)”. É raiz de (miss, C, V) com as três desigualdades válidas na caixa, e não “solução física certificada”: a ausência de colisões extras não foi provada em intervalos, e vem do verificador (COMPUTED).
      - **Origem de s:** no R8, s = 0 é o ponto em que a linha inicial entra na boca 2. Na Fig. 7(b), s é a altura acima da linha dos centros: s = 0 fica a meia distância das bocas. A solução da Fig. 7(b) com |s| < R é a de troca de velocidades da Fig. 6(b) do artigo (p. 1083), longe das bocas (`docs/phase3_borda.md`).
      - **Caixas examinadas** inclui as que foram divididas ao meio; cada uma delas termina excluída, dividida, indecidida na largura mínima ou não examinada.
      - **Krawczyk:** foi aplicado só às caixas indecididas que a máquina examinou. As caixas que ficaram na pilha por orçamento não foram testadas.
      - **Orçamento:** pedaços sem solução conhecida param em `ekt_cert_control_budget` caixas, e pedaços com solução em `ekt_cert_control_sol_budget`.
    - **Não coberto pela prova:** duas ou mais travessias entre as visitas, colisão depois de uma travessia, e posições de colisão fora desse intervalo de s. Para esses vale o que a busca por bracketing diz acima e nas variações abaixo (COMPUTED).
    - Nenhuma caixa ficou indecidida, então nenhuma raiz precisou do teste de Krawczyk neste caso; o Krawczyk aparece no controle positivo acima.

### Faixa junto à boca do caso vazio (máquina por direção)

- **Caixa:** s ∈ [-0,00135; 0], com s = 0 incluído como superconjunto; ν inteiro; direção de partida α ∈ [φ − π/2 − pad; φ + π/2 + pad] (`ekt_cert_alpha_pad`). Mesma precisão (53 bits) e mesma largura mínima relativa (`ekt_cert_min_width`), medida contra a faixa inteira; 8 pedaços de s, só para paralelizar. Enunciado em `docs/phase3_borda.md`.
- **Resultado:** 3132 subcaixas, todas excluídas (miss ≠ 0: 706; C ≠ 0: 272; impulso ≤ 0: 252; V ≠ 0: 224; bola velha ainda não saiu da boca 1 (along ≤ 0): 116).
- **Verificação cruzada** em s ∈ [-0,01; -0,00135], pedaço que a máquina atual já prova: máquina pela mira: 456195 subcaixas, todas excluídas; máquina por direção: 321 subcaixas, todas excluídas. **As duas concordam:** nenhuma solução.

- **Controle positivo dentro da faixa** (decisão (a)): raízes construídas a partir da história (colisão em s fixo na faixa; v₁ obtido da história e escrito como λ₁ decimal), com B e R do caso vazio e outra linha, (θ; φ) = (0,5; 0,9), com λ₁ fora de [0; 1] (v₁ acima do intervalo da eq. 4.10); mesma máquina e mesma construção de caixa da faixa (s e ν iguais, α em torno do φ do caso). Critério: nenhuma caixa que contenha a raiz pode ser excluída, e o Krawczyk a certifica.

  | controle (θ; φ; λ₁; B; R) | s | verificador independente | d·n_X; along; impulso | subcaixas | em caixa excluída | Krawczyk na caixa da raiz | outras caixas: sem raiz / indecidido / raiz |
  |---|---|---|---|---|---|---|---|
  | 0,5; 0,9; 3,5848820852079515833; 0,3; 0,135 | -0,001 | rejeita: contact point inside a mouth (excluded by the model, p. 1081 item i); return leg enters a mouth before reaching the collision | -0,3314; 0,7581; 3,039 | 785491 | 0 | raiz única, raiz certificada do modelo geométrico (billiards.py) (CERTIFIED) | 4090 / 0 / 0 |
  | 0,5; 0,9; 3,8788614447816208674; 0,3; 0,135 | -1,0e-5 | rejeita: contact point inside a mouth (excluded by the model, p. 1081 item i); return leg enters a mouth before reaching the collision | -0,3211; 0,7581; 3,303 | 617151 | 0 | raiz única, raiz certificada do modelo geométrico (billiards.py) (CERTIFIED) | 5850 / 0 / 0 |
  | 0,5; 0,9; 3,8821746427593530025; 0,3; 0,135 | -1,0e-7 | rejeita: contact point inside a mouth (excluded by the model, p. 1081 item i); return leg enters a mouth before reaching the collision | -0,321; 0,7581; 3,306 | 519839 | 0 | raiz única, raiz certificada do modelo geométrico (billiards.py) (CERTIFIED) | 3572 / 0 / 0 |

  - **O controle passou:** 3 de 3 raízes fora de caixas excluídas; 3 de 3 certificadas pelo Krawczyk (CERTIFIED, CERTIFIED, CERTIFIED).
  - **Não são histórias físicas:** o verificador independente rejeita as 3 raízes polidas (contact point inside a mouth (excluded by the model, p. 1081 item i); return leg enters a mouth before reaching the collision). Os motivos são condições do modelo que a máquina não usa (o item (i) da p. 1081 entre elas). São raízes das equações da máquina com as suas três desigualdades, e testam a máquina, não a física. Nenhuma história física com a colisão na faixa foi construída; o porquê está em `docs/phase3_borda.md`, seção 3.
- **Provado, em aritmética intervalar, que não existe solução de uma colisão para esta condição inicial com a colisão em s ∈ [-0,00135; 0), no modelo de [uma colisão; a bola jovem colide na linha inicial (Fig. 11) antes da boca 2, em s ∈ [-0,00135; 0]; exatamente 1 travessia entre as duas visitas (entra na boca 2 e sai da boca 1 uma unidade de tempo antes); entrada transversal na boca 2; a bola velha já fora da boca 1 no instante da colisão; impulso normal > 0; modelo geométrico de billiards.py (EKT 1991, Sec. II A, p. 1081); mpmath.iv a 53 bits, arredondamento para fora; máquina por direção de partida (certify.conditions_dir), com s = 0 incluído na caixa como superconjunto (em s = 0 a bola jovem já está na boca, fora do modelo); entrada na boca 2 pela primeira interseção do raio].** Origem: CERTIFIED, rótulo “exclusão certificada”. Com a prova do intervalo principal acima, cobre s ∈ [−3; 0) para uma travessia e colisão antes da boca 2.
