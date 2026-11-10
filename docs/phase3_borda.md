# Fase 3 · borda: solução de troca de velocidades da Fig. 7(b) e faixa junto à boca do caso vazio

Registro de 2026-11-09.

- **Seção 1:** verificação da solução da Fig. 7(b), feita.
- **Seção 2:** o plano, aprovado pelo usuário com as decisões (a)–(d) registradas em 2.5.
- **Seção 3:** a implementação e os resultados.

A certificação roda em `python -m chronosim.run --certify` e escreve `outputs/phase3_certify.md`.

## 0. Uma correção de premissa: o s da Fig. 7(b) não é o s do caso vazio

O nome é o mesmo, mas as origens diferem:

- **R8 e caso vazio** (`phase3.line_theta_phi`): p₀ é o ponto em que a linha inicial entra na boca 2, sobre a superfície dela. s = 0 é a boca, e s < 0 é antes dela.
- **Fig. 7(b)** (`phase3.SYM_P0` = (0,5; 0), direção (0; 1)): s é a altura acima da linha dos centros. s = 0 fica a meia distância das bocas, a 0,5 − B de cada uma.

A solução da Fig. 7(b) com s ≈ −4,6·10⁻⁸ fica, portanto, longe das bocas (seção 1). Ela não indica que existam soluções do caso vazio em s ∈ [−0,00135; 0). A faixa continua merecendo exame por outro motivo: o corte em s = −R/100 vem de uma heurística (`ekt_s_near_over_R`), e hoje nada a cobre, nem a busca por bracketing nem a exclusão intervalar.

O relatório do controle positivo passa a dizer qual é a origem de s na Fig. 7(b).

## 1. A solução s ≈ −0,46R da Fig. 7(b)

**Como foi verificada:**
- As soluções foram recalculadas com `phase3.symmetric(B, "1e-7", "0.5", 1)`, a 40 dígitos, para B = 0,04; 0,02; 0,01 e 0,005.
- Cada uma passou por `billiards.verify`, que reconstrói a história inteira com `billiards.history` e confere toda lei do modelo, sem usar as equações do resolvedor.
- Além disso, foram conferidas a 60 dígitos as grandezas abaixo.

Todos os valores são COMPUTED.

| B | s/R | `verify` | folga mín. \|centro − Cᵢ\| − B − R (as duas bolas, as duas bocas) | entrada na boca 2: d·n_X; ângulo com a radial | sai da boca 1 em t | v₁′ − (1 − 2B) |
|---|---|---|---|---|---|---|
| 0,04 | −0,4775 | ok, sem falhas | 0,45999972 | −1,0; 1,04·10⁻⁷ rad | −0,4999999045 | −1,757·10⁻⁷ |
| 0,02 | −0,4619 | ok, sem falhas | 0,47999972 | −1,0; 9,62·10⁻⁸ rad | −0,4999999076 | −1,774·10⁻⁷ |
| 0,01 | −0,4545 | ok, sem falhas | 0,48999972 | −1,0; 9,28·10⁻⁸ rad | −0,4999999091 | −1,782·10⁻⁷ |
| 0,005 | −0,4508 | ok, sem falhas | 0,49499972 | −1,0; 9,11·10⁻⁸ rad | −0,4999999098 | −1,785·10⁻⁷ |

**O que a história mostra** (é a mesma nos quatro casos):
- **Colisão:** a bola jovem está em P = (0,5; s), e a velha em Q ≈ (0,4999998; −s). O ponto de contato fica na linha dos centros. As duas bolas estão inteiramente fora das duas bocas no instante da colisão, com folga ≥ 0,46.
- **Entrada:** a bola jovem segue quase radialmente até a boca 2, a 10⁻⁷ rad da radial, e entra transversalmente (d·n_X = −1) em X ≈ (1 − B; 0).
- **Saída:** ela sai da boca 1 em E ≈ (B; 0), meia unidade de tempo antes da colisão, e volta como bola velha.
- **Depois:** a bola velha segue para cima, na linha inicial.
- **Validade:** |P − Q| − 2R ≈ 10⁻⁶², e não há sobreposição em nenhum par de segmentos. A história depois da colisão não fica truncada.

**No artigo** (EKT 1991, p. 1083, lido no PDF): é a solução de troca de velocidades da Fig. 6(b), e a Fig. 6(c) mostra o caso de bola finita.
- O texto diz, sobre esta mesma trajetória inicial da Fig. 7: "Yet another solution is that of Fig. 6(b), which involves velocity exchange by contrast with the mirror exchange of Fig. 7".
- E a descreve assim: "The ball is hit as it crosses the mouths' line of centers and gets knocked radially into mouth 2 … it is hit with just the right impulse to give it a speed v₁′ = (D − 2B)/T_d ≡ 1 − 2B".
- A diferença v₁′ − (1 − 2B) ≈ −1,8R medida acima é a modificação por raio finito, que o artigo descreve qualitativamente na Fig. 6(c).
- O R3 do relatório já reproduz essa solução, para R = 10⁻⁴ a 10⁻⁶.

**Conclusão:**
- No modelo de `billiards.py`, é uma história consistente: o verificador independente não acusa nenhuma falha. A origem é COMPUTED, e não "solução física certificada", porque a ausência de colisões extras não foi provada em intervalos.
- Não é nova: está no artigo.

**Um cuidado registrado na verificação:** um `mpf` que volta de um processo filho, via pickle, é arredondado para a precisão do processo que o recebe.
- Uma primeira passada, feita com soluções guardadas por um script auxiliar que rodava a 15 dígitos, falhou no verificador por esse motivo, e não por defeito das soluções.
- `run.py` fixa `mp.dps = 80` no processo pai antes de iniciar a fase 3, então as execuções oficiais não são afetadas.

## 2. Faixa junto à boca do caso vazio: s ∈ [−0,00135; 0) (plano, aprovado)

### 2.1 O que o método atual faz ali

A sondagem rodou em 2026-11-09 com o código de `certify.exclude` sem mudanças, no caso vazio (θ = 0; φ = 0,4; λ₁ = 0,75; B = 0,3; R = 0,135).
- Os pedaços são diádicos: [−0,00135·2⁻ʲ; −0,00135·2⁻ʲ⁻¹].
- Cada um tem largura mínima relativa ao próprio pedaço e 1M de caixas de orçamento.
- Rodaram 13 pedaços em paralelo, cerca de 18 min cada.

| j | \|s\| | caixas | indecididas na largura mínima | não examinadas | exclusões por impulso ≤ 0 |
|---|---|---|---|---|---|
| 0–7, 12 | 1,35·10⁻³ a 1,6·10⁻⁷ | 1M (orçamento) | 0 | 19–47 | 378k–500k |
| 20, 30, 40, 50 | 1,3·10⁻⁹ a 6·10⁻¹⁹ | 1M (orçamento) | 264.646 (a mesma contagem nos quatro) | 63 | 235.085 |

**Nenhum pedaço termina.** Há três mecanismos. O primeiro é provável e dominante, e os outros dois são certos:

1. **Divisão por L₀ = |X − P|.** `conditions` calcula V e o impulso com d/L₀.
   - No caso vazio a entrada é radial (θ = 0), então L₀ ≥ |s|, e a meia-abertura do arco visível de P é ≈ √(2|s|/B).
   - Uma caixa de largura Δ em mira tem d com largura ~BΔ e L₀ ≥ |s|. O invólucro de d/L₀ fica então com largura ~BΔ/|s|.
   - Com a largura mínima atual (Δ ≈ 6·10⁻⁹), isso já passa de 1 quando |s| ≲ 2·10⁻⁹. É o que se vê nos quatro pedaços de |s| pequeno, todos com a mesma contagem de caixas indecididas: abaixo dessa escala, s não importa mais.
2. **Superfície do movimento livre.** Na mira de p₀ (d paralelo a e), V = (W − v₁)(e × n) e impulso = (W − v₁)(e·n) se anulam juntos onde W = v₁.
   - Perto dessa superfície, só miss e C podem excluir uma caixa. Isso explica o volume de exclusões por impulso e o esgotamento do orçamento nos pedaços maiores.
   - Distinguir 1 de 2 é o primeiro passo do plano.
3. **Precisão de 53 bits:** P = p₀ + s·e deixa de distinguir P de p₀ quando |s| ≲ 10⁻¹⁶. A aritmética continua rigorosa, porque o arredondamento é para fora, mas perde toda a resolução em s.

### 2.2 Proposta: reparametrizar pela direção de partida e incluir s = 0

**Incógnitas.** Passam a ser (s, ν, α), com α a direção da bola jovem depois da colisão, u = W(cos α, sin α). A mira deixa de ser incógnita: X é a primeira interseção do raio P + L·(cos α, sin α) com a boca 2.
- L₀ = −b − √(b² − c), com b = d̂·(P − C₂) e c = |P − C₂|² − B².
- Como |p₀ − C₂| = B exatamente, c = s·(2 e·(p₀ − C₂) + s) é calculado nessa forma, sem cancelamento. Isso resolve o mecanismo 3 sem subir a precisão.

**Função.** As mesmas condições: miss, C, V, impulso e along. d/L₀ vira d̂, que é a própria incógnita, e não há mais divisão por L₀ (mecanismo 1). W = L₀ + along continua contínua quando L₀ → 0. Entra uma desigualdade nova: o raio atinge a boca (b < 0 e b² − c ≥ 0); a entrada transversal continua sendo d̂·n_X < 0.

**Limite s → 0.**
- A caixa passa a ser fechada, s ∈ [−R/100; 0]. Em s = 0, X = p₀ para toda direção que entra, e as condições são contínuas.
- s = 0 não pertence ao modelo, porque ali a bola jovem chega à boca. Incluí-lo só aumenta o conjunto examinado, então uma exclusão completa da caixa fechada prova a faixa aberta.
- **Se sobrarem caixas indecididas em s = 0:** o sistema-limite tem 3 equações em 2 incógnitas (ν, α). Genericamente ele não tem solução, e o passo seguinte seria examiná-lo à parte, com exclusão e Krawczyk na face s = 0. Uma solução ali seria uma colisão degenerada na superfície da boca, e seria relatada, não descartada.

**Mecanismo 2.** Se a sondagem mostrar que ele persiste depois da reparametrização, troca-se o impulso por uma incógnita própria: u = v₁e + j·n, com j > 0. Isso torna V identicamente nulo. O lado j → 0 (movimento livre) seria então tratado como uma face da caixa.

**Equivalência.** É preciso um teste como `test_conditions_are_the_loop_quantities_of_billiards_up_to_positive_factors`: a função nova contra `billiards._loop_state`, em muitos pontos da faixa, a 10⁻³⁰.

**Controle positivo.** A reparametrização é uma máquina nova. Pela decisão (b) de 2026-11-08, ela precisa do seu próprio controle positivo: os mesmos 9 casos, agora no intervalo inteiro de s, que inclui a faixa. Não se conhece nenhuma solução dentro da faixa em nenhum caso. O controle mostra que a máquina nova recupera as soluções conhecidas, mas não testa a faixa em si contra uma solução conhecida, e isso precisa ser dito.

### 2.3 Resolução e custo estimado

- **Resolução.** Na parametrização por direção, a escala das condições não depende mais de 1/|s|. A largura mínima relativa atual (10⁻⁹) deve bastar. Isso é estimativa e precisa ser medido na primeira sondagem depois da implementação.
- **Custo da faixa.** O caso vazio inteiro usa hoje 656.220 subcaixas em 8 pedaços, cerca de 82 mil por pedaço. Espera-se que a faixa custe a mesma ordem de um pedaço, alguns minutos em 12 núcleos. Também é estimativa. O método atual, medido, não termina nem com 13M de caixas.
- **Custo do controle com a máquina nova.** O controle atual leva cerca de 100 min no `--full`. Se a máquina nova substituir a atual em todo o intervalo, o custo deve ser da mesma ordem. Se rodar só na faixa, o controle novo pode ficar restrito aos pedaços junto à boca dos 9 casos, com custo menor. As duas opções estão nas decisões abaixo.

### 2.4 Passos, depois de aprovado

1. Diagnóstico com o código atual: medir as larguras dos invólucros de V e do impulso nas caixas indecididas da sondagem, para separar os mecanismos 1 e 2.
2. `certify.conditions_dir`, com a desigualdade "o raio atinge a boca" e o teste de equivalência.
3. Sondagem de custo na faixa fechada [−R/100; 0], com relato.
4. Controle positivo da máquina nova, com Krawczyk nas caixas indecididas.
5. Integração no `--full` e redação do relatório. A frase de prova só se estende a [−3; 0) se a faixa for excluída por inteiro e o controle passar. Se não, o relatório diz o que ficou indecidido.

### 2.5 Decisões que cabem ao usuário

- **(a)** Aprovar a parametrização por direção como uma variante da máquina, com o seu próprio controle positivo.
- **(b)** Incluir s = 0 na caixa, como superconjunto.
- **(c)** Escolher se a máquina nova roda só na faixa ou substitui a atual em todo o intervalo, o que muda o custo do controle.
- **(d)** Definir o orçamento de tempo aceitável a mais no `--full`.

**Decisões do usuário (2026-11-09):**
- **(a)** Aprovada a parametrização pela direção de partida, como variante da máquina, com controle positivo **dentro** da faixa:
  - pelo menos três histórias autoconsistentes com a colisão em s ∈ [−0,00135; 0), por exemplo s = −1e-3, −1e-5 e −1e-7;
  - construídas a partir da história, verificadas pelo verificador independente;
  - critério: nenhuma excluída, e o Krawczyk a certifica;
  - se a construção não for possível em algum s, relatar por quê.
- **(b)** s = 0 incluído na caixa; o sistema-limite é examinado à parte, se sobrarem caixas ali.
- **(c)** A máquina nova roda na faixa, com verificação cruzada em [−0,01; −0,00135].
- **(d)** Um comando separado, `--certify`, roda a certificação, que sai do `--full`. O CLAUDE.md foi atualizado.

## 3. Implementação e resultados

Os números abaixo vêm das sondagens de 2026-11-09. O relatório oficial é `outputs/phase3_certify.md`.

### 3.1 A máquina por direção

O código está em `certify.conditions_dir` e `certify._reason_dir`, e a escolha de máquina é o parâmetro `machine="dir"` de `exclude`, `jacobian` e `krawczyk`. As fórmulas são as de 2.2, com duas correções que as sondagens exigiram:

- **Quina da direção tangente em s = 0.** Ali b → 0 e b² − c → 0, e a forma racionalizada de L₀ vira 0/0, com invólucro infinito.
  - Com ela, o pedaço que contém s = 0 não terminava: 21.900 caixas ficavam indecididas em 200 mil.
  - Em intervalos, L₀ passou a ser a interseção da forma racionalizada com −b − √(b² − c) e com [0; ∞). Isso vale porque as três envolvem o mesmo valor em todo ponto da família, onde b < 0 e b² ≥ c.
  - Com essa mudança, [−1e-6; 0] é excluída com 237 caixas.
- **Ponto exatamente em s = 0 na direção tangente.** A avaliação pontual usa −b − √(b² − c) = 0.

Os testes estão em `tests/test_certify_band.py`. Eles cobrem:
- a equivalência com `billiards._loop_state` a 10⁻³⁰ em três casos;
- a continuidade em s = 0;
- invólucros que contêm os valores pontuais, inclusive na quina tangente;
- o Krawczyk na raiz de controle;
- a separação entre `--full` e `--certify`;
- a frase de prova da faixa, que só sai sob as três condições.

### 3.2 A faixa do caso vazio e a verificação cruzada

- **Faixa:** s ∈ [−0,00135; 0], com s = 0 incluído, em 8 pedaços. A largura mínima é relativa à faixa inteira. São 3.132 subcaixas, todas excluídas, e nenhuma fica indecidida em s = 0. O sistema-limite não precisou ser examinado à parte.
- **Verificação cruzada** em [−0,01; −0,00135], na resolução da prova do caso vazio:

  | máquina | subcaixas | resultado | tempo |
  |---|---|---|---|
  | pela mira | 456.195 | todas excluídas | 171 s |
  | por direção | 321 (1.295 na sondagem anterior à interseção em L₀) | todas excluídas | 1 s |

  As duas concordam.

### 3.3 Controle positivo dentro da faixa: por que não há história física, e o que se usou

**A construção, a partir da história:**
1. Fixam-se a linha (θ, φ), B, R e s.
2. Resolvem-se miss = 0 e C = 0 em (ν, mira).
3. v₁ sai da equação V e é escrito como λ₁ decimal pela eq. 4.10.
4. A raiz é polida com esse λ₁ fixo e passa por `billiards.verify`.

**Tentativas de história física:**
- Foram 2 linhas com B = 0,3, R = 0,135 e mais 32 geometrias, com (θ, φ) entre 8 linhas e (B, R) entre 4 pares, em s = −1e-5.
- Em toda geometria em que há solução na faixa, as soluções vêm em pares, e nos dois a normal de contato é radial à boca 2. Para R ≪ B, vale |n·n_out| = 1,0000 com 4 algarismos.
- **Ramo com a bola velha fora da boca:** o impulso é < 0, ou seja, as bolas se afastam, e não há colisão.
- **Ramo com impulso > 0:** a bola velha e o ponto de contato ficam dentro da boca 2, o que o item (i) da p. 1081 exclui, e a bola velha entra na boca antes da colisão.
- O verificador independente rejeita todos esses casos.

**Por que acontece.** Com o centro da bola jovem a |s| ≪ R da superfície da boca, a bola jovem entra na boca logo após a colisão. O laço, pela travessia, só fecha com a normal de contato alinhada ao raio da boca em p₀; o ramo em que a bola velha empurra a jovem para dentro, vindo de fora, não aparece. Isso foi observado em todas as geometrias amostradas, e não está provado em geral.

**Não foi possível construir história física com a colisão na faixa**, em nenhum dos s pedidos (−1e-3, −1e-5, −1e-7).

**O que se usou como controle.** O ramo com impulso > 0 é raiz exata das equações da máquina (miss, C, V), com as três desigualdades dela satisfeitas: d·n_X < 0, along > 0 e impulso > 0. A máquina não usa o item (i). Esses casos testam a máquina, não a física, e o relatório diz isso.

Os três controles usam (θ; φ; B; R) = (0,5; 0,9; 0,3; 0,135), com os B e R do caso vazio. Em cada um, a máquina roda na faixa inteira, e o Krawczyk, nas caixas indecididas.

| s | λ₁ | subcaixas | em caixa excluída | Krawczyk na caixa da raiz | outras caixas indecididas |
|---|---|---|---|---|---|
| −1e-3 | 3,5848820852079515833 | 785.491 | 0 | raiz única | 4.090, todas "sem raiz" |
| −1e-5 | 3,8788614447816208674 | 617.151 | 0 | raiz única | 5.850, todas "sem raiz" |
| −1e-7 | 3,8821746427593530025 | 519.839 | 0 | raiz única | 3.572, todas "sem raiz" |

Nos três, λ₁ fica fora de [0; 1], o intervalo da eq. 4.10. O v₁ fica acima do intervalo do triângulo perigoso, o que é permitido para um controle da máquina.

**Observação sobre o sistema-limite.** Seguindo a raiz de controle de s = −1e-5 até −1e-7, ela converge suavemente: mira − mira₀ ∝ s, e v₁ → 3,8234. Há, portanto, raiz-limite em s = 0 para esse v₁. No caso vazio (λ₁ = 0,75), a faixa fechada foi excluída por inteiro, então não há raiz ali, nem em s = 0.

### 3.4 Custo

Medido nas sondagens:
- faixa: segundos;
- verificação cruzada: cerca de 3 min;
- controles da faixa: 5 a 7 min cada, em paralelo.

O custo dominante do `--certify` continua sendo o controle positivo do intervalo principal, que leva cerca de 90 min.
