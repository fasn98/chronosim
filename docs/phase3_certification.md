# Fase 3 · certificação intervalar (exclusão)

Este documento registra o enunciado, a função, a caixa e o critério do teste intervalar, como pede o plano (`docs/phase3_plan.md`, decisão 2 da primeira rodada).

**Onde roda.** Desde 2026-11-09, por decisão (d) do usuário, toda a certificação roda num comando à parte: `python -m chronosim.run --certify`, que escreve `outputs/phase3_certify.md`. O `--full` não a roda mais. A faixa junto à boca, com a máquina por direção, a verificação cruzada e o controle dentro da faixa, está em `docs/phase3_borda.md`. O código está em `chronosim/certify.py`, e os testes em `tests/test_certify.py` e `tests/test_certified_origin.py`.

## Origem CERTIFIED

- **Ordem:** CITED < DERIVED < CERTIFIED < COMPUTED < HEURISTIC, aprovada pelo usuário em 2026-11-07.
- **Rótulo obrigatório:** toda quantidade CERTIFIED tem um rótulo, e só ela tem.
  - "raiz certificada", quando só o teste de Krawczyk passou. Leva sempre o nome do sistema: "raiz certificada das eqs. 4.2–4.3" ou "raiz certificada do modelo geométrico (billiards.py)" (decisão de 2026-11-09). Um rótulo "raiz certificada" sem o nome do sistema é recusado, e há teste para isso. Uma raiz certificada não é "solução física certificada".
  - "solução física certificada": `certified()` recusa o rótulo sem as duas provas intervalares, (a) a, a′, c > 0 e (b) ausência de colisões extras (aprovado). Sem (b), `physical_solution()` devolve COMPUTED, com a lista do que foi e do que não foi provado.
  - "exclusão certificada": nenhuma raiz das equações declaradas na caixa declarada. Aprovado em 2026-11-08 como rótulo próprio, distinto dos outros dois.
- **Combinação:** `combine()` nunca devolve CERTIFIED. Com COMPUTED dá COMPUTED, como decidido. Com CITED ou DERIVED também dá COMPUTED, porque a aritmética de `combine()` não faz parte de uma prova intervalar. Essa segunda regra é escolha deste branch, mais restritiva que a decisão.
- **Sem promoção:** `degrade()` não leva nada para CERTIFIED, nem a partir de DERIVED. Só `certified()` cria a origem.

## O teste de exclusão

**Família examinada.** Uma autocolisão, no modelo geométrico de `billiards.py` (EKT 1991, Sec. II A, p. 1081, Fig. 2):

- A bola jovem vem pela linha inicial da Fig. 11, p₀ + s e (`phase3.line_theta_phi`). Ela colide em P = p₀ + s e, com s no intervalo declarado e antes de chegar à boca 2 (s < 0).
- Depois da colisão, vai em linha reta até a boca 2, faz **exatamente uma** travessia (sai da boca 1 uma unidade de tempo antes) e volta como bola velha. No instante da colisão, a bola velha está em Q = P − 2R n, com n = (cos ν, sin ν).

**Incógnitas.** (s, ν, mira); a mira é o ângulo polar do ponto de entrada X = C₂ + B (cos mira, sin mira).

**Função** (`certify.conditions`). São as mesmas condições de `billiards._loop_state` e `billiards.solve_mp`, multiplicadas por fatores positivos:

| condição | significado |
|---|---|
| miss = cross(w, Q − E) = 0 | o laço passa por Q |
| C = (v₂ − u)·n = 0 | as componentes normais são trocadas |
| V = cross(u, n) − v₁ cross(e, n) = 0 | a velocidade de chegada é o v₁ da eq. 4.10 |
| d·n_X < 0 | o raio entra na boca 2 em X, transversalmente |
| along = (Q − E)·w > 0 | a bola velha já saiu da boca 1 na colisão |
| kick = (u − v₁ e)·n > 0 | o impulso afasta as bolas |

As três desigualdades são condições necessárias de uma solução do modelo. As outras condições do verificador (a bola não passa por outra boca antes, nenhuma sobreposição ao longo da história) não entram: omiti-las só aumenta o conjunto examinado.

**Caixa.**
- s: o intervalo da busca, [−S_FAR; −R·`ekt_s_near_over_R`], arredondado para fora.
- ν ∈ [−0,01; 2π + 0,01] e mira ∈ [−3,15; 3,15]: superconjuntos do círculo. As miras invisíveis de P são excluídas pela desigualdade d·n_X < 0.
- Entradas: θ, φ, λ₁, B e R são as strings decimais do caso. p₀, e e v₁ são calculados em intervalos a partir delas, pelas eqs. 4.8 e 4.10.

**Critério.**
- Uma caixa está excluída quando o invólucro intervalar de miss, C ou V não contém 0, ou quando uma das desigualdades falha na caixa inteira.
- Uma caixa que não é excluída é dividida ao meio no lado mais largo, medido em relação ao lado da caixa raiz.
- A divisão para quando a caixa é excluída ou fica mais estreita que `ekt_cert_min_width`. O que sobra é relatado como indecidido, nunca descartado.
- Se não sobra nada, nenhum ponto da caixa raiz satisfaz as condições.
- Aritmética: mpmath.iv a `ekt_cert_prec` = 53 bits, com arredondamento para fora.

**Em que a prova confia.**
1. A correção do `mpmath.iv`, versão 1.4.1, fixada em `requirements.txt`.
2. A função `conditions` reescreve o modelo de `billiards.py` sem erro. O teste `test_conditions_are_the_loop_quantities_of_billiards_up_to_positive_factors` compara as duas em mais de 50 pontos, para dois casos, a 10⁻³⁰.
3. A leitura do modelo de EKT feita em `billiards.py`, a mesma de toda a fase 3.
4. Duas hipóteses que excluem casos de fronteira:
   - a entrada na boca 2 não é tangente (d·n_X < 0 estrito);
   - o impulso é estritamente positivo (kick > 0; impulso nulo é o movimento livre, que neste caso se sobrepõe a si mesmo).

   A terceira desigualdade, along > 0, é a mesma de `billiards.loop`.

**O teste não é vácuo.** Uma caixa de ±10⁻⁴ em torno de uma solução verificada do caso central não é excluída, e a caixa indecidida que sobra contém a solução. Uma caixa longe dela é excluída.

**O que não está coberto:**
- duas ou mais travessias entre as visitas;
- colisão depois de uma travessia;
- s fora do intervalo. A faixa [−R/100; 0) junto à boca é examinada à parte, pela máquina por direção (`docs/phase3_borda.md`); a cobertura depende do resultado em `outputs/phase3_certify.md`;
- várias colisões.

Para essas famílias vale só a busca por bracketing, cuja origem é COMPUTED.

**Krawczyk** (`certify.krawczyk`): K(X) = x̃ − Y·F(x̃) + (I − Y·J(X))·(X − x̃). J(X) é envolvido por diferenciação automática em intervalos, e Y é o inverso do jacobiano no ponto médio. Se K ⊂ int X, há uma única raiz em X; se K ∩ X = ∅, não há nenhuma. No caso vazio, nenhuma caixa ficou indecidida, e não houve raiz a certificar.

**Controle positivo** (decisão do usuário de 2026-11-08). A mesma máquina, com o mesmo código, a mesma precisão e as mesmas hipóteses, roda sobre os casos com soluções conhecidas no mesmo submodelo: R8 #0, #1, #3, #4 e #5, e a Fig. 7(b) com B = 0,04, 0,02, 0,01 e 0,005. Isso é feito por `phase3.control_piece` e `control_summary`, na execução `--full`.
- Critério: nenhuma caixa que contenha uma solução conhecida pode ser excluída, nem uma cópia periódica dela (ν ± 2π, mira ± 2π) que caia na caixa. Cada cópia é vigiada, mas a solução é contada uma vez.
- **Resolução:** o espaço de cada caso é dividido em pedaços de s, que só servem para paralelizar e relatar. A largura mínima das caixas é medida contra 1/`ekt_cert_chunks` do intervalo de s e contra ν e mira inteiros, a mesma escala dos pedaços do caso vazio (parâmetro `scale` de `exclude`).
- **Orçamento:** pedaços com solução conhecida param em `ekt_cert_control_sol_budget` caixas, e os outros em `ekt_cert_control_budget`. O que sobra é relatado como não examinado, com a fração do volume. `ekt_cert_control_sol_budget` vale 20M desde 2026-11-09, o dobro das 9.963.011 caixas que o pedaço com solução do R8 #3 usou sob o limite anterior de 10M. É uma margem contra mudanças futuras do método, não um ajuste de resultado: qualquer orçamento que complete o pedaço dá o mesmo resultado.
- **Krawczyk:** roda nas caixas indecididas que a máquina examinou. Uma solução certificada em todas as suas cópias recebe CERTIFIED com o rótulo "raiz certificada do modelo geométrico (billiards.py)".

O resultado está em `outputs/phase3_full.md`.
