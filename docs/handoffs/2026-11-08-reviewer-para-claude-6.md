# Handoff · reviewer → claude (item B, certificação)

Gravado por claude a partir do relatório final do reviewer, que não tem ferramenta de escrita. Conteúdo resumido, sem alteração de veredito.

## Entregue

- Revisão somente leitura de `feat/fase3-certificacao` (sem commit) contra `main`. O reviewer não rodou python nem pytest; leu o código-fonte do mpmath.iv 1.4.1.
- **A prova é válida dentro do enunciado.** As condições são necessárias e os escalonamentos estão corretos (`certify.py:74-101`). A caixa cobre o espaço declarado, com pontas para fora e pedaços contíguos. Nenhum caminho exclui sem prova: NaN levanta, divisão por 0 dá (−inf, inf), `x**2` ≥ 0. As regras de CERTIFIED conferem com a decisão.
- **B1 (decisão do usuário):**
  - o rótulo "exclusão certificada" não foi aprovado;
  - a exclusão veio antes do Krawczyk (o plano, linha 13, a põe na etapa 2);
  - a redação "Provado…" toma o lugar de "nenhuma solução encontrada sob [critério]" (decisão 4).
- **I1** `provenance.py:59-67`: o construtor aceita CERT_PHYSICAL sem provas.
- **I2** `phase3.py:777-807`: a origem do relatório é texto fixo, não sai de uma `Quantity`.
- **I3** `phase3.py:790-806`: "Ω inteiro" não confere o caso, o n e o conjunto de pedaços.
- **I4** `billiards.py:191`, `phase3.py:829`: as pontas indefinidas na busca da mira são puladas sem registro, e o texto diz "nunca descartados em silêncio".
- **I5** `certify.py:46-47`: AIM_PAD e NU_PAD estão fora de `heuristic()`.
- **Notas:**
  - N1: "0 raízes certificadas" é fixo e não diz que o Krawczyk não foi implementado.
  - N2: o valor da `Quantity` CERTIFIED é o número de caixas.
  - N3: a restrição de família deveria estar na oração principal.
  - N4: silêncios anteriores ao branch: deduplicação só por s em `phase3.py:82-83`, arco encolhido em `billiards.py:307`, zero num nó contado como negativo em `billiards.py:481`.
  - N5: `refine_nu` registra duas vezes os colchetes das pontas originais.
  - N6: faltam as hipóteses "não tangente" e "impulso > 0" na lista do documento.

## O que falta

1. Decisão do usuário sobre B1, registrada em `docs/phase3_plan.md`.
2. Corrigir I1–I5. Testes do scalpel para I1, I3 e I4.
3. Rodar pytest e `--full` depois das correções.

## Riscos

- A prova depende do mpmath.iv 1.4.1 e da leitura do modelo de EKT em `billiards.py`.
- Os números do run não foram conferidos pelo reviewer.

Veredito: **BLOQUEADO, 1 bloqueante (B1, decisão pendente do usuário; nenhum defeito técnico na prova).**
