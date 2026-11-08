# Handoff: reviewer → claude (fase 3, reprodução R1–R11)

Gravado por claude a partir do relato do reviewer (somente-leitura). Veredito: **BLOQUEADO: 2 achados bloqueantes**, mais 8 notas.

## Entregue

- **Escopo:** revisão somente-leitura de `main..e7bd1d8` (commits 7712bfb e e7bd1d8) na branch `feat/fase3-bilhar`, lida com `git show e7bd1d8:<caminho>` e `git diff main..e7bd1d8`.
- **Conferido sem achado:**
  - **Transcrição:** `ekt.py` contra `docs/equation_transcription.md`, termo a termo. Conferem 2.1b, 3.1, 3.2, os instantes da p. 1084, 3.6–3.8, 4.2–4.4, 4.8, 4.10, 4.12–4.15, B3–B11, B13–B15, B20 (as duas formas), B30–B32, A6, A7, A18 e A19.
  - **B20:** a razão 4 entre os termos em B das duas formas é a identidade tan²φ − tan²θ = sin(φ−θ) sin(φ+θ)/(cos²φ cos²θ), aplicada ao que está transcrito.
  - **Travessia e itinerário:** `billiards.traverse` (`billiards.py:105-109`), `first_entry` com `skip`, `itinerary_windows` (`:290-324`) e `_gap` (`:241-247`).
  - **Condições do solver:** tempo líquido nulo, troca normal (v2 − u)·n = 0, v1 por paralelismo e sinal de `kick`.
  - **Verificador (`billiards.py:613-670`):** reconstrói a história de forma independente, testa todos os pares de segmentos com sobreposição temporal e isenta só os três pares contíguos em t = 0. Colisões extras são rejeitadas.
  - **Caso vazio:** "distância mínima 0,135" conferida à mão. Com θ = 0, a separação em x vale |2R(2λ₁ − 1)| = R.
  - **Filtro das eqs. 4.2–4.3:** eq. 4.6 e distâncias do Apêndice A, como descrito.
  - **Precisão:** nenhum `mp.dps` global fora de `_init_worker`. As grades de `r8x` são restauradas em `finally`. As bisseções têm guarda.
  - **Anti-RNG:** limpo. **Linguagem:** limpa. **Origens:** COMPUTED para o modelo; DERIVED para as curvas transcritas na figura.

### BLOQUEANTE B1: o critério de "nenhuma solução encontrada" está incompleto, e uma das afirmações dele é falsa

- **Onde:** `phase3.py:530`, `:618-624`; `billiards.py:391-392`, `:405-407`, `:410-412`, `:519-521`, `:528-529`.
- **O quê:**
  - "São contados, nunca descartados" só vale para o eixo s.
  - Em `middle_roots`, intervalos de ν com mudança de contagem são pulados com `continue`.
  - Bisseções que perdem a raiz e cruzamentos rejeitados como polo somem sem registro.
  - Os intervalos de s ([−3, −R/100] em `dangerous`; [10⁻⁶, 3] em `after_traversal`) não estão declarados.
  - "Sem candidato geométrico" é texto fixo.
- **Correção:** contar e propagar esses descartes, declarar o intervalo de s e condicionar o texto aos dados.

### BLOQUEANTE B2: o relatório afirma verificações que a execução não faz

- **Onde:** `phase3.py:609-610`, `:683`.
- **O quê, parte 1:** `sign_rule` é calculado mas nunca relatado. A eq. A6 só é comparada em teste, e só no caso central. "sign(d) = s" é tautológico.
- **O quê, parte 2:** "Estabelece" afirma limite B → 0 para B32, 4.14 e 4.15, que foram medidas num único B = 0,01.
- **Correção:** acrescentar às linhas de R8 o desvio de A6 e `sign_rule` e condicionar o texto a eles. Restringir "limite B → 0" a 3.1, 3.2 e 4.12–4.13.

### Notas

- **N1:** constantes de domínio e tolerância sem `heuristic()` (lista no relato, `billiards.py` e `phase3.py`).
- **N2:** `ekt._q` rotula DERIVED mesmo quando recebe valor COMPUTED (`phase3.py:241`, `:278`, `:330`). Não chega à saída, mas quebra a garantia de não promoção na API.
- **N3:** `solutions` descarta os dígitos de `backend_for` (`phase3.py:59`), e `_circ` usa `math.pi` (`billiards.py:456-458`). Latente: nenhuma tarefa atual entra no modo mpmath.
- **N4:** texto fixo no gerador do relatório (`phase3.py:607-608`, `:654-657`, `:662`, `:684-685`) e quebras com listas vazias (`_order`, `max(near)`, `min(offs)`, `min(found_n)`).
- **N5:** `r8x` está fixo em `R8_CASES[2]`, mas o relatório o imprime sob qualquer caso vazio.
- **N6:** o domínio 0 < α < π, |β| < π/2 de `eq42_43_roots` não está declarado.
- **N7:** a decisão sobre a eq. B30 (`docs/equation_transcription.md:273`) não é reproduzível a partir do repositório. Um teste de ordem em B discriminaria.
- **N8:** o título "Trechos do texto que parecem trocados" (`docs/equation_transcription.md:287`) deveria ser "Trechos do texto com inconsistência aparente".

## O que falta

1. Corrigir B1 e B2.
2. N1, N3 e N4 merecem tarefa antes do mapa.
3. Rodar `pytest` e `run`, e mandar para uma revisão curta só de B1 e B2.
4. Pendências humanas sem mudança: B20 no scan e citação Q12 (p. 1088).

## Riscos

- O reviewer não executou testes (hook somente-leitura).
- O scan não foi conferido.
- B1 é risco de completude, não defeito demonstrado: o código não registrava os descartes.
