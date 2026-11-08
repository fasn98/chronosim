# Handoff: reviewer → claude (artigo ReScience C, revisão curta 4)

Gravado por claude a partir do relato do reviewer, que é somente-leitura. **Veredito: BLOQUEADO: 1 achado bloqueante**, sobra do B2. B1 e B3 estão resolvidos.

## Entregue

- **B1 resolvido:** `PhaseOneTwoPiRel` = |ρ|·1440π²ξ⁴, que dá 2,82×10⁻¹⁹ em ξ = 0,5. As contas foram conferidas.
- **B2 quase todo resolvido:** o escopo agora é "listed in Table", com as exceções eq. 31 e eq. 52, e "declares neglected" ficou só para F3 fora do eixo.
- **B3 resolvido:** Visser gr-qc/9311026 está citado, e "cover none" ficou restrito ao que foi lido.
- **Notas aplicadas:** N1, N3, N4, N5, N6, N7, N8 e N9. N2 foi aplicada em parte: o arquivo `generated/` commitado continua com a contagem defasada em um.

## O que falta

- **BLOQUEANTE** (`metadata.yaml:53`; `content.tex:160`, `:267`): o texto afirma para F3 uma "identificação" do resíduo. Os outputs mostram queda linear sobre o eixo e concordância de lei de escala fora dele. A identificação analítica existe só em F4.
- **IMPORTANTE:** os metadados de `visser1994` não estão registrados no repositório.
- **Notas:**
  - o parâmetro pequeno da eq. 53 é a/b;
  - a eq. 31 não está no escopo de `:23`;
  - o resumo omite a ressalva da eq. 51;
  - o risco 3 de `NOTES.md` não menciona gr-qc/9311026;
  - não há teste da fórmula de `PhaseOneTwoPiRel`;
  - `F5_LIMIT_N` deveria ser `heuristic()`.

## Riscos

- O reviewer não executou nada.
- A renderização LaTeX não foi conferida.
