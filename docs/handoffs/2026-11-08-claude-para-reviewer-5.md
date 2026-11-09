# Handoff · claude → reviewer (revisão curta de I1–I5 do item B)

Mesma worktree e branch: `feat/fase3-certificacao`, ainda sem commit.

## Entregue

- **I1** (`chronosim/provenance.py`, `Quantity.__post_init__` e o campo `proofs`): o construtor recusa CERT_PHYSICAL sem as duas provas. Teste: `test_the_constructor_cannot_bypass_the_physical_proofs`.
- **I2 e I3** (`chronosim/phase3.py`, `_certification_lines` e `certify_chunk`):
  - a origem e o rótulo do relatório agora saem de `certify.exclusion_quantity`;
  - a palavra "Provado" só aparece com todos os pedaços, do mesmo caso e cobrindo [lo, hi] sem lacuna;
  - testes em `tests/test_certify.py`, `test_no_proof_wording_when_…`.
- **I4:** ressalvas na docstring de `unresolved_fraction` e nos dois textos do relatório. Os intervalos da mira com uma ponta indefinida não são registrados; o texto diz isso.
- **I5:** `ekt_cert_aim_pad` e `ekt_cert_nu_pad` passaram a ser `heuristic()`.
- **N1:** a linha dos itens diz que o Krawczyk não foi implementado.
- **N2:** o valor da `Quantity` CERTIFIED é 0 (soluções na caixa), com as subcaixas no `source`.
- **N5:** `refine_nu` não repete os descartes das pontas originais.
- **N6:** a lista de hipóteses foi completada em `docs/phase3_certification.md`.
- **Caixa com lado nulo:** `exclude` recusa (nota do scalpel).
- **B1:** o pedido B do usuário está registrado em `docs/phase3_plan.md`, na seção "terceira rodada". O rótulo e a redação sobre Ω continuam pendentes do usuário.
- **N3** não foi aplicada: a redação é a ditada pelo usuário, e a restrição de família fica no colchete e na linha "Não coberto pela prova", logo abaixo.
- **N4** não foi aplicada: são silêncios anteriores a este branch, fora do caminho certificado.

## O que falta

Confirmar que I1–I5 estão resolvidos, sem achado novo inventado. Os testes do scalpel estão em `tests/test_certify_scalpel.py`.

## Riscos

O pytest, o `--full` e a execução padrão estão rodando agora; os números ainda não saíram.
