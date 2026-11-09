# Handoff · reviewer → claude (revisão curta de I1–I5, item B)

Gravado por claude a partir do relatório final do reviewer, que não tem ferramenta de escrita. O conteúdo vai resumido e o veredito não muda.

## Entregue

- **Resolvidos, verificados no código:**
  - I1 (`provenance.py:73-74`);
  - I2 (`phase3.py:808-818`);
  - I3 (`phase3.py:792-796`, `:279`);
  - I4 (`phase3.py:139-141`, `:839`, `:969`; `billiards.py:186`);
  - I5 (`certify.py:46-47`);
  - N1, N2, N5 e N6;
  - a caixa de lado nulo (`certify.py:128-129`).
- **N3 e N4 sem aplicação:** as justificativas foram aceitas.
- **B1 continua pendente do usuário:** (a) o rótulo "exclusão certificada" e (b) a redação "provado…" apoiada em Ω inteiro (`docs/phase3_plan.md:43-46`).
- **I6, novo:** `phase3.py:827`. A seção abre com "Todas as grandezas desta seção são COMPUTED", mas o caso vazio imprime "Origem: CERTIFIED".
- O reviewer não rodou python e não leu `tests/test_certify_scalpel.py`.

## O que falta

1. O usuário decide B1.
2. Claude corrige I6.
3. Os resultados reais do pytest, do `--full` e da execução padrão entram no relato.

## Riscos

- Os números do caso vazio não foram conferidos pelo reviewer.
- A prova depende do mpmath.iv 1.4.1, do modelo de `billiards.py` e das duas hipóteses de fronteira.

Veredito: **BLOQUEADO, 1 bloqueante (B1, decisão do usuário); I1–I5 resolvidos; I6 novo.**

Nota de claude: I6 foi corrigido em `phase3.py:827`. O teste está em `tests/test_certify.py`, em `test_the_section_header_does_not_deny_the_certified_line`.
