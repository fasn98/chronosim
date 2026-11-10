# reviewer → claude · faixa junto à boca, máquina por direção, --certify (2026-11-09)

## Entregue
Revisão somente leitura do diff sem commit em `.claude/worktrees/fase3-faixa-boca` (branch `feat/fase3-faixa-boca`).
- **Arquivos lidos:** `chronosim/certify.py`, `chronosim/phase3.py`, `chronosim/run.py`, `CLAUDE.md`, `docs/phase3_certification.md`, `tests/test_certify.py`, `tests/test_certify_band.py`.
- **Contexto:** `chronosim/billiards.py` (`verify`, `history`) e `chronosim/provenance.py:40-51`.
- **Governam a tarefa:** `docs/phase3_borda.md` e as decisões (a)–(d) do usuário de 2026-11-09.
- **Rigor intervalar** de `conditions_dir`, `_meet`, `_nonneg`, `_reason_dir` e `alpha_range`: conferido, sem achado.

Achados:
- Nenhum bloqueante.
- Três importantes:
  - `phase3.py:1258`: a frase sobre o verificador estava fixa no código;
  - `certify.py`, em `krawczyk` e `jacobian`: b < 0 e b² − c > 0 não eram conferidos na máquina "dir";
  - `phase3.py:1254` e `provenance.py:48`: o rótulo CERT_ROOT_MODEL é dado a raízes que o verificador rejeita.
- Cinco notas: `phase3.py:1237`, `phase3.py:1264`, `phase3.py:1281`, `docs/phase3_certification.md`, e o desempenho e o teste de `band_control`.

Veredito: APROVADO COM RESSALVAS.

## O que falta
1. claude: frase do verificador a partir da execução, com teste.
2. claude: em "dir", exigir b < 0 e b² − c > 0 no Krawczyk, com teste.
3. Usuário: decidir o rótulo CERTIFIED dos controles da faixa.
4. claude, opcional: as notas 4 a 8.
5. Antes do merge: pytest, execução padrão, `--full` (o diff muda `full_tasks`, `report_section` e `_task r8r`) e `--certify`.
6. Conferir o resultado da faixa em `outputs/phase3_certify.md`.

## Riscos
- O reviewer não rodou nada, porque o hook é somente leitura.
- A caixa [0, ∞) em `_meet` só vale com c ≥ 0.
- A verificação cruzada só detecta caixas indecididas por engano. Uma exclusão por engano só seria detectada pelos controles, que usam outra geometria.
