# claude → scalpel · faixa junto à boca, máquina por direção, --certify (2026-11-09)

## Entregue
Mesmo diff de `docs/handoffs/2026-11-09-claude-para-reviewer-2.md` (worktree `.claude/worktrees/fase3-faixa-boca`). Os testes atuais estão em `tests/test_certify_band.py`.

## O que falta
1. Escrever testes num arquivo novo, `tests/test_certify_band_scalpel.py`, sem alterar produção nem testes existentes. Cobrir:
   - `conditions_dir` em iv, envolvendo os valores pontuais perto da quina tangente em s = 0, em mais caixas que as atuais;
   - `_meet` com intervalos disjuntos;
   - `_reason_dir` em cada motivo;
   - `krawczyk(machine="dir")` numa caixa sem raiz;
   - `exclude(machine="dir")` com `watch`, sem normalizar;
   - `start(mode="certify")` com a lista de tarefas;
   - `run.main(["--certify"])` com saída para um diretório temporário, se for rápido o bastante (provavelmente não; nesse caso, simular).
2. Rodar `python -m pytest` inteiro com o `.venv` e relatar a saída.

## Riscos
- Há um `--certify` rodando em segundo plano (cerca de 100 min). Não matar.
