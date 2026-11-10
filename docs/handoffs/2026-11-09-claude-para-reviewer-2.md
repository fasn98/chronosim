# claude → reviewer · faixa junto à boca, máquina por direção, --certify (2026-11-09)

## Entregue
Worktree `.claude/worktrees/fase3-faixa-boca`, branch `feat/fase3-faixa-boca`. O diff ainda não tem commit.

Arquivos:
- `chronosim/certify.py`:
  - `conditions_dir`, `_reason_dir`, `_meet`, `_nonneg`, `alpha_range`, `MACHINES`;
  - `exclude`, `jacobian` e `krawczyk` ganharam o parâmetro `machine`;
  - `HYPOTHESES_DIR`;
  - `pc` em `case_inputs`.
- `chronosim/phase3.py`:
  - `BAND_CONTROLS`, `band_range`, `band_scale`, `band_chunk`, `cross_check`, `band_control`, `_alpha_of`;
  - `certify_tasks`, `start(mode=)`, `r8rc`;
  - `_band_lines`, `certify_section`;
  - o cabeçalho do relatório do `--full`.
- `chronosim/run.py`: `--certify` e `certify()`.
- `CLAUDE.md`, `docs/phase3_certification.md`.
- Testes: `tests/test_certify_band.py` (novo) e `tests/test_certify.py` (dois testes atualizados para a decisão (d)).

Decisões: `docs/phase3_borda.md`, mais as decisões (a)–(d) de 2026-11-09 do usuário.
- (a) Máquina por direção, com controle positivo dentro da faixa: três histórias construídas a partir da história.
- (b) s = 0 incluído na caixa.
- (c) Verificação cruzada em [−0,01; −0,00135].
- (d) `--certify` separado do `--full`.

## O que falta
Revisar, somente leitura. Pontos a examinar:
1. **Rigor intervalar de `conditions_dir`.** Interessam três coisas:
   - a forma de c = s(2e·pc + s), que supõe |pc| = B;
   - a interseção `_meet` das duas formas de L0 e de [0, ∞);
   - o `_nonneg` em b² − c.

   Cada passo precisa envolver o valor verdadeiro em todo ponto da família, ou seja, onde b < 0 e b² ≥ c.
2. **`_reason_dir`.** As exclusões "b ≥ 0", "raio não atinge" e "P dentro da boca" estão corretas?
3. **Controles da faixa.** O verificador rejeita os três pelo item (i). Eles são raízes da máquina, não histórias físicas. O relatório diz isso sem exagero?
4. **Frase de prova da faixa** (`_band_lines`). Ela só sai com exclusão completa, verificação cruzada concordante e controle aprovado. A cobertura de [−3; 0) só é afirmada com a prova principal presente.
5. **Separação `--full`/`--certify`.** Nada de certificação continua no `--full`, e o relatório não perde informação.

## Riscos
- A primeira versão da forma racionalizada de L0 era 0/0 na direção tangente em s = 0, e o pedaço que contém s = 0 não terminava. A correção é a interseção com −b − √(b² − c).
- Os controles não são físicos. Não se construiu nenhuma história física na faixa: em todas as geometrias amostradas, a solução tem normal de contato radial e cai num destes casos:
  - bola velha fora, com impulso < 0;
  - bola velha dentro da boca, com impulso > 0.
