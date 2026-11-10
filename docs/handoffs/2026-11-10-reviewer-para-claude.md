# reviewer → claude · 2026-11-10

(Gravado por claude: o reviewer é somente-leitura e não tem Write. Conteúdo do relatório dele, resumido sem alteração de achados.)

## Entregue
Revisão somente-leitura do diff sobre `de09722` em `feat/fase3-faixa-boca`. Veredito: **APROVADO COM RESSALVAS**, nada bloqueante.
- Rótulo CERT_ROOT_DIR correto (`chronosim/provenance.py:52`, `:55`; `chronosim/phase3.py:1263`, `:1268`); CERT_ROOT_MODEL só no intervalo principal; recusa de "raiz certificada" sem sistema intacta.
- Seção 4 de `docs/phase3_borda.md`: citação do item (i) confere com o PDF (p. 1081); leitura de `billiards.py:734-737` e as desigualdades conferem; nem `certify.conditions` nem `certify.conditions_dir` usam o item (i).
- Testes cobrem o rótulo novo.

## O que falta
- IMPORTANTE: números da seção 4 (tabela 4.2, grade 4.3, ponto de exemplo) sem script/teste que os reproduza; 21.120 ≠ 3 × 120 × 60 = 21.600 sem explicação.
- NOTA: `docs/phase3_borda.md:241` (dizer que 0,225/0,234 são da linha θ = 0; "implicada pela primeira"; margem da boca 1 contra 2R).
- NOTA: nomes de teste desatualizados (`tests/test_certify_scalpel_2.py:227`, `tests/test_certified_origin.py:23`); guarda frágil `"billiards.py)" not in text` em `tests/test_certify_band.py:168`.

## Riscos
- O reviewer não rodou nada (hook nega python); grade e ponto de exemplo não verificados por ele.
