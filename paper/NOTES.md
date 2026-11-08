# Notas de submissão: ReScience C

Consultado em 2026-11-08. Só está aqui o que as páginas dizem; cada ponto traz a URL de onde veio. Onde uma página não carregou, isso está dito.

## Fontes consultadas

| URL | Estado |
|---|---|
| https://rescience.github.io/ | lida |
| https://rescience.github.io/faq/ | lida |
| https://rescience.github.io/write/ ("Overview of the submission process") | lida |
| https://rescience.github.io/edit/ (guia de editores e revisores) | lida |
| https://rescience.github.io/edition/, /submit/, /about/ | **404**; não há instrução nessas URLs |
| https://github.com/ReScience/template (README, `metadata.yaml`, `content.tex`, `Makefile`, `article.tex`, `rescience.cls`) | lidos, branch `master` |
| https://github.com/ReScience/submissions (README e `.github/ISSUE_TEMPLATE/*.md`) | lidos |
| https://github.com/rescience/submissions/issues/new/choose | só a tela de login do GitHub; os modelos foram lidos direto do repositório |
| ReScience/submissions#143 (optosim), ReScience/ReScience#121 e #27 | lidos |
| fasn98/optosim (`metadata.yaml`, `optosim_rescience.pdf`) | lidos, como referência de estrutura e tom |

As leituras foram feitas por subagentes de busca na web. Os trechos entre aspas são citações literais relatadas por eles; não conferi cada um contra a página.

## Tipo de artigo e título

- **Tipos de artigo** (FAQ):
  - "Re means replication. Rp means reproduction. Either one can be negated in case of failure to give: ¬Re and ¬Rp."
  - Replicação é "repeating a published protocol, respecting its spirit and intentions but varying the technical details".
  - Reprodução é "running the same computation on the same input data".
- **Prefixo do título** (comentário do `metadata.yaml` do template): "[Re]" para replicação bem-sucedida e "[¬Re]" para replicação que falhou.
- **Replicação parcial** (/edit/): "A full replication covers all the results of the original work, whereas a partial replication covers only a subset of the results." O modelo de issue `replication.md` diz "success or failure, full or partial".
- **Resultados negativos** (/write/): "ReScience C accepts negative results as well".
- **Tipo nos metadados:** o `metadata.yaml` aceita `Editorial`, `Letter` ou `Replication`; este artigo é `Replication`.

## Estrutura e extensão

- **Nenhuma página exige seções nem limite de tamanho.**
- O `content.tex` do template é o corpo de um editorial de 2018, sem modelo de seções.
- **O que /write/ pede:**
  - "the main focus should be on the replicated results and on any obstacles encountered during the replication."
  - O artigo deve ser "sufficiently self-contained, but not to repeat all aspects of the original article or your documentation".
  - É útil um resumo de alto nível do modelo original, porque o original pode estar "behind a paywall". O PRD 43, 3929 é de acesso fechado.
- **O que /edit/ pede:** os autores "should clearly state why they think they have replicated the paper (same figures, same graphics, same behavior, etc.) and explain any obstacles".
- **Formato:** o template LaTeX é o preferido (FAQ). Outro formato é aceito "if the result is comparable to the PDF produced by our template".
- **Referência de tom:** o artigo do optosim (#143) tem 5 páginas, com as seções Introduction, Methods, Results, "What this does and does not establish", Conclusion e "AI assistance".

## Metadados (`metadata.yaml`, do template)

- **O autor preenche na submissão:**
  - `title`, `authors` (name, orcid, email, affiliations) e `affiliations` (code, name, address). O código `*` marca o autor de contato e é exigido mesmo com autor único.
  - `keywords`.
  - `code`: `url` (obrigatório para replicação), `doi` e `swh`. Comentário do template: "url is mandatory for replication, doi after acceptance".
  - `data`: url e doi, opcional.
  - `replication`: `cite`, `bib` (chave no `.bib`), `url` ("try to link to a non-paywall version") e `doi`.
  - `abstract`, entre aspas duplas.
  - `bibliography`, `type`, `domain` (um só, "not overly specific") e `language`.
- **O autor preenche depois do aceite:** `review.url` e `contributors` (editor e revisores).
- **O editor fornece:** `dates`, `article` (número, DOI do Zenodo e URL) e `journal`.
- **Usados no optosim (#143), como referência:** `affiliations: Independent researcher, São Paulo, Brazil`; `type: Replication`; `language: Python`.

## Código, licença e arquivamento

- **Código público** (/write/): "Upload your code to a public repository (e.g. GitHub)".
- **DOI no Zenodo** (/write/): depois do aceite, "upload your code to Zenodo in order to get a DOI".
- **Software Heritage:** nenhuma página exige; aparece só como campo `swh` opcional.
- **Licença aberta** (FAQ): obrigatória, "entirely under open licenses". O editor confere "that there is a license in the code repository" (/edit/).
- **Licença do artigo:** CC BY (/write/).
- **Linguagem:** "We favor a Python implementation" (FAQ).
- **Critérios de publicação** (/write/): "Rigorous methodology", "Original source code" e "Substantial evidence for replication of the original results or explanations why original results cannot be replicated".
- **Critério dos revisores** (/edit/): eles "must be able to run the proposed implementation on their computer, and obtain the same results". Também: "Uncommented or obfuscated code is as bad as no code at all."

## Compilação do template

- **Comando:** `make`. O alvo gera `metadata.tex` a partir de `metadata.yaml` com `./yaml-to-latex.py`, depois roda `latexmk` com **XeLaTeX** e **biber** (biblatex, estilo numérico).
- **Dependências:** Python 3 com PyYAML, e as pastas de fontes do template (`source-serif-pro/`, `roboto/`, `source-code-pro/`) ao lado dos fontes.
- **Bibliografia:** `bibliography.bib`.
- **Faixa "UNDER REVIEW":** aparece no PDF enquanto o DOI do artigo está vazio.

## Submissão

- **Forma** (/write/): uma issue em https://github.com/rescience/submissions/issues/new/choose, revisada em público.
- **Campos do modelo `replication.md`:** "Original article", "PDF URL", "Metadata URL", "Code URL", "Scientific domain", "Programming language" e "Suggested editor".
- **Checklist do README de submissions:** URL do PDF, URL do YAML, URL do código ("mandatory for replication"), dados (opcional) e sugestão de editor e revisores do board (opcional).
- **Processo** (/write/, /edit/):
  - Um editor aceita a submissão ou a rejeita, com justificativa, e convida pelo menos dois revisores na mesma issue.
  - Os revisores clonam o repositório e precisam obter os mesmos resultados.
  - Se os dois concordam, o artigo é aceito.

## Regras sobre o artigo original

- **Publicado:** o original precisa estar publicado ("already published research", página inicial). O PRD 43, 3929 (1991) está.
- **Trabalho próprio:** não se pode replicar trabalho próprio nem de colaboradores próximos (FAQ; /write/). Não se aplica aqui.
- **Já replicado** (FAQ): "If the original research has already been replicated previously (in ReScience C or elsewhere), it is unlikely we'll consider another replication". Ver o risco 3 abaixo.
- **Contato com os autores originais**, quando se alega falha:
  - /write/: "expected ... to have attempted to contact the authors of the original article."
  - /edit/: contato "strongly encouraged", mas "not a requirement".
- **Idioma:** nenhuma página traz regra de idioma.

## IA

- **Política:** nenhuma página do site fala de uso de IA.
- **ReScience/ReScience#121** ("AI submissions + AI reviews", aberta por Nicolas Rougier em 3 de novembro de 2026; aberta e sem decisão):
  - Ele diz que recebe muitas submissões com "AI code/AI prose".
  - Diz que vem "asking author to review at least two submissions".
  - Diz que algumas submissões foram "postponed" quando um mesmo autor enviou várias em poucos dias.
  - As opções que ele lista vão de "Pre-submission with desk reject" a "We asked for some human work proof".
- **#143 (optosim):**
  - Rótulos: `AI: Code`, `AI: Prose`.
  - O editor pediu duas revisões de outras submissões antes de prosseguir.
  - O primeiro revisor questionou (ponto 1) se uma re-implementação conta como replicação; o ponto ficou para o editor decidir.
  - Ainda não há segundo revisor confirmado. Houve uma oferta em 6 de novembro.

## Riscos e decisões que são do autor

1. **Segunda submissão com a #143 em andamento.** Pelo que o editor relatou em #121, várias submissões do mesmo autor levam ao adiamento. Ele também pede duas revisões por autor. Isso define o momento da submissão; não está nas regras publicadas.
2. **Escopo "computational research".** Kim & Thorne (1991) é sobretudo analítico. O que se replica são as derivações e avaliações numéricas das fórmulas fechadas a partir das expressões de partida do próprio artigo. O artigo precisa dizer isso com clareza, porque o ponto 1 do revisor em #143 levanta a mesma questão.
3. **"Já replicado."** Klinkhammer (1992) e Visser (1993) recalculam a divergência perto da geodésica nula fechada (o análogo da eq. 67) e declaram concordância. Um editor pode ler isso como replicação prévia. O argumento contrário é que os dois obtêm metade da eq. 67 e não reproduzem F1–F7. Isso precisa ser enfrentado no texto, não omitido.
4. **[Re] e o fator 2.**
   - Os valores impressos são reproduzidos pelo operador como impresso, e o fator 2 é uma inconsistência de convenção diante dos árbitros.
   - Isso sustenta "[Re]" com discrepâncias documentadas. O título pedido já usa "[Re]".
   - Se o editor achar que é falha parcial, o prefixo seria "[¬Re]". A decisão é do autor e do editor.
5. **Contato com Kim e Thorne.** Não é exigido quando não se alega falha, mas /write/ e /edit/ o recomendam ao apontar discrepâncias. Decisão do autor; nenhum agente faz esse contato.
6. **Repositório público e licença.**
   - Em 2026-11-08, `https://api.github.com/repos/fasn98/chronosim` devolve "Not Found": o repositório está privado ou não foi publicado.
   - Não há arquivo `LICENSE` na raiz.
   - As duas coisas são obrigatórias antes da submissão.
7. **Compilação local.** `xelatex`, `latexmk` e `biber` não estão instalados nesta máquina; PyYAML 6.0.3 está. Sem eles, o PDF não pode ser gerado aqui.
8. **DOI do código.** Pelo template, o DOI vem depois do aceite. A release v0.2.0 e o DOI do Zenodo ficam PENDENTES.
