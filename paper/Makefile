# Builds the PDF with the official ReScience C template, pinned to one commit.
# The template files (GPL-3+; fonts under OFL) are fetched into this directory
# and are not committed: see .gitignore and README.md.
TEMPLATE_URL    = https://github.com/ReScience/template.git
TEMPLATE_COMMIT = 0978d0e5d9a86246789fd3468aec6fc71e39c944
TEMPLATE_FILES  = rescience.cls header.tex yaml-to-latex.py article.py roboto source-serif-pro source-sans-pro source-code-pro

.PHONY: all template clean
all: article.pdf

.template:
	git clone -q $(TEMPLATE_URL) .template
	git -C .template checkout -q $(TEMPLATE_COMMIT)

template: .template
	cp -r $(addprefix .template/,$(TEMPLATE_FILES)) .

metadata.tex: metadata.yaml | template
	LC_ALL=en_US.UTF-8 python3 ./yaml-to-latex.py -i $< -o $@

article.pdf: article.tex content.tex bibliography.bib metadata.tex generated/numbers.tex | template
	latexmk -pdf -pdflatex="xelatex -interaction=nonstopmode -halt-on-error" article.tex

clean:
	-latexmk -CA
	rm -rf .template metadata.tex *.bbl *.run.xml $(TEMPLATE_FILES)
