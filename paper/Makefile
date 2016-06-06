
MARKDOWN := title-page.md introduction.md results.md discussion.md methods.md references.md
TARGETS := title.yaml $(MARKDOWN)
NAME := paper
SUPP := supplement
BIBLIO := refs.bib
PANDOC_OPTS := --bibliography=${BIBLIO} --csl=plos.csl

all: html word pdf

html: $(NAME).html $(SUPP).html

word: $(NAME).docx $(SUPP).docx

pdf: $(NAME).pdf $(SUPP).pdf

## Build paper

$(NAME).html: $(TARGETS)
	pandoc $^ --to html --output $@ --standalone --section-divs --toc $(PANDOC_OPTS) --mathjax

$(NAME).docx: $(TARGETS)
	pandoc $^ --to docx --output $@ $(PANDOC_OPTS) --reference-docx=reference.docx

$(NAME).pdf: $(TARGETS)
	pandoc $^ --to latex --output $@ --toc $(PANDOC_OPTS)

# For debugging purposes, create intermediate tex file:
$(NAME).tex: $(TARGETS)
	pandoc $^ --to latex --output $@ --toc $(PANDOC_OPTS) --standalone

## Build supplement

$(SUPP).html: $(SUPP).md
	pandoc $^ --to html --output $@ --standalone --section-divs --toc $(PANDOC_OPTS) --mathjax

$(SUPP).docx: $(SUPP).md
	pandoc $^ --to docx --output $@ --self-contained --standalone --section-divs --toc $(PANDOC_OPTS)

$(SUPP).pdf: $(SUPP).md
	pandoc $^ --to latex --output $@ --self-contained --standalone --section-divs --toc $(PANDOC_OPTS)

## Knit R Markdown

%.md: %.Rmd
	R --slave -e "set.seed(100);knitr::knit('$<')"

## Format references

# Need to use Python 3
bibtex:
	python format-bibtex.py $(BIBLIO)

## Miscellaneous

.PHONY: clean bibtex

clean:
	rm $(MARKDOWN) $(NAME).* $(SUPP).md $(SUPP).html $(SUPP).docx $(SUPP).pdf

# /mnt/lustre/data/tools/rstudio_new/bin/pandoc/pandoc introduction.utf8.md --to html --from markdown+autolink_bare_uris+ascii_identifiers+tex_math_single_backslash-implicit_figures --output introduction.html --smart --email-obfuscation none --self-contained --standalone --section-divs --template /mnt/lustre/home/jdblischak/lib64/R/library/rmarkdown/rmd/h/default.html --variable 'theme:bootstrap' --include-in-header /tmp/RtmpY86uML/rmarkdown-str37db4af30af6.html --mathjax --variable 'mathjax-url:https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML' --no-highlight --variable highlightjs=/mnt/lustre/home/jdblischak/lib64/R/library/rmarkdown/rmd/h/highlight
