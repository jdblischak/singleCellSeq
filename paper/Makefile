
MARKDOWN := introduction.md results.md methods.md references.md
TARGETS := title.yaml $(MARKDOWN)
NAME := paper
SUPP := supplement
BIBLIO := refs.bib
PANDOC_OPTS := --bibliography=${BIBLIO}

all: html word pdf

html: $(NAME).html $(SUPP).html

word: $(NAME).docx $(SUPP).docx

pdf: $(NAME).pdf $(SUPP).pdf

## Build paper

$(NAME).html: $(TARGETS)
	pandoc $^ --to html --output $@ --self-contained --standalone --section-divs --toc $(PANDOC_OPTS)

$(NAME).docx: $(TARGETS)
	pandoc $^ --to docx --output $@ $(PANDOC_OPTS)

$(NAME).pdf: $(TARGETS)
	pandoc $^ --to latex --output $@ --toc $(PANDOC_OPTS)

## Build supplement

$(SUPP).html: $(SUPP).md
	pandoc $^ --to html --output $@ --self-contained --standalone --section-divs --toc $(PANDOC_OPTS)

$(SUPP).docx: $(SUPP).md
	pandoc $^ --to docx --output $@ --self-contained --standalone --section-divs --toc $(PANDOC_OPTS)

$(SUPP).pdf: $(SUPP).md
	pandoc $^ --to latex --output $@ --self-contained --standalone --section-divs --toc $(PANDOC_OPTS)

## Knit R Markdown

%.md: %.Rmd
	R --slave -e "set.seed(100);knitr::knit('$<')"

## Format references

bibtex:
	python format-bibtex.py $(BIBLIO)

## Miscellaneous

.PHONY: clean bibtex

clean:
	rm $(MARKDOWN) $(NAME).*

# /mnt/lustre/data/tools/rstudio_new/bin/pandoc/pandoc introduction.utf8.md --to html --from markdown+autolink_bare_uris+ascii_identifiers+tex_math_single_backslash-implicit_figures --output introduction.html --smart --email-obfuscation none --self-contained --standalone --section-divs --template /mnt/lustre/home/jdblischak/lib64/R/library/rmarkdown/rmd/h/default.html --variable 'theme:bootstrap' --include-in-header /tmp/RtmpY86uML/rmarkdown-str37db4af30af6.html --mathjax --variable 'mathjax-url:https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML' --no-highlight --variable highlightjs=/mnt/lustre/home/jdblischak/lib64/R/library/rmarkdown/rmd/h/highlight
