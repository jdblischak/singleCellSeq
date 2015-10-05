
MARKDOWN := introduction.md methods.md results.md references.md
TARGETS := title.yaml $(MARKDOWN)
NAME := paper
PANDOC_OPTS := "--bibliography=refs.bib"

all: html word pdf

html: $(NAME).html

word: $(NAME).docx

pdf: $(NAME).pdf

$(NAME).html: $(TARGETS)
	pandoc $^ --to html --output $(NAME).html --self-contained --standalone --section-divs --toc $(PANDOC_OPTS)

$(NAME).docx: $(TARGETS)
	pandoc $^ --to docx --output $(NAME).docx $(PANDOC_OPTS)

$(NAME).pdf: $(TARGETS)
	pandoc $^ --to latex --output $(NAME).pdf --toc $(PANDOC_OPTS)

%.md: %.Rmd
	R --slave -e "set.seed(100);knitr::knit('$<')"

.PHONY: clean
clean:
	rm $(MARKDOWN) $(NAME).*

# /mnt/lustre/data/tools/rstudio_new/bin/pandoc/pandoc introduction.utf8.md --to html --from markdown+autolink_bare_uris+ascii_identifiers+tex_math_single_backslash-implicit_figures --output introduction.html --smart --email-obfuscation none --self-contained --standalone --section-divs --template /mnt/lustre/home/jdblischak/lib64/R/library/rmarkdown/rmd/h/default.html --variable 'theme:bootstrap' --include-in-header /tmp/RtmpY86uML/rmarkdown-str37db4af30af6.html --mathjax --variable 'mathjax-url:https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML' --no-highlight --variable highlightjs=/mnt/lustre/home/jdblischak/lib64/R/library/rmarkdown/rmd/h/highlight
