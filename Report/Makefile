TEX = $(wildcard *.tex)
BIB = $(wildcard *.bib)
EPS = $(wildcard *.eps)
FIG = $(patsubst %.eps, %-eps-converted-to.pdf, $(EPS))
PAPER = paper
VPATH = figure

$(PAPER): $(PAPER).pdf
	open $<

$(PAPER).pdf: $(TEX) $(BIB) $(FIG)
	@echo $(FIG)
	pdflatex $(basename $@)
	bibtex $(basename $@)
	pdflatex $(basename $@)
	pdflatex $(basename $@)

$(FIG): %-eps-converted-to.pdf: %.eps:
	epstopdf $<

clean:
	rm -f *.pdf *.aux *.bbl *.blg *.log *.out figure/*-eps-converted-to.pdf
