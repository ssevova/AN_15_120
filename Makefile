LATEX = pdflatex
LATEXFLAGS =

DOCFILE = note.pdf

RERUN = '(No file .*\.toc|There were undefined references)'

.PHONY: all
all: $(DOCFILE)

%.pdf: %.tex $(wildcard *.tex) 
	@echo "====> LaTeX first pass"
	$(LATEX) $(LATEXFLAGS) $(<:.tex=)
	@if egrep -q $(RERUN) $*.log ; then echo "====> LaTeX rerun" && $(LATEX) $<; fi
	@if egrep -q $(RERUN) $*.log ; then echo "====> LaTeX rerun" && $(LATEX) $<; fi

.PHONY: clean
clean:
	rm -f *.pdf *.aux *.dvi *.log *.qsl *.sol *.lof *.lot *.toc *~

