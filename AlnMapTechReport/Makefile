MAIN = main

.PHONY: all clean

all: $(MAIN).pdf

$(MAIN).pdf: $(MAIN).tex
	pdflatex $(MAIN).tex
	pdflatex $(MAIN).tex   # second pass for TOC / cross-references

clean:
	rm -f *.aux *.log *.out *.toc *.pdf
