.PHONY: all clean

all: protocol.pdf clean

protocol.pdf: protocol.tex biblio.bib 
	latexmk -bibtex -synctex=1 -pdf -pdflatex="pdflatex -interactive=nonstopmode" -use-make -halt-on-error $<
	$(MAKE) clean

clean:
	rm -f *~ *.log *.aux *.out Rplots.pdf comment.cut *.backup *.toc *.fls  *fdb_latexmk *~ *.synctex.gz *.brf *.blg *.bbl




