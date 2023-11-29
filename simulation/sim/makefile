TARGET = README

all: $(TARGET).pdf

$(TARGET).pdf: $(TARGET).Rnw
	Rscript -e "knitr::knit('$(TARGET).Rnw')"
	latexmk -pdf $(TARGET).tex

clean:
	latexmk -c

cleanall:
	latexmk -C

.PHONY: all preview clean cleanall
