.DEFAULT_GOAL := data
README = README
data = RData/simulate_all.RData

data: simulation.R
	Rscript install_packages.R && Rscript simulation.R 2>&1 | tee simulate_all.Rout

all: $(README).pdf

$(README).pdf: $(README).Rnw
	Rscript -e "knitr::knit('$(README).Rnw')"
	latexmk -pdf $(README).tex

clean:
	latexmk -c

cleanall:
	latexmk -C

.PHONY: all data preview clean cleanall
