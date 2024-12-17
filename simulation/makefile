.DEFAULT_GOAL := data
README = README
data = RData/simulate_all.RData

data: simulation.R
	Rscript install_packages.R && Rscript simulation.R

figs: data plot_results.R
	Rscript plot_results.R 2>&1 | tee plot_results.Rout

$(README).pdf: $(README).Rnw
	Rscript -e "knitr::knit('$(README).Rnw')"
	latexmk -pdf $(README).tex

readme: $(README).pdf

all: data figs readme

clean:
	latexmk -c

cleanall:
	latexmk -C

.PHONY: data figs readme all clean cleanall
