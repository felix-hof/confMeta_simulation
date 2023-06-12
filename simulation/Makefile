.DEFAULT_GOAL := RData/simulate_all.RData

RData/simulate_all.RData: simulate_all.R
	Rscript install_packages.R && Rscript simulate_all.R 2>&1 | tee simulate_all.Rout

figs: RData/simulate_all.RData
	if [ ! -f error.txt ]; then \
	  Rscript simulate_all_plot_results.R > simulate_all_plot_results.Rout 2>&1; \
	fi

.PHONY: figs
