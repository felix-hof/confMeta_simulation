.DEFAULT_GOAL := figs

RData/simulate_all.RData: simulate_all.R ReplicationSuccess_extension_LH.R
	Rscript install_packages.R && Rscript simulate_all.R > simulate_all.Rout 2>&1

figs: RData/simulate_all.RData
	if [ ! -f error.txt ]; then \
	  Rscript simulate_all_plot_results.R > simulate_all_plot_results.Rout 2>&1; \
	fi

.PHONY: figs
