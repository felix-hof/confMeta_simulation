.DEFAULT_GOAL := RData/simulate_all.RData

RData/simulate_all.RData: simulate_all.R ReplicationSuccess_extension.R
	Rscript install_packages.R && Rscript simulate_all.R > simulate_all.Rout 2>&1

figs: RData/simulate_all.RData
	Rscript simulate_all_plot_results.R > simulate_all_plot_results.Rout 2>&1

.PHONY: figs
