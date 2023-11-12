.DEFAULT_GOAL := data

RData/simulate_all.RData: simulate_all.R
	Rscript install_packages.R && Rscript simulate_all.R 2>&1 | tee simulate_all.Rout

data: RData/simulate_all.RData


figs: RData/simulate_all.RData
	if [ ! -f error.txt ]; then \
	  Rscript simulate_all_plot_results.R > simulate_all_plot_results.Rout 2>&1; \
	fi

clean:
	rm -rf error.* pars.rds

.PHONY: figs clean data
