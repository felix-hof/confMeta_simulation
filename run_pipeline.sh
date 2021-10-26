#!/bin/bash

cd simulation

make

cd ..

git add simulation/RData/* simulation/figs/* simulation/simulate_all.Rout simulation/simulate_all_plot_results.Rout

git commit -m "Push simulation results"

git push

