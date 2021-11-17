#!/bin/bash

cd simulation

make

cd ..

git add simulation/RData/* simulation/simulate_all.Rout

git commit -m "Push simulation results"

git push

