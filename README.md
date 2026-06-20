# SSDr Analysis

This repository contains analysis materials supporting the SSDr manuscript.
The R package itself is maintained separately at
[ChengyuRicardoDu/SSDr](https://github.com/ChengyuRicardoDu/SSDr).

## Contents

- `simulations/simulation1_sharp.R`: generates the sharp laminar Negative
  Binomial simulation.
- `simulations/simulation2_bootstrap_voronoi.R`: generates the DLPFC bootstrap
  Voronoi simulation.
- `data/reference/dlpfc_bootstrap_reference.rds`: reference pool used by
  `simulation2_bootstrap_voronoi.R`.
- `data/dlpfc_151673/`: raw filtered 10X h5 count matrix and matching
  `tissue_positions_list.txt` for DLPFC sample `151673`.
- `parameters/`: SSDr-F and SSDr-NN parameter settings used for the manuscript
  analyses, including the final SSDr-NN tuning panels.
- `examples/dlpfc_real_data.Rmd`: runnable DLPFC walkthrough starting from raw
  counts and positions, then fitting SSDr-F and SSDr-NN with spatial plots.
- `examples/dlpfc_real_data.html`: rendered version of the DLPFC walkthrough
  with the analysis code and output figures.

The simulation scripts write generated datasets to `data/simulations/`.
