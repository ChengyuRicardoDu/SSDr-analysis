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
- `parameters/`: SSDr-F and SSDr-NN parameter settings used for the manuscript
  analyses, including the final SSDr-NN tuning panels.
- `examples/basic_usage.R`: minimal example showing how to run SSDr-F, extract
  the embedding, compute a UMAP representation and plot it.

The simulation scripts write generated datasets to `data/simulations/`.
