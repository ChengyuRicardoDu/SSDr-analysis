# SSDr Analysis

This repository contains analysis materials supporting the SSDr manuscript.
The R package itself is maintained separately at
[ChengyuRicardoDu/SSDr](https://github.com/ChengyuRicardoDu/SSDr).

## Contents

- `simulations/`: simulation generators, seeds and SSDr fitting examples.
  Start with [the simulation instructions](simulations/README.md).
  `reference.r` downloads the source data and builds the Scenario 2 reference
  locally; `run.r` generates the evaluation datasets.
- `data/dlpfc_151673/`: raw filtered 10X h5 count matrix and matching
  `tissue_positions_list.txt` for DLPFC sample `151673`.
- `parameters/`: SSDr-F and SSDr-NN parameter settings used for the manuscript
  analyses.
- `examples/dlpfc_real_data.Rmd`: runnable DLPFC walkthrough starting from raw
  counts and positions, then fitting SSDr-F and SSDr-NN with spatial plots.
