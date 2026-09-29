# SSDr Analysis

This repository contains analysis materials supporting the SSDr manuscript.
The R package itself is maintained separately at
[ChengyuRicardoDu/SSDr](https://github.com/ChengyuRicardoDu/SSDr).

## Contents

- `simulations/`: simulation generators, seeds and SSDr fitting examples.
  Start with [the simulation instructions](simulations/README.md).
  `reference.r` downloads the source data and builds the Scenario 2 reference
  locally; `run.r` generates the evaluation datasets.
- `parameters/`: SSDr-F and SSDr-NN parameter settings used for the manuscript
  analyses.

Data are downloaded from the original authors by the scripts, not bundled
in this repository. The revised real-data examples will be added after review.
