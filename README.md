# SSDr Analysis

This repository contains analysis materials supporting the SSDr manuscript.
The R package itself is maintained separately at
[ChengyuRicardoDu/SSDr](https://github.com/ChengyuRicardoDu/SSDr).

## Contents

- `simulations/`: simulation generators, seeds and SSDr fitting examples.
  Start with [the simulation instructions](simulations/README.md).
  `reference.r` downloads the source data and builds the Scenario 2 reference
  locally; `run.r` generates the evaluation datasets.
- `realdata/`: DLPFC SSDr-F and SSDr-NN fitting, evaluation and S1-S3 plotting.
  Follow [the real-data instructions](realdata/README.md).

Data are downloaded from the original authors by the scripts, not bundled
in this repository.

Use `ssdr` 0.2.0 and `mclust` 6.1.1, with Seurat, Matrix, hdf5r, FNN and
gtools. SSDr-NN also needs R torch and its installed backend.
