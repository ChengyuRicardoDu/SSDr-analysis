# Simulations

Run the commands below from this directory.

## Prepare the reference

Simulation 2 resamples DLPFC expression. `reference.r` downloads section
151507 from the original authors and creates `reference.rds` locally:

```sh
Rscript reference.r
```

It takes the first 3,500 annotated spots in H5 order, samples 2,500 profiles
by layer with seed 2025, and selects 500 genes with `scran::modelGeneVar()`.
The resulting matrix contains raw counts. The downloads and generated RDS
files are not included in this upload.

This step needs Seurat, hdf5r, SingleCellExperiment, scuttle and scran.
The matrix and gene order were verified with Seurat 5.2.1,
SingleCellExperiment 1.28.1, scuttle 1.16.0 and scran 1.34.0.

Data: [Maynard et al., 2021](https://doi.org/10.1038/s41593-020-00787-0),
[official download list](https://github.com/LieberInstitute/HumanPilot/blob/master/AWS_File_locations.tsv).
The original reference was prepared through
[spatialLIBD](https://doi.org/10.1186/s12864-022-08601-w); this script recreates
the same input from the authors' H5 and annotation files.

## Generate data

```sh
Rscript run.r
```

The default generates the first evaluation replicate of each scenario in
`data/`. Set `reps <- 1:200` in `run.r` for all 200 replicates per scenario.
After generating the first replicate, use `reps <- 2:200` to continue.
Existing results are not overwritten. `seeds.csv` contains the original
generation, fitting and clustering seeds.

To generate a single dataset without saving:

```r
source("simulation1.r")
source("simulation2.r")
ref <- readRDS("reference.rds")
data1 <- simulation1(seed = 2026)
data2 <- simulation2(seed1 = 2031, seed2 = 2032, ref = ref)
```

## Fit SSDr

With `ssdr` 0.2.0 and the R torch backend installed:

```r
source("simulation_ssdr_f.r")
source("simulation_ssdr_nn.r")
seeds <- read.csv("seeds.csv")
s <- seeds[seeds$simulation == "simulation1" & seeds$rep == 1, ]
data <- readRDS("data/s1_evaluation_001.rds")
f <- simulation_ssdr_f(data, "simulation1", seed = s$kmeans)
nn <- simulation_ssdr_nn(data, "simulation1", seed1 = s$nn, seed2 = s$kmeans)
```

Both functions return `fit` and `cluster` without saving them. For Scenario 2,
use its seed row, dataset and `"simulation2"` argument. The simplified F
example does not cap scaled expression; the reported Scenario 2 analysis
used a cap of 10. The NN example uses raw counts.
