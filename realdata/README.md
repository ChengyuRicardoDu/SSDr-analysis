# DLPFC analysis

Run from this directory:

1. Choose `samples` in `run.r`, set `method <- "F"`, and run `Rscript run.r`.
2. Set `method <- "NN"` and run `Rscript run.r` again.
3. Run `Rscript evaluate.r` to calculate the metrics.
4. Run `Rscript figures.r` to draw S1-S3.


For the paper comparison, section 151510 excluded `AGAAGAGCGCCGTTCC-1`
before GMM, but retained it during fitting. 
