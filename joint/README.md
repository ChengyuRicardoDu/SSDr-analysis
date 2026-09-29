# Joint DLPFC analysis

Reuse the 12 sections downloaded by `realdata/run.r`: counts and coordinates
in `../data/`, and annotations in `../realdata/layers.tsv`.

Run from this directory:

1. Run `Rscript run.r` to select common genes and fit Joint SSDr-F.
2. Run `Rscript figures.r` to draw S4-S6.

The fit is saved as `fits/joint.rds`; the figures are `figs4.pdf`,
`figs5.pdf` and `figs6.pdf`.

This uses `ssdr` 0.2.0. The tested run selects one different gene from the
paper's 2,000-gene list. Landmark sampling and spectral truncation also
differ from the paper implementation, so the numerical results are not identical.
