library(ssdr)
library(uwot)

source("simulations/simulation1_sharp.R")

sim <- readRDS("data/simulations/simulation1_sharp.rds")

set.seed(1)
spot_id <- sort(sample(seq_len(nrow(sim$counts)), 600))
gene_id <- seq_len(120)

counts <- sim$counts[spot_id, gene_id]
coords <- sim$coords[spot_id, c("x", "y")]

pca_input <- stats::prcomp(log1p(counts), rank. = 15)$x

fit <- ssdr_f(
  X = pca_input,
  coords = coords,
  rank = 2,
  bandwidth = 0.15,
  max_iter = 20
)

embedding <- fit$U
head(embedding)

umap_embedding <- uwot::umap(embedding, n_neighbors = 30, min_dist = 0.3)
colnames(umap_embedding) <- c("UMAP1", "UMAP2")
head(umap_embedding)

rescale01 <- function(x) {
  rng <- range(x)
  (x - rng[1]) / (rng[2] - rng[1])
}

umap_colour <- grDevices::hcl(
  h = 300 * rescale01(umap_embedding[, "UMAP1"]),
  c = 70,
  l = 35 + 45 * rescale01(umap_embedding[, "UMAP2"])
)

if (interactive()) {
  old_par <- par(mfrow = c(1, 2), mar = c(3, 3, 2, 1))

  plot(
    umap_embedding[, "UMAP1"],
    umap_embedding[, "UMAP2"],
    col = umap_colour,
    pch = 16,
    cex = 0.6,
    xlab = "UMAP1",
    ylab = "UMAP2",
    main = "UMAP"
  )

  plot(
    coords$x,
    -coords$y,
    col = umap_colour,
    pch = 16,
    cex = 0.6,
    xlab = "x",
    ylab = "y",
    main = "Spatial map",
    asp = 1
  )

  par(old_par)
}
