library(ssdr)
library(uwot)

if (!file.exists("data/simulations/simulation1_sharp.rds")) {
  source("simulations/simulation1_sharp.R")
}

sim <- readRDS("data/simulations/simulation1_sharp.rds")

set.seed(1)
spot_id <- sort(sample(seq_len(nrow(sim$counts)), 600))
gene_id <- seq_len(120)

counts <- sim$counts[spot_id, gene_id]
coords <- sim$coords[spot_id, c("x", "y")]
labels <- sim$labels[spot_id]

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

dir.create("outputs", showWarnings = FALSE)

pdf("outputs/basic_usage_umap.pdf", width = 4, height = 3.5)
plot(
  umap_embedding,
  col = as.integer(labels),
  pch = 16,
  cex = 0.6,
  xlab = "UMAP1",
  ylab = "UMAP2",
  main = "SSDr-F embedding"
)
legend(
  "topright",
  legend = levels(labels),
  col = seq_along(levels(labels)),
  pch = 16,
  cex = 0.7,
  bty = "n"
)
dev.off()
