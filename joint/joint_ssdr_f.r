library(Seurat)
library(Matrix)
library(ssdr)

joint_ssdr_f <- function(counts, coords, genes) {
  G <- lapply(counts, function(X) {
    object <- CreateSeuratObject(counts = X)
    object <- NormalizeData(object, normalization.method = "LogNormalize",
                            scale.factor = 10000, verbose = FALSE)
    normalized <- LayerData(object[["RNA"]], layer = "data")
    t(normalized[genes, , drop = FALSE])
  })

  C <- matrix(0, length(genes), length(genes), dimnames = list(genes, genes))
  for (X in G) {
    C <- C + cov(as.matrix(X)) / length(G)
  }
  d <- sqrt(diag(C))
  C <- cov2cor(C)
  L <- eigen(C, symmetric = TRUE)$vectors[, 1:15]
  dimnames(L) <- list(genes, paste0("PC", 1:15))
  for (j in 1:15) {
    if (L[which.max(abs(L[, j])), j] < 0) L[, j] <- -L[, j]
  }

  # Keep section means in the common PC scores.
  loadings <- L / d
  Y <- lapply(G, function(X) as.matrix(X %*% loadings))
  fit <- ssdr_joint_f(
    Y, coords, rank = 5, bandwidth = 0.10, lambda = 1,
    solver = "nystrom", n_landmarks = 200, landmark_seed = 20260711
  )
  fit$pca <- list(genes = genes, loadings = L, scale = d)
  fit
}
