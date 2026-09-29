library(Seurat)
library(ssdr)
library(mclust)

realdata_ssdr_f <- function(counts, coords, genes, groups) {
  object <- CreateSeuratObject(counts = counts)
  object <- NormalizeData(object, normalization.method = "LogNormalize",
                          scale.factor = 10000, verbose = FALSE)
  object <- ScaleData(object, features = genes, verbose = FALSE)
  object <- RunPCA(object, features = genes, npcs = 15,
                   seed.use = 42, verbose = FALSE)
  X <- Embeddings(object, "pca")

  fit <- ssdr_f(
    X, coords, rank = 5, bandwidth = 0.15,
    lambda = NULL, center = FALSE, max_iter = 100, tol = 0.01
  )
  rownames(fit$U) <- rownames(X)
  U <- fit$U
  mclust::mclust.options(subset = nrow(U))
  set.seed(20260724)
  cluster <- Mclust(U, G = groups, modelNames = "EEE",
                    verbose = FALSE)$classification
  names(cluster) <- rownames(U)
  list(fit = fit, cluster = cluster)
}
