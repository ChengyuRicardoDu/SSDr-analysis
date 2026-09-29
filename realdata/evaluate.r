library(mclust)
library(FNN)
library(Matrix)

files <- list.files("fits", pattern = "[.]rds$", full.names = TRUE)
if (!length(files)) stop("No fit files found in fits/.")
dir.create("evaluation", showWarnings = FALSE)
metrics <- list()
components <- list()

for (file in files) {
  name <- tools::file_path_sans_ext(basename(file))
  id <- strsplit(name, "_", fixed = TRUE)[[1]]
  groups <- if (id[1] %in% as.character(151669:151672)) 5 else 7
  data <- readRDS(file)
  xy <- as.matrix(data$coords)
  barcodes <- rownames(xy)
  cluster <- data$cluster[barcodes]
  stopifnot(!anyNA(cluster))
  n <- nrow(xy)

  # ARI
  truth <- data$annotation[barcodes]
  annotated <- !is.na(truth)
  ari <- adjustedRandIndex(cluster[annotated], truth[annotated])

  # CHAOS
  scaled <- scale(xy)
  distance <- numeric(n)
  for (i in split(1:n, cluster)) {
    if (length(i) > 1) distance[i] <- get.knn(scaled[i, , drop = FALSE], k = 1)$nn.dist[, 1]
  }
  chaos <- mean(distance)

  # PAS
  neighbours <- get.knn(xy, k = 10)$nn.index
  neighbour_labels <- matrix(cluster[neighbours], ncol = 10)
  pas <- 100 * mean(rowSums(neighbour_labels != cluster) >= 6)

  # Moran's I
  U <- data$fit$U[barcodes, 1:5, drop = FALSE]
  nn <- get.knn(xy, k = 6)$nn.index
  W <- sparseMatrix(i = rep(1:n, each = 6), j = as.vector(t(nn)), x = 1,
                    dims = c(n, n))
  W <- drop0(((W + t(W)) > 0) * 1)
  W <- Diagonal(x = 1 / rowSums(W)) %*% W
  z <- sweep(U, 2, colMeans(U), "-")
  moran <- n / sum(W) * colSums(z * (W %*% z)) / colSums(z^2)

  metrics[[name]] <- data.frame(
    sample_id = id[1], method = paste0("SSDr-", id[2]), n_spots = n,
    n_annotated = sum(annotated), G = groups,
    ARI = ari, CHAOS = chaos, PAS_pct = pas, Moran_first5 = mean(moran)
  )
  components[[name]] <- data.frame(
    sample_id = id[1], method = paste0("SSDr-", id[2]), component_index = 1:5,
    component_variance = apply(U, 2, var), moran_i = moran
  )
  message(name)
}
write.csv(do.call(rbind, metrics), "evaluation/metrics.csv", row.names = FALSE)
write.csv(do.call(rbind, components), "evaluation/moran.csv", row.names = FALSE)
