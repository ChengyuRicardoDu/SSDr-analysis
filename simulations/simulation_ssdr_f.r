simulation_ssdr_f <- function(data, simulation, seed) {
  X <- data$counts
  total <- rowSums(X)
  total[total <= 0] <- 1

  if (simulation == "simulation1") {
    X <- log1p(X / total * median(total))
    bandwidth <- 0.15
    k <- 7
  } else {
    X <- log1p(X / total * 10000)
    bandwidth <- 0.10
    k <- 4
  }

  X <- scale(X)
  X[!is.finite(X)] <- 0

  X <- prcomp(X, center = FALSE, scale. = FALSE, rank. = 15)$x[, 1:15]

  fit <- ssdr::ssdr_f(
    X, data$coords, rank = 7, bandwidth = bandwidth,
    lambda = NULL, center = FALSE, max_iter = 100, tol = 0.01
  )
  rownames(fit$U) <- rownames(X)

  set.seed(seed)
  cluster <- kmeans(fit$U, centers = k, nstart = 20,
                    iter.max = 100, algorithm = "Hartigan-Wong")$cluster

  list(fit = fit, cluster = cluster)
}
