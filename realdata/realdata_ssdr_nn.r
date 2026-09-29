realdata_ssdr_nn <- function(data, groups, seeds = 1:3) {
  seeds <- sort(seeds)
  suppressPackageStartupMessages(library(mclust))
  torch::torch_set_num_threads(1)
  fits <- lapply(seeds, function(seed) {
    ssdr::ssdr_nn(
      data$counts, data$coords, rank = 5, bandwidth = 0.05, lambda = 0.01,
      hidden_dim = 64, hidden_layers = 2, learning_rate = 0.01,
      max_iter = 2000, patience = 50, tol = 1e-5, seed = seed
    )
  })

  objectives <- sapply(fits, function(fit) fit$objective)
  best <- which.min(objectives)
  fit <- fits[[best]]

  rownames(fit$U) <- rownames(data$counts)
  U <- fit$U
  mclust::mclust.options(subset = nrow(U))
  set.seed(20260724)
  cluster <- mclust::Mclust(U, G = groups, modelNames = "EEE",
                            verbose = FALSE)$classification
  names(cluster) <- rownames(U)
  list(fit = fit, cluster = cluster,
       objectives = data.frame(seed = seeds, objective = objectives))
}
