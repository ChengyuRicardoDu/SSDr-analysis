simulation_ssdr_nn <- function(data, simulation, seed1, seed2) {
  if (simulation == "simulation1") {
    bandwidth <- 0.10
    lambda <- 0.005
    hidden_dim <- 128
    max_iter <- 500
    k <- 7
  } else {
    bandwidth <- 0.20
    lambda <- 0.03
    hidden_dim <- 64
    max_iter <- 1200
    k <- 4
  }

  torch::torch_set_num_threads(1)
  fit <- ssdr::ssdr_nn(
    data$counts, data$coords,
    rank = 5, bandwidth = bandwidth, lambda = lambda,
    hidden_dim = hidden_dim, hidden_layers = 2,
    learning_rate = 0.005, max_iter = max_iter,
    patience = 50, tol = 1e-5, seed = seed1
  )
  rownames(fit$U) <- rownames(data$counts)

  set.seed(seed2)
  cluster <- kmeans(fit$U, centers = k, nstart = 20,
                    iter.max = 100, algorithm = "Hartigan-Wong")$cluster

  list(fit = fit, cluster = cluster)
}
