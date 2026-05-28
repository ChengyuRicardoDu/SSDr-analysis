library(ssdr)

sim_path <- file.path("data", "simulations", "simulation1_sharp.rds")
if (!file.exists(sim_path)) {
  source(file.path("simulations", "simulation1_sharp.R"))
}

sim <- readRDS(sim_path)

set.seed(1)
spot_index <- sort(sample(seq_len(nrow(sim$counts)), 600))
gene_index <- seq_len(120)

counts <- sim$counts[spot_index, gene_index]
coords <- sim$coords[spot_index, c("x", "y")]

pca_input <- stats::prcomp(log1p(counts), rank. = 15)$x

fit_f <- ssdr_f(
  X = pca_input,
  coords = coords,
  rank = 2,
  bandwidth = 0.15,
  max_iter = 20
)

embedding_f <- fit_f$U
loadings_f <- fit_f$V
scale_f <- fit_f$sigma
parameters_f <- fit_f$parameters
diagnostics_f <- fit_f$diagnostics

run_ssdr_nn <- FALSE

if (run_ssdr_nn) {
  fit_nn <- ssdr_nn(
    X = counts,
    coords = coords,
    rank = 2,
    bandwidth = 0.15,
    lambda = 0.01,
    hidden_dim = 32,
    hidden_layers = 2,
    learning_rate = 0.005,
    max_iter = 50,
    patience = 10,
    seed = 1
  )

  embedding_nn <- fit_nn$U
  loadings_nn <- fit_nn$V
  scale_nn <- fit_nn$sigma
  parameters_nn <- fit_nn$parameters
  diagnostics_nn <- fit_nn$diagnostics
}
