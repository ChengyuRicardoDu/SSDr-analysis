out_dir <- file.path("data", "simulations")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

reference_path <- file.path("data", "reference", "dlpfc_bootstrap_reference.rds")
if (!file.exists(reference_path)) {
  stop("Missing reference pool: ", reference_path, call. = FALSE)
}

reference <- readRDS(reference_path)
X_pool <- reference$X
labels_pool <- as.character(reference$meta$label_factor)

geometry_seed <- 2031L
bootstrap_seed <- 2032L
grid_size <- 50L
domain_labels <- c("Layer1", "Layer3", "Layer5", "WM")
n_domains <- length(domain_labels)

grid <- expand.grid(row = seq_len(grid_size), col = seq_len(grid_size))
grid$spot_id <- sprintf("S%d", seq_len(nrow(grid)))
coords_mat <- as.matrix(grid[, c("col", "row")])

set.seed(geometry_seed)
base_seeds <- rbind(
  c(grid_size * 0.25, grid_size * 0.25),
  c(grid_size * 0.75, grid_size * 0.25),
  c(grid_size * 0.25, grid_size * 0.75),
  c(grid_size * 0.75, grid_size * 0.75)
)
jitter <- matrix(runif(n_domains * 2L, -grid_size * 0.05, grid_size * 0.05), n_domains, 2L)
seeds <- base_seeds + jitter

dist_sq <- sapply(seq_len(n_domains), function(i) {
  (coords_mat[, 1] - seeds[i, 1])^2 + (coords_mat[, 2] - seeds[i, 2])^2
})
domain_id <- max.col(-dist_sq, ties.method = "first")
domain_label <- domain_labels[domain_id]

set.seed(bootstrap_seed)
counts <- matrix(0L, nrow = nrow(coords_mat), ncol = ncol(X_pool))
colnames(counts) <- colnames(X_pool)
sample_index <- integer(nrow(coords_mat))
for (i in seq_len(nrow(coords_mat))) {
  pool_index <- which(labels_pool == domain_label[i])
  sample_index[i] <- sample(pool_index, 1L)
  counts[i, ] <- X_pool[sample_index[i], ]
}
storage.mode(counts) <- "integer"
rownames(counts) <- grid$spot_id

coords <- data.frame(
  spot_id = grid$spot_id,
  col = grid$col,
  row = grid$row,
  x = grid$col,
  y = grid$row
)

labels <- factor(domain_label, levels = domain_labels)
names(labels) <- grid$spot_id

spot_info <- data.frame(
  spot_id = grid$spot_id,
  row = grid$row,
  col = grid$col,
  domain_id = domain_id,
  domain_label = labels,
  reference_index = sample_index,
  reference_label = labels_pool[sample_index]
)

gene_info <- data.frame(
  gene_id = colnames(counts),
  source = "DLPFC_HVG",
  stringsAsFactors = FALSE
)

params <- list(
  scenario = "simulation2_bootstrap_voronoi",
  geometry_seed = geometry_seed,
  bootstrap_seed = bootstrap_seed,
  grid_size = grid_size,
  n_spots = nrow(counts),
  n_genes = ncol(counts),
  n_domains = n_domains,
  domain_labels = domain_labels,
  reference_pool_labels = levels(reference$meta$label_factor),
  reference_pool_spots = nrow(X_pool),
  reference_pool_genes = ncol(X_pool),
  seeds = seeds,
  geometry = "four Voronoi regions from jittered quarter-grid seeds",
  sampling = "one complete expression vector sampled with replacement from the matched reference label"
)

simulation <- list(
  description = "Simulation 2: four-domain DLPFC bootstrap Voronoi design",
  counts = counts,
  coords = coords,
  labels = labels,
  spot_info = spot_info,
  gene_info = gene_info,
  params = params
)

rds_path <- file.path(out_dir, "simulation2_bootstrap_voronoi.rds")
saveRDS(simulation, rds_path, compress = "xz", version = 2)

message("Wrote ", rds_path)
message("Counts: ", nrow(counts), " spots x ", ncol(counts), " genes")
message("Labels: ", paste(names(table(labels)), as.integer(table(labels)), collapse = ", "))
