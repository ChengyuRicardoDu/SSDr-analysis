out_dir <- file.path("data", "simulations")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

seed <- 2026L
grid_size <- 50L
n_layers <- 7L
markers_per_layer <- 15L
n_housekeeping <- 195L
mu_background <- 1.0
mu_marker <- 2.0
nb_size <- 0.5

layer_labels <- paste0("L", seq_len(n_layers))

grid <- expand.grid(row = seq_len(grid_size), col = seq_len(grid_size))
grid$spot_id <- sprintf("S%d", seq_len(nrow(grid)))
grid$domain_label <- cut(
  grid$col,
  breaks = n_layers,
  labels = layer_labels,
  include.lowest = TRUE
)
grid$domain_id <- as.integer(grid$domain_label)

set.seed(seed)

count_columns <- vector("list", n_layers * markers_per_layer + n_housekeeping)
gene_info <- data.frame(
  gene_id = character(),
  gene_type = character(),
  marker_layer_id = integer(),
  marker_layer_label = character(),
  stringsAsFactors = FALSE
)

idx <- 0L
for (layer_id in seq_len(n_layers)) {
  layer_label <- layer_labels[[layer_id]]
  for (marker_id in seq_len(markers_per_layer)) {
    idx <- idx + 1L
    gene_id <- sprintf("M_L%d_%d", layer_id, marker_id)
    mu <- ifelse(grid$domain_label == layer_label, mu_marker, mu_background)
    count_columns[[idx]] <- rnbinom(nrow(grid), mu = mu, size = nb_size)
    gene_info[idx, ] <- list(gene_id, "marker", layer_id, layer_label)
  }
}

for (gene_id_num in seq_len(n_housekeeping)) {
  idx <- idx + 1L
  gene_id <- sprintf("HK_%d", gene_id_num)
  count_columns[[idx]] <- rnbinom(nrow(grid), mu = mu_background, size = nb_size)
  gene_info[idx, ] <- list(gene_id, "non_marker", NA_integer_, NA_character_)
}

counts <- do.call(cbind, count_columns)
storage.mode(counts) <- "integer"
colnames(counts) <- gene_info$gene_id
rownames(counts) <- grid$spot_id

coords <- data.frame(
  spot_id = grid$spot_id,
  col = grid$col,
  row = grid$row,
  x = grid$col,
  y = grid$row
)

labels <- factor(grid$domain_label, levels = layer_labels)
names(labels) <- grid$spot_id

spot_info <- data.frame(
  spot_id = grid$spot_id,
  row = grid$row,
  col = grid$col,
  domain_id = grid$domain_id,
  domain_label = labels
)

params <- list(
  scenario = "simulation1_sharp",
  seed = seed,
  grid_size = grid_size,
  n_spots = nrow(counts),
  n_genes = ncol(counts),
  n_layers = n_layers,
  markers_per_layer = markers_per_layer,
  n_marker_genes = n_layers * markers_per_layer,
  n_housekeeping = n_housekeeping,
  mu_background = mu_background,
  mu_marker = mu_marker,
  fold_change = mu_marker / mu_background,
  nb_size = nb_size,
  nb_parameterization = "R rnbinom(mu, size); Var(Y) = mu + mu^2 / size",
  layer_rule = "cut(col, breaks = 7, include.lowest = TRUE)"
)

simulation <- list(
  description = "Simulation 1: sharp laminar Negative Binomial design",
  counts = counts,
  coords = coords,
  labels = labels,
  spot_info = spot_info,
  gene_info = gene_info,
  params = params
)

rds_path <- file.path(out_dir, "simulation1_sharp.rds")
saveRDS(simulation, rds_path, compress = "xz", version = 2)

message("Wrote ", rds_path)
message("Counts: ", nrow(counts), " spots x ", ncol(counts), " genes")
message("Labels: ", paste(names(table(labels)), as.integer(table(labels)), collapse = ", "))
