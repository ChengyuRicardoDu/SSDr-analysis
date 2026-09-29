library(SingleCellExperiment)

options(timeout = 300)
dir.create("data", showWarnings = FALSE)
if (!file.exists("data/151507.h5")) {
  download.file(
    "https://spatial-dlpfc.s3.us-east-2.amazonaws.com/h5/151507_filtered_feature_bc_matrix.h5",
    "data/151507.h5", mode = "wb"
  )
}
if (!file.exists("data/layers.tsv")) {
  download.file(paste0(
    "https://raw.githubusercontent.com/LieberInstitute/HumanPilot/",
    "044446d6bd8fc154aa74f7be62ec67effb1ec376/10X/barcode_level_layer_map.tsv"
  ), "data/layers.tsv", mode = "wb")
}

X <- Seurat::Read10X_h5("data/151507.h5", use.names = FALSE)
a <- read.delim("data/layers.tsv", header = FALSE,
                col.names = c("barcode", "section", "layer"))
a <- a[a$section == 151507, ]
pool <- a$layer[match(colnames(X), a$barcode)]
pool <- sub("^L([1-6])$", "Layer\\1", pool)
i <- which(!is.na(pool))[1:3500]
X <- X[, i]
pool <- pool[i]

grid <- expand.grid(row = 1:50, col = 1:50)
layers <- c(paste0("Layer", 1:6), "WM")
labels <- cut(grid$col, breaks = seq(0, 50, length.out = 8),
              labels = layers, include.lowest = TRUE)
set.seed(2025)
index <- integer(2500)
for (i in 1:2500) index[i] <- sample(which(pool == labels[i]), 1)
X <- X[, index]
colnames(X) <- paste0("Spot_", 1:2500)

sce <- SingleCellExperiment(assays = list(counts = X))
sce <- scuttle::logNormCounts(sce)
genes <- scran::getTopHVGs(scran::modelGeneVar(sce), n = 500)
ref <- list(
  X = t(as.matrix(X[genes, ])),
  coords = as.matrix(grid[, c("col", "row")]),
  meta = data.frame(x = grid$col, y = grid$row, label_factor = labels),
  k = nlevels(labels), top_hvgs = genes
)
saveRDS(ref, "reference.rds")
