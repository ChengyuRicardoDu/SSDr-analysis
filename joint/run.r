source("joint_ssdr_f.r")

samples <- c("151507", "151508", "151509", "151510", "151669", "151670",
             "151671", "151672", "151673", "151674", "151675", "151676")
file <- "fits/joint.rds"
if (file.exists(file)) stop("File already exists: ", file)
dir.create("fits", showWarnings = FALSE)

layers <- read.delim("../realdata/layers.tsv", header = FALSE,
                     col.names = c("barcode", "section", "layer"))
layers$layer <- sub("L", "Layer", layers$layer, fixed = TRUE)

counts <- coords <- annotation <- scores <- list()
for (section in samples) {
  path <- file.path("..", "data", section)
  X <- Read10X_h5(file.path(path, "counts.h5"), use.names = FALSE)
  p <- read.csv(file.path(path, "positions.csv"), header = FALSE,
                 col.names = c("barcode", "tissue", "row", "col", "y", "x"))
  p <- p[p$tissue == 1, ]
  p <- p[match(colnames(X), p$barcode), ]
  counts[[section]] <- X
  coords[[section]] <- data.frame(x = p$x, y = p$y, row.names = p$barcode)

  object <- CreateSeuratObject(counts = X)
  object <- FindVariableFeatures(object, selection.method = "vst",
                                 nfeatures = 2000, verbose = FALSE)
  info <- HVFInfo(object, method = "vst", status = FALSE)
  scores[[section]] <- info[rownames(X), "variance.standardized"]

  a <- layers[layers$section == section, ]
  labels <- a$layer[match(colnames(X), a$barcode)]
  annotation[[section]] <- setNames(labels, colnames(X))
}

scores <- do.call(cbind, scores)
features <- rownames(counts[[1]])
percentiles <- apply(scores, 2, function(x)
  (rank(x, ties.method = "average") - 1) / (length(x) - 1))
i <- order(-rowMeans(percentiles), -rowMeans(scores), features, method = "radix")
genes <- features[i[1:2000]]

fit <- joint_ssdr_f(counts, coords, genes)
fit$coords <- coords
fit$annotation <- annotation
saveRDS(fit, file)
