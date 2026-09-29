source("realdata_ssdr_f.r")
source("realdata_ssdr_nn.r")

samples <- c("151507", "151508", "151509", "151510", "151669", "151670",
             "151671", "151672", "151673", "151674", "151675", "151676")
method <- "F"
folder <- "fits"
dir.create(folder, showWarnings = FALSE)

url <- "https://raw.githubusercontent.com/LieberInstitute/HumanPilot/044446d6bd8fc154aa74f7be62ec67effb1ec376/10X/"
download.file(paste0(url, "barcode_level_layer_map.tsv"), "layers.tsv", mode = "wb")
layers <- read.delim("layers.tsv", header = FALSE,
                     col.names = c("barcode", "section", "layer"))
layers$layer <- sub("L", "Layer", layers$layer, fixed = TRUE)

for (section in samples) {
  file <- file.path(folder, paste0(section, "_", method, ".rds"))
  if (file.exists(file)) stop("File already exists: ", file)

  path <- file.path("..", "data", section)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  download.file(paste0(
    "https://spatial-dlpfc.s3.us-east-2.amazonaws.com/h5/",
    section, "_filtered_feature_bc_matrix.h5"
  ), file.path(path, "counts.h5"), mode = "wb")
  download.file(paste0(url, section, "/tissue_positions_list.txt"),
                file.path(path, "positions.csv"), mode = "wb")

  counts <- Read10X_h5(file.path(path, "counts.h5"), use.names = FALSE)
  p <- read.csv(file.path(path, "positions.csv"), header = FALSE,
                 col.names = c("barcode", "tissue", "row", "col", "y", "x"))
  p <- p[p$tissue == 1, ]
  counts <- counts[, p$barcode]
  coords <- data.frame(x = p$x, y = p$y, row.names = p$barcode)

  object <- CreateSeuratObject(counts = counts)
  object <- FindVariableFeatures(object, selection.method = "vst",
                                 nfeatures = 2000, verbose = FALSE)
  genes <- VariableFeatures(object)

  a <- layers[layers$section == section, ]
  annotation <- a$layer[match(p$barcode, a$barcode)]
  names(annotation) <- p$barcode
  groups <- if (section %in% c("151669", "151670", "151671", "151672")) 5 else 7

  if (method == "F") {
    result <- realdata_ssdr_f(counts, coords, genes, groups)
  } else {
    data <- list(counts = as.matrix(t(counts[genes[1:1000], ])),
                 coords = coords)
    result <- realdata_ssdr_nn(data, groups)
  }
  result$coords <- coords
  result$annotation <- annotation
  saveRDS(result, file)
}
