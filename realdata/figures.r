library(gtools)

folder <- "fits"
sections <- list(151507:151510, 151669:151672, 151673:151676)
layers <- c(paste("Layer", 1:6), "WM")
cols <- c(
  "0" = "#BDBDBD", "1" = "#332288", "2" = "#67B7D1", "3" = "#369985",
  "4" = "#117733", "5" = "#89892C", "6" = "#CBB85D", "7" = "#C55A70"
)

draw <- function(section, labels) {
  f <- readRDS(file.path(folder, paste0(section, "_F.rds")))
  nn <- readRDS(file.path(folder, paste0(section, "_NN.rds")))
  barcodes <- rownames(f$coords)
  reference <- match(f$annotation[barcodes], c(paste0("Layer", 1:6), "WM"))
  reference[is.na(reference)] <- 0
  values <- cbind(reference, f$cluster[barcodes], nn$cluster[barcodes])
  G <- length(labels)
  orders <- permutations(G, G)
  for (j in 2:3) {
    counts <- table(factor(values[, j], levels = 1:G),
                    factor(reference, levels = labels))
    scores <- apply(orders, 1, function(order) sum(counts[cbind(1:G, order)]))
    values[, j] <- labels[orders[which.max(scores), ]][values[, j]]
  }

  plot.new()
  text(0.5, 0.5, section, cex = 0.9)
  xy <- cbind(f$coords$x, -f$coords$y)
  xlim <- range(xy[, 1])
  ylim <- range(xy[, 2])
  xlim <- xlim + c(-1, 1) * diff(xlim) * 0.01
  ylim <- ylim + c(-1, 1) * diff(ylim) * 0.01
  for (j in 1:3) {
    plot(xy, pch = 16, cex = 0.35, col = cols[as.character(values[, j])],
         asp = 1, xlim = xlim, ylim = ylim, axes = FALSE,
         xlab = "", ylab = "", xaxs = "i", yaxs = "i")
  }
}

for (i in seq_along(sections)) {
  labels <- if (i == 2) 3:7 else 1:7
  pdf(paste0("figs", i, ".pdf"), width = 166 / 25.4, height = 208 / 25.4,
      pointsize = 8.5, family = "Helvetica", useDingbats = FALSE)
  layout(rbind(c(0, 1:3), 4:7, 8:11, 12:15, 16:19, rep(20, 4)),
         widths = c(16, 50, 50, 50), heights = c(8, 48, 48, 48, 48, 8))
  par(mar = rep(0.05, 4), cex = 1)
  for (title in c("Reference", "SSDr-F", "SSDr-NN")) {
    plot.new()
    text(0.5, 0.5, title, font = 2)
  }

  for (section in as.character(sections[[i]])) draw(section, labels)

  plot.new()
  legend("center", legend = c(layers[labels], "Unlabelled"),
         col = cols[c(as.character(labels), "0")], pch = 16,
         bty = "n", horiz = TRUE, cex = 0.85, x.intersp = 0.5)
  dev.off()
}
