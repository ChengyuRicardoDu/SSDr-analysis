fit <- readRDS("fits/joint.rds")
sections <- list(151507:151510, 151669:151672, 151673:151676)
cols <- c(
  Layer1 = "#332288", Layer2 = "#67B7D1", Layer3 = "#369985",
  Layer4 = "#117733", Layer5 = "#89892C", Layer6 = "#CBB85D",
  WM = "#C55A70"
)
pal <- colorRampPalette(
  c("#2C7FB8", "#F7F7F7", "#D98E04"), space = "Lab"
)(513)
limits <- apply(abs(do.call(rbind, fit$U)), 2, quantile,
                probs = 0.99, type = 8)

draw <- function(section) {
  plot.new()
  text(0.5, 0.5, section, srt = 90, cex = 0.9)
  xy <- cbind(fit$coords[[section]][, 1], -fit$coords[[section]][, 2])
  xlim <- range(xy[, 1])
  ylim <- range(xy[, 2])
  xlim <- xlim + c(-1, 1) * diff(xlim) * 0.01
  ylim <- ylim + c(-1, 1) * diff(ylim) * 0.01
  labels <- fit$annotation[[section]]
  col <- cols[labels]
  col[is.na(labels)] <- "#BDBDBD"
  plot(xy, pch = 16, cex = 0.27, col = adjustcolor(col, 0.98),
       asp = 1, xlim = xlim, ylim = ylim, axes = FALSE,
       xlab = "", ylab = "", xaxs = "i", yaxs = "i")

  dx <- diff(par("usr")[1:2]) / par("pin")[1] * 0.12 / 25.4
  dy <- diff(par("usr")[3:4]) / par("pin")[2] * 0.12 / 25.4
  d <- expand.grid(y = c(-dy, 0, dy), x = c(-dx, 0, dx))
  for (layer in names(cols)) {
    i <- which(labels == layer)
    if (length(i) == 0) next
    j <- i[which.min((xy[i, 1] - median(xy[i, 1]))^2 +
                     (xy[i, 2] - median(xy[i, 2]))^2)]
    x <- xy[j, 1]
    y <- xy[j, 2]
    if (section == "151510" && layer == "WM") x <- x - 0.06 * diff(xlim)
    lab <- sub("Layer", "L", layer)
    text(x + d$x, y + d$y, lab, cex = 0.75, font = 2, col = "white")
    text(x, y, lab, cex = 0.75, font = 2, col = "#222222")
  }

  for (k in 1:5) {
    z <- fit$U[[section]][, k] / limits[k]
    z <- pmax(-1, pmin(1, z))
    col <- pal[1 + round((z + 1) * 256)]
    plot(xy, pch = 16, cex = 0.27, col = adjustcolor(col, 0.98),
         asp = 1, xlim = xlim, ylim = ylim, axes = FALSE,
         xlab = "", ylab = "", xaxs = "i", yaxs = "i")
  }
}

for (page in seq_along(sections)) {
  pdf(paste0("figs", page + 3, ".pdf"), width = 180 / 25.4, height = 144 / 25.4,
      pointsize = 8.5, family = "Helvetica", useDingbats = FALSE)
  layout(rbind(c(0, 1:6), 7:13, 14:20, 21:27, 28:34, rep(35, 7)),
         widths = c(12, rep(28, 6)), heights = c(8, 32, 32, 32, 32, 8))
  par(mar = rep(0.05, 4), cex = 1)
  for (title in c("Reference", paste("Component", 1:5))) {
    plot.new()
    text(0.5, 0.5, title, font = 2)
  }
  for (section in as.character(sections[[page]])) draw(section)

  plot.new()
  plot.window(xlim = c(-3, 3), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  x <- seq(-1, 1, length.out = 514)
  rect(x[-514], 0.55, x[-1], 0.85, col = pal, border = NA)
  text(c(-1, 0, 1), 0.2, c("-q99", "0", "q99"), cex = 0.85)
  dev.off()
}
