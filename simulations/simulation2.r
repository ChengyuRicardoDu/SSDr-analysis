simulation2 <- function(seed1 = 2031, seed2 = 2032, ref) {
  pool <- as.character(ref$meta$label_factor)

  coords <- expand.grid(y = 1:50, x = 1:50)[, c("x", "y")]
  rownames(coords) <- paste0("S", 1:2500)
  groups <- c("Layer1", "Layer3", "Layer5", "WM")

  set.seed(seed1)
  centers <- rbind(c(12.5, 12.5), c(37.5, 12.5),
                   c(12.5, 37.5), c(37.5, 37.5))
  centers <- centers + matrix(runif(8, -2.5, 2.5), 4, 2)
  d <- sapply(1:4, function(k)
    (coords$x - centers[k, 1])^2 + (coords$y - centers[k, 2])^2)
  labels <- groups[max.col(-d, ties.method = "first")]

  set.seed(seed2)
  index <- integer(2500)
  for (i in 1:2500) {
    j <- which(pool == labels[i])
    index[i] <- sample(j, 1)
  }

  X <- ref$X[index, , drop = FALSE]
  rownames(X) <- rownames(coords)
  labels <- factor(labels, levels = groups)
  names(labels) <- rownames(coords)

  list(
    counts = X,
    coords = coords,
    labels = labels,
    params = list(geometry_seed = seed1, bootstrap_seed = seed2)
  )
}
