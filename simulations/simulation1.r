simulation1 <- function(seed = 2026) {
  set.seed(seed)

  coords <- expand.grid(y = 1:50, x = 1:50)[, c("x", "y")]
  rownames(coords) <- paste0("S", 1:2500)
  labels <- cut(coords$x, breaks = 7, labels = paste0("L", 1:7),
                include.lowest = TRUE)
  names(labels) <- rownames(coords)

  X <- matrix(0, 2500, 300)
  j <- 0

  for (k in 1:7) {
    mu <- ifelse(labels == paste0("L", k), 2, 1)

    for (i in 1:15) {
      j <- j + 1
      X[, j] <- rnbinom(2500, mu = mu, size = 0.5)
    }
  }

  for (j in 106:300) {
    X[, j] <- rnbinom(2500, mu = 1, size = 0.5)
  }

  rownames(X) <- rownames(coords)
  colnames(X) <- c(
    paste0("M_L", rep(1:7, each = 15), "_", rep(1:15, 7)),
    paste0("HK_", 1:195)
  )

  list(
    counts = X,
    coords = coords,
    labels = labels,
    params = list(seed = seed)
  )
}
