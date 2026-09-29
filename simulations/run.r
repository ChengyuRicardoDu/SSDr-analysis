source("simulation1.r")
source("simulation2.r")
ref <- readRDS("reference.rds")
seeds <- read.csv("seeds.csv")

reps <- 1
seeds <- seeds[seeds$rep %in% reps, ]

dir.create("data", showWarnings = FALSE)

for (i in seq_len(nrow(seeds))) {
  s <- seeds[i, ]
  id <- if (s$simulation == "simulation1") "s1" else "s2"
  file <- sprintf("data/%s_evaluation_%03d.rds", id, s$rep)
  if (file.exists(file)) stop("File already exists: ", file)

  if (s$simulation == "simulation1") {
    data <- simulation1(s$seed1)
  } else {
    data <- simulation2(s$seed1, s$seed2, ref)
  }

  saveRDS(data, file)
  message(file)
}
