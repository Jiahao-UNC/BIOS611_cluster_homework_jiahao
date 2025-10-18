#!/usr/bin/env Rscript
# Generate hypercube-cluster datasets across (n, side_length)
# - Each vertex of an n-D hypercube (0 or side_length per axis) is a cluster center.
# - For each center, sample k_per_cluster points with isotropic Gaussian noise (sd = noise_sd).
# - Save a list of conditions; each element contains X (matrix) and y (labels).

suppressPackageStartupMessages({
  library(optparse); library(tibble); library(purrr)
})

# ---- CLI options ----
option_list <- list(
  make_option("--out", type="character", help="Output RDS path"),
  make_option("--ns", type="character", default="6,5,4,3,2",
              help="Comma-separated dimensions, e.g. '6,5,4,3,2'"),
  make_option("--side_from", type="integer", default=10, help="Max side length"),
  make_option("--side_to",   type="integer", default=1,  help="Min side length"),
  make_option("--side_step", type="integer", default=1,  help="Side length step"),
  make_option("--k_per_cluster", type="integer", default=100,
              help="Points per vertex/cluster"),
  make_option("--noise_sd",  type="double", default=1.0,
              help="Std. dev. for isotropic Gaussian noise")
)
opt <- parse_args(OptionParser(option_list=option_list))
dir.create(dirname(opt$out), showWarnings = FALSE, recursive = TRUE)

# ---- helpers ----
# All vertices of a hypercube with side length L in n-D: {0, L}^n
hypercube_vertices <- function(n, L) {
  # Return a matrix of size (2^n x n)
  B <- expand.grid(rep(list(c(0, L)), n))
  as.matrix(B)
}

generate_hypercube_clusters <- function(n, side_length, k_per_cluster, noise_sd) {
  # Build centers at hypercube vertices
  centers <- hypercube_vertices(n, side_length)
  K <- nrow(centers)  # = 2^n
  N <- K * k_per_cluster
  
  # Generate points around each center
  X <- matrix(NA_real_, nrow=N, ncol=n)
  y <- integer(N)
  idx <- 1L
  for (k in 1:K) {
    mu <- centers[k, ]
    block <- matrix(rnorm(k_per_cluster * n, mean=0, sd=noise_sd), ncol=n)
    block <- sweep(block, 2, mu, "+")
    X[idx:(idx + k_per_cluster - 1), ] <- block
    y[idx:(idx + k_per_cluster - 1)] <- k
    idx <- idx + k_per_cluster
  }
  list(X=X, y=y, n=n, side_length=side_length, K_true=K)
}

# ---- build condition grid ----
ns <- as.integer(strsplit(opt$ns, ",")[[1]])
side_vals <- seq(opt$side_from, opt$side_to, by=-opt$side_step)
if (side_vals[length(side_vals)] != opt$side_to) side_vals <- c(side_vals, opt$side_to)

# ---- generate for all (n, side_length) ----
set.seed(611)
conds <- list()
for (n in ns) {
  for (L in side_vals) {
    conds[[paste0("n", n, "_L", L)]] <-
      generate_hypercube_clusters(
        n=n,
        side_length=L,
        k_per_cluster=opt$k_per_cluster,
        noise_sd=opt$noise_sd
      )
  }
}

saveRDS(conds, opt$out)
cat("Task1 simulations saved ->", opt$out, "\n")
