#!/usr/bin/env Rscript
# Generate concentric-shell (ring) datasets in 2D.
# - n_shells rings; each ring is a true cluster (class label = shell id)
# - For shell s in {1..n_shells}, radius r_s = s * (max_radius / n_shells)
# - Points sampled around each ring with small radial/angular noise.

suppressPackageStartupMessages({
  library(optparse); library(tibble)
})

option_list <- list(
  make_option("--out", type="character", help="Output RDS path"),
  make_option("--n_shells", type="integer", default=4, help="Number of shells (true clusters)"),
  make_option("--k_per_shell", type="integer", default=100, help="Points per shell"),
  make_option("--max_from", type="double", default=10, help="Max radius (start)"),
  make_option("--max_to",   type="double", default=1,  help="Max radius (end)"),
  make_option("--step",     type="double", default=1,  help="Radius step (decrease)"),
  make_option("--noise_sd", type="double", default=0.1, help="Radial noise sd")
)
opt <- parse_args(OptionParser(option_list=option_list))
dir.create(dirname(opt$out), showWarnings = FALSE, recursive = TRUE)

# Build a single 2D shell dataset at given max_radius
gen_shell_once <- function(n_shells, k_per_shell, max_radius, noise_sd) {
  # Protect against too small radius
  Rmax <- max(max_radius, 1e-6)
  r_base <- Rmax / n_shells
  X <- matrix(NA_real_, nrow = n_shells * k_per_shell, ncol = 2)
  y <- integer(n_shells * k_per_shell)
  idx <- 1L
  for (s in seq_len(n_shells)) {
    r <- s * r_base
    theta <- runif(k_per_shell, 0, 2*pi)
    # Radial jitter (Gaussian), angular jitter (small)
    rr <- r + rnorm(k_per_shell, 0, noise_sd)
    xx <- rr * cos(theta)
    yy <- rr * sin(theta)
    X[idx:(idx + k_per_shell - 1), ] <- cbind(xx, yy)
    y[idx:(idx + k_per_shell - 1)] <- s
    idx <- idx + k_per_shell
  }
  list(X = X, y = y, n_shells = n_shells, max_radius = max_radius, K_true = n_shells)
}

# Sweep max_radius from max_from down to max_to
r_vals <- seq(opt$max_from, opt$max_to, by = -opt$step)
if (tail(r_vals, 1) != opt$max_to) r_vals <- c(r_vals, opt$max_to)

set.seed(611)
conds <- list()
for (R in r_vals) {
  key <- paste0("R", gsub("\\.", "_", as.character(R)))
  conds[[key]] <- gen_shell_once(
    n_shells = opt$n_shells,
    k_per_shell = opt$k_per_shell,
    max_radius = R,
    noise_sd = opt$noise_sd
  )
}

saveRDS(conds, opt$out)
cat("Task2 simulations saved ->", opt$out, "\n")
