#!/usr/bin/env Rscript
# Generate concentric-shell (2D rings) datasets with optional density compensation.
# - If compensate = "none": each shell has exactly k_per_shell points (original behavior).
# - If compensate = "ring": points per shell scale ~ r (constant angular density on a circle).
# - If compensate = "sphere": points per shell scale ~ r^2 (constant density on a sphere surface).
#
# English comments only in code.

suppressPackageStartupMessages({
  library(optparse); library(tibble)
})

option_list <- list(
  make_option("--out", type="character", help="Output RDS path"),
  make_option("--n_shells", type="integer", default=4, help="Number of shells (true clusters)"),
  make_option("--k_per_shell", type="integer", default=100,
              help="Points per shell when compensate='none'"),
  make_option("--k_per_base", type="integer", default=100,
              help="Target points on the OUTER-most shell when compensation is on"),
  make_option("--min_per_shell", type="integer", default=40,
              help="Minimum points per shell under compensation"),
  make_option("--compensate", type="character", default="none",
              help="Density compensation: 'none' | 'ring' | 'sphere'"),
  make_option("--max_from", type="double", default=10, help="Max radius (start)"),
  make_option("--max_to",   type="double", default=1,  help="Max radius (end)"),
  make_option("--step",     type="double", default=1,  help="Radius step (decrease)"),
  make_option("--noise_sd", type="double", default=0.1, help="Radial noise sd")
)
opt <- parse_args(OptionParser(option_list=option_list))
dir.create(dirname(opt$out), showWarnings = FALSE, recursive = TRUE)

# Decide exponent for compensation
get_exp <- function(mode) {
  mode <- tolower(mode)
  if (mode == "ring") return(1)     # circumference ~ r
  if (mode == "sphere") return(2)   # surface area ~ r^2
  return(NA_real_)                   # none
}

# Build a single 2D ring dataset at given max_radius
gen_shell_once <- function(n_shells, max_radius, noise_sd,
                           compensate = "none", k_per_shell = 100,
                           k_per_base = 100, min_per_shell = 40) {
  Rmax <- max(max_radius, 1e-8)
  r_base <- Rmax / n_shells

  # Decide per-shell counts
  exp <- get_exp(compensate)
  if (is.na(exp)) {
    # No compensation: equal counts
    k_vec <- rep.int(k_per_shell, n_shells)
  } else {
    # Compensation: points per shell scale with (r / Rmax)^exp
    k_vec <- integer(n_shells)
    for (s in seq_len(n_shells)) {
      r_s <- s * r_base
      if (Rmax <= 1e-8) {
        k_s <- k_per_base
      } else {
        ratio <- r_s / Rmax           # in (0,1]
        k_s <- round(k_per_base * (ratio ^ exp))
      }
      k_vec[s] <- max(min_per_shell, k_s)
    }
  }

  total_n <- sum(k_vec)
  X <- matrix(NA_real_, nrow = total_n, ncol = 2)
  y <- integer(total_n)
  idx <- 1L
  for (s in seq_len(n_shells)) {
    r <- s * r_base
    k_s <- k_vec[s]
    theta <- runif(k_s, 0, 2*pi)
    rr <- r + rnorm(k_s, 0, noise_sd)
    X[idx:(idx + k_s - 1), ] <- cbind(rr * cos(theta), rr * sin(theta))
    y[idx:(idx + k_s - 1)] <- s
    idx <- idx + k_s
  }

  list(X = X, y = y, n_shells = n_shells,
       max_radius = max_radius, K_true = n_shells,
       k_per_shell_vec = k_vec,
       compensate = compensate)
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
    max_radius = R,
    noise_sd = opt$noise_sd,
    compensate = opt$compensate,
    k_per_shell = opt$k_per_shell,
    k_per_base = opt$k_per_base,
    min_per_shell = opt$min_per_shell
  )
}

saveRDS(conds, opt$out)
cat("Task2 simulations saved ->", opt$out, "\n")
cat("Compensation mode:", opt$compensate, "| outer-shell target:", opt$k_per_base,
    "| min_per_shell:", opt$min_per_shell, "\n")
