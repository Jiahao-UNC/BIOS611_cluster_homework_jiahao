#!/usr/bin/env Rscript
# For each (max_radius) shell dataset:
# - Run clusGap with spectral clustering (epsilon graph threshold = 1.0)
# - Select k by Tibshirani's rule
# - Save summary table, failure points (first R where k_hat < K_true), and a figure.

suppressPackageStartupMessages({
  library(optparse); library(readr); library(dplyr); library(ggplot2); library(tidyr)
  source("scripts/utils_gap_helpers.R")
  source("scripts/spectral_clustering.R")
})

option_list <- list(
  make_option("--in",  type="character", help="Input RDS from task2_generate"),
  make_option("--tab", type="character", help="Output CSV summary"),
  make_option("--fig", type="character", help="Output figure path"),
  make_option("--fail_tab", type="character", default="results/tables/t2_failure_points.csv",
              help="Output CSV for first-underestimate points"),
  make_option("--kmax_cap", type="integer", default=10, help="K.max for clusGap"),
  make_option("--B", type="integer", default=50, help="Bootstrap samples for clusGap"),
  make_option("--d_threshold", type="double", default=1.0, help="Epsilon threshold for graph")
)
opt <- parse_args(OptionParser(option_list=option_list))
dir.create(dirname(opt$tab), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(opt$fig), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(opt$fail_tab), showWarnings = FALSE, recursive = TRUE)

# Load datasets
conds <- readRDS(opt$`in`)

# Spectral FUN for clusGap
FUNs <- spectral_k_wrapper(d_threshold = opt$d_threshold, normalized = TRUE)

summaries <- lapply(names(conds), function(key) {
  obj <- conds[[key]]
  X   <- obj$X
  R   <- obj$max_radius
  K_true <- obj$K_true

  Kmax <- min(K_true + 3L, opt$kmax_cap)

  set.seed(611)
  gap <- run_gap(X, FUNcluster = FUNs, k.max = Kmax, B = opt$B, verbose = FALSE)
  df  <- tidy_gap(gap)
  k_hat <- select_k_gap(df)

  tibble(
    key = key,
    max_radius = R,
    K_true = K_true,
    Kmax = Kmax,
    k_hat = k_hat,
    underestimate = as.integer(k_hat < K_true),
    diff_true_minus_hat = K_true - k_hat,
    d_threshold = opt$d_threshold
  )
})

tab <- bind_rows(summaries) %>%
  arrange(desc(max_radius))
readr::write_csv(tab, opt$tab)

# Failure points scanning R from large to small
fail_df <- tab %>%
  arrange(desc(max_radius)) %>%
  summarize(
    R_fail = {
      idx <- which(underestimate == 1L)
      if (length(idx) == 0) NA_real_ else max_radius[min(idx)]
    },
    .by = d_threshold
  )
readr::write_csv(fail_df, opt$fail_tab)

# Plot
p <- ggplot(tab, aes(x = max_radius, y = k_hat)) +
  geom_line() + geom_point(aes(color = underestimate == 1L)) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red"), guide = "none") +
  geom_hline(aes(yintercept = K_true), linetype = "dashed", color = "grey40") +
  geom_vline(
    data = fail_df %>% filter(!is.na(R_fail)),
    aes(xintercept = R_fail), color = "red", linetype = "dotted", linewidth = 0.6
  ) +
  geom_text(
    data = fail_df %>% filter(!is.na(R_fail)),
    aes(x = R_fail, y = Inf, label = paste0("R_fail=", R_fail)),
    vjust = 1.2, color = "red", size = 3
  ) +
  scale_x_reverse(breaks = sort(unique(tab$max_radius), decreasing = TRUE)) +
  labs(
    title = "Task 2: Estimated k vs Max Radius (Spectral + Gap)",
    subtitle = paste0("d_threshold = ", unique(tab$d_threshold),
                      " | K_true = ", unique(tab$K_true)),
    x = "Max radius (outer shell)",
    y = "Estimated k (k_hat)"
  ) +
  theme_minimal(base_size = 12)

ggsave(opt$fig, p, width = 10, height = 6, dpi = 300)

cat("Task2 summary saved ->", opt$tab, "\n")
cat("Task2 failure points saved ->", opt$fail_tab, "\n")
cat("Task2 figure  saved ->", opt$fig, "\n")
