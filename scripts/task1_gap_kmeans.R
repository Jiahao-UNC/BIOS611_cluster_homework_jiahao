#!/usr/bin/env Rscript
# For each (n, side_length) dataset:
# - Run clusGap with K-means (nstart=20, iter.max=50)
# - Select k via Tibshirani's rule
# - Save:
#   (1) a summary table of k_hat across conditions
#   (2) a figure with k_hat vs side_length (faceted by n)
#   (3) a "failure points" table: first side_length where k_hat < K_true (per n)
#
# Note: code comments are in English per project convention.

suppressPackageStartupMessages({
  library(optparse); library(readr); library(dplyr); library(ggplot2); library(tidyr)
  source("scripts/utils_gap_helpers.R")
})

option_list <- list(
  make_option("--in",  type="character", help="Input RDS from task1_generate"),
  make_option("--tab", type="character", help="Output CSV summary"),
  make_option("--fig", type="character", help="Output figure path"),
  make_option("--fail_tab", type="character", default="results/tables/t1_failure_points.csv",
              help="Output CSV for first-underestimate points per n"),
  make_option("--kmax_cap", type="integer", default=12,
              help="Cap for K.max in clusGap to control runtime"),
  make_option("--B", type="integer", default=50,
              help="Number of bootstrap samples for clusGap")
)
opt <- parse_args(OptionParser(option_list=option_list))
dir.create(dirname(opt$tab), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(opt$fig), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(opt$fail_tab), showWarnings = FALSE, recursive = TRUE)

# ---- load datasets ----
conds <- readRDS(opt$`in`)

# ---- run clusGap for each condition ----
summaries <- lapply(names(conds), function(key) {
  obj <- conds[[key]]
  X   <- obj$X
  n   <- obj$n
  L   <- obj$side_length
  K_true <- obj$K_true
  
  # Choose K.max: up to min(2^n + 3, kmax_cap)
  Kmax <- min(K_true + 3L, opt$kmax_cap)
  
  # K-means wrapper with required settings
  k_fun <- function(xx, k) kmeans(xx, centers=k, nstart=20, iter.max=50)
  
  set.seed(611)
  gap <- run_gap(X, FUNcluster=k_fun, k.max=Kmax, B=opt$B, verbose=FALSE)
  df  <- tidy_gap(gap)
  k_hat <- select_k_gap(df)
  
  tibble(
    key = key,
    n = n,
    side_length = L,
    K_true = K_true,
    Kmax = Kmax,
    k_hat = k_hat,
    underestimate = as.integer(k_hat < K_true),
    diff_true_minus_hat = K_true - k_hat
  )
})

tab <- bind_rows(summaries) %>%
  arrange(desc(n), desc(side_length))

# Save summary table
readr::write_csv(tab, opt$tab)

# ---- compute failure points per n ----
# Definition: scanning side_length from large to small, the first L such that k_hat < K_true
fail_df <- tab %>%
  group_by(n) %>%
  arrange(desc(side_length), .by_group = TRUE) %>%
  summarize(
    L_fail = {
      idx <- which(underestimate == 1L)
      if (length(idx) == 0) NA_real_ else side_length[min(idx)]
    },
    .groups = "drop"
  )

readr::write_csv(fail_df, opt$fail_tab)

# ---- plot k_hat vs side_length, facet by n; add failure vertical line ----
p <- ggplot(tab, aes(x=side_length, y=k_hat, group=n)) +
  geom_line() +
  geom_point(aes(color = underestimate == 1L)) +
  scale_color_manual(values = c("FALSE"="black", "TRUE"="red"), guide = "none") +
  facet_wrap(~ n, scales="free_y") +
  geom_hline(aes(yintercept=K_true), linetype="dashed", color="grey40") +
  # add vertical lines per facet for L_fail
  geom_vline(
    data = fail_df %>% filter(!is.na(L_fail)),
    aes(xintercept = L_fail), color = "red", linetype = "dotted", linewidth = 0.6
  ) +
  # annotate text near the vertical line
  geom_text(
    data = fail_df %>% filter(!is.na(L_fail)),
    aes(x = L_fail, y = Inf, label = paste0("L_fail=", L_fail)),
    vjust = 1.2, color = "red", size = 3
  ) +
  scale_x_reverse(breaks=sort(unique(tab$side_length), decreasing = TRUE)) +
  labs(
    title = "Task 1: Estimated number of clusters (k_hat) vs Side length",
    subtitle = "K-means (nstart=20, iter.max=50) + clusGap (Tibshirani's rule)\nRed dotted line = first underestimation point (L_fail)",
    x = "Side length",
    y = "Estimated k (k_hat)"
  ) +
  theme_minimal(base_size = 12)

ggsave(opt$fig, p, width=10, height=6, dpi=300)

cat("Task1 summary saved ->", opt$tab, "\n")
cat("Task1 failure points saved ->", opt$fail_tab, "\n")
cat("Task1 figure  saved ->", opt$fig, "\n")
