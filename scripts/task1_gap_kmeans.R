#!/usr/bin/env Rscript
# For each (n, side_length) dataset:
# - Run clusGap with K-means (nstart=20, iter.max=50)
# - Select k via Tibshirani's rule
# - Save a summary table and a figure: k_hat vs side_length (faceted by n)

suppressPackageStartupMessages({
  library(optparse); library(readr); library(dplyr); library(ggplot2)
  source("scripts/utils_gap_helpers.R")
})

option_list <- list(
  make_option("--in",  type="character", help="Input RDS from task1_generate"),
  make_option("--tab", type="character", help="Output CSV summary"),
  make_option("--fig", type="character", help="Output figure path"),
  make_option("--kmax_cap", type="integer", default=12,
              help="Cap for K.max in clusGap to control runtime"),
  make_option("--B", type="integer", default=50,
              help="Number of bootstrap samples for clusGap")
)
opt <- parse_args(OptionParser(option_list=option_list))
dir.create(dirname(opt$tab), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(opt$fig), showWarnings = FALSE, recursive = TRUE)

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
    k_hat = k_hat
  )
})

tab <- bind_rows(summaries) %>%
  arrange(desc(n), desc(side_length))

readr::write_csv(tab, opt$tab)

# ---- plot k_hat vs side_length, facet by n ----
p <- ggplot(tab, aes(x=side_length, y=k_hat, group=n)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ n, scales="free_y") +
  geom_hline(aes(yintercept=K_true), linetype="dashed", color="grey40") +
  scale_x_reverse(breaks=seq(max(tab$side_length), min(tab$side_length), by=-1)) +
  labs(
    title = "Task 1: Estimated number of clusters (k_hat) vs Side length",
    subtitle = "K-means (nstart=20, iter.max=50) + clusGap (Tibshirani's rule)",
    x = "Side length",
    y = "Estimated k (k_hat)"
  ) +
  theme_minimal(base_size = 12)

ggsave(opt$fig, p, width=10, height=6, dpi=300)
cat("Task1 summary saved ->", opt$tab, "\n")
cat("Task1 figure  saved ->", opt$fig, "\n")
