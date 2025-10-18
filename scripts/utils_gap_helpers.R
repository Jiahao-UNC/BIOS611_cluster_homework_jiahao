# Helpers for clusGap pipeline

suppressPackageStartupMessages({
  library(cluster)  # for clusGap
  library(tibble)   # for tibble return
})

# Run clusGap with a provided clustering function
run_gap <- function(x, FUNcluster, k.max = 10, B = 50, verbose = TRUE, ...) {
  # x: numeric matrix/data.frame
  # FUNcluster: function(x, k) -> object with $cluster
  set.seed(1)
  cluster::clusGap(x, FUNcluster = FUNcluster, K.max = k.max, B = B, verbose = verbose, ...)
}

# Tidy clusGap object into a tibble for plotting and selection
tidy_gap <- function(gap_obj) {
  # gap_obj: result of cluster::clusGap
  k <- seq_len(nrow(gap_obj$Tab))
  tab <- as.data.frame(gap_obj$Tab)
  tibble::tibble(
    k       = k,
    logW    = tab$logW,
    E_logW  = tab$E.logW,
    gap     = tab$gap,
    se_sim  = tab$SE.sim
  )
}

# Select k by Tibshirani's rule: gap(k) >= gap(k+1) - se(k+1)
select_k_gap <- function(df) {
  # df: tibble from tidy_gap()
  ks <- df$k
  for (i in seq_along(ks)) {
    k <- ks[i]
    if (k == max(ks)) return(k)
    if (df$gap[i] >= df$gap[i + 1] - df$se_sim[i + 1]) return(k)
  }
  return(max(ks))
}
