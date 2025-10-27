#!/usr/bin/env Rscript
# Task 2: clusGap ON A SPECTRAL-CLUSTERING FUN (embed each dataset, incl. references).
# Parameters per assignment: n_shells=4, k_per_shell=100, noise_sd=0.1, d_threshold=1.0,
# max_radius from 10 down to 0.

suppressPackageStartupMessages({
  library(optparse); library(readr); library(dplyr); library(ggplot2); library(tidyr)
  source("scripts/utils_gap_helpers.R")
  source("scripts/spectral_clustering.R")   # spectral_embed()
})

# Elbow selector (backup to Tibshirani)
select_k_gap_elbow <- function(df) {
  dg <- c(NA_real_, diff(df$gap))
  for (i in 2:nrow(df)) if (is.finite(dg[i]) && dg[i] <= df$se_sim[i]) return(df$k[i-1])
  df$k[which.max(df$gap)]
}

option_list <- list(
  make_option("--in",  type="character", help="Input RDS from task2_generate"),
  make_option("--tab", type="character", help="Output CSV summary"),
  make_option("--fig", type="character", help="Output figure path"),
  make_option("--fail_tab", type="character", default="results/tables/t2_failure_points.csv",
              help="Output CSV for first-underestimate points"),
  make_option("--kmax_cap", type="integer", default=6, help="K.max for clusGap (<= m_embed)"),
  make_option("--B", type="integer", default=100, help="Bootstrap samples for clusGap"),
  make_option("--d_threshold", type="double", default=1.0, help="Epsilon threshold"),
  make_option("--normalized", type="logical", default=TRUE, help="Use normalized Laplacian"),
  make_option("--m_embed", type="integer", default=4, help="Fixed embedding dimension (≈ true k)"),
  make_option("--selector", type="character", default="tib", help="tib|elbow")
)
opt <- parse_args(OptionParser(option_list=option_list))
dir.create(dirname(opt$tab), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(opt$fig), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(opt$fail_tab), showWarnings = FALSE, recursive = TRUE)

conds <- readRDS(opt$`in`)
pick_k <- if (tolower(opt$selector) == "elbow") select_k_gap_elbow else select_k_gap

# --- define FUNcluster: spectral embed + kmeans, applied to ANY dataset X (incl. references) ---
spectral_FUN <- function(X, k) {
  U <- spectral_embed(X, m = opt$m_embed, d_threshold = opt$d_threshold, normalized = opt$normalized)
  kmeans(U, centers = k, nstart = 20, iter.max = 50)
}

summaries <- lapply(names(conds), function(key) {
  o <- conds[[key]]
  X <- o$X; R <- o$max_radius; K_true <- o$K_true
  Kmax <- min(opt$kmax_cap, opt$m_embed)

  set.seed(611)
  gp <- run_gap(
    X,
    FUNcluster = spectral_FUN,   # <<=== IMPORTANT: Fun embeds both real AND reference datasets
    k.max = Kmax, B = opt$B, verbose = FALSE
  )
  df <- tidy_gap(gp)
  k_hat <- pick_k(df)

  tibble(
    key=key, max_radius=R, K_true=K_true, Kmax=Kmax,
    m_embed=opt$m_embed, selector=tolower(opt$selector),
    d_threshold=opt$d_threshold, k_hat=k_hat,
    underestimate=as.integer(k_hat < K_true)
  )
})

tab <- bind_rows(summaries) %>% arrange(desc(max_radius))
write_csv(tab, opt$tab)

fail_df <- tab %>%
  arrange(desc(max_radius)) %>%
  summarize(R_fail = { i<-which(underestimate==1L); if(length(i)==0) NA_real_ else max_radius[min(i)] },
            .by = c(d_threshold, selector, m_embed))
write_csv(fail_df, opt$fail_tab)

p <- ggplot(tab, aes(max_radius, k_hat)) +
  geom_line() + geom_point(aes(color = underestimate==1L)) +
  scale_color_manual(values=c("FALSE"="black","TRUE"="red"), guide="none") +
  geom_hline(aes(yintercept=K_true), linetype="dashed", color="grey40") +
  geom_vline(data=fail_df %>% filter(!is.na(R_fail)), aes(xintercept=R_fail),
             color="red", linetype="dotted", linewidth=0.6) +
  geom_text(data=fail_df %>% filter(!is.na(R_fail)),
            aes(x=R_fail, y=Inf, label=paste0("R_fail=", R_fail)),
            vjust=1.2, color="red", size=3) +
  scale_x_reverse(breaks=sort(unique(tab$max_radius), decreasing=TRUE)) +
  labs(title="Task 2: clusGap on spectral clustering (ε-graph)",
       subtitle=paste0("d_th=", unique(tab$d_threshold),
                       " | m_embed=", unique(tab$m_embed),
                       " | selector=", unique(tab$selector)),
       x="Max radius", y="Estimated k (k_hat)") +
  theme_minimal(base_size=12)
ggsave(opt$fig, p, width=10, height=6, dpi=300)

cat("Task2 summary saved ->", opt$tab, "\n")
cat("Task2 failure points saved ->", opt$fail_tab, "\n")
cat("Task2 figure  saved ->", opt$fig, "\n")
