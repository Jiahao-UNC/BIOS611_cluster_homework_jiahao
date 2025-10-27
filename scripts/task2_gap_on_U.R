#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse); library(readr); library(dplyr); library(ggplot2)
  source("scripts/utils_gap_helpers.R")
  source("scripts/spectral_clustering.R")
})
select_k_gap_elbow <- function(df){
  dg <- c(NA_real_, diff(df$gap))
  for(i in 2:nrow(df)) if(is.finite(dg[i]) && dg[i] <= df$se_sim[i]) return(df$k[i-1])
  df$k[which.max(df$gap)]
}
option_list <- list(
  make_option("--in",  type="character"),
  make_option("--tab", type="character"),
  make_option("--fig", type="character"),
  make_option("--fail_tab", type="character", default="results/tables/t2_failure_points_onU.csv"),
  make_option("--kmax_cap", type="integer", default=5),
  make_option("--B", type="integer", default=200),
  make_option("--d_threshold", type="double", default=1.0),
  make_option("--normalized", type="logical", default=TRUE),
  make_option("--m_embed", type="integer", default=3),
  make_option("--selector", type="character", default="elbow") # elbow or tib
)
opt <- parse_args(OptionParser(option_list=option_list))
dir.create(dirname(opt$tab), showWarnings=FALSE, recursive=TRUE)
dir.create(dirname(opt$fig), showWarnings=FALSE, recursive=TRUE)
dir.create(dirname(opt$fail_tab), showWarnings=FALSE, recursive=TRUE)

pick_k <- if(tolower(opt$selector)=="tib") select_k_gap else select_k_gap_elbow
conds <- readRDS(opt$`in`)
summ <- lapply(names(conds), function(key){
  o <- conds[[key]]; X <- o$X; R <- o$max_radius; Kt <- o$K_true
  Kmax <- min(opt$kmax_cap, opt$m_embed)
  U <- spectral_embed(X, m=opt$m_embed, d_threshold=opt$d_threshold, normalized=opt$normalized)
  set.seed(611)
  gp <- run_gap(U, FUNcluster=function(z,k) kmeans(z, centers=k, nstart=20, iter.max=50),
                k.max=Kmax, B=opt$B, verbose=FALSE)
  df <- tidy_gap(gp)
  kh <- pick_k(df)
  tibble(key=key, max_radius=R, K_true=Kt, k_hat=kh, Kmax=Kmax,
         m_embed=opt$m_embed, selector=tolower(opt$selector),
         d_threshold=opt$d_threshold, underestimate=as.integer(kh < Kt))
})
tab <- dplyr::bind_rows(summ) |> dplyr::arrange(dplyr::desc(max_radius))
readr::write_csv(tab, opt$tab)

fail_df <- tab |> dplyr::arrange(dplyr::desc(max_radius)) |>
  dplyr::summarize(R_fail = {i<-which(underestimate==1L); if(length(i)==0) NA_real_ else max_radius[min(i)]},
                   .by = c(d_threshold, selector, m_embed))
readr::write_csv(fail_df, opt$fail_tab)

p <- ggplot(tab, aes(max_radius, k_hat)) +
  geom_line() + geom_point(aes(color=underestimate==1L)) +
  scale_color_manual(values=c("FALSE"="black","TRUE"="red"), guide="none") +
  geom_hline(aes(yintercept=K_true), linetype="dashed", color="grey40") +
  geom_vline(data=fail_df |> dplyr::filter(!is.na(R_fail)),
             aes(xintercept=R_fail), color="red", linetype="dotted") +
  scale_x_reverse(breaks=sort(unique(tab$max_radius), decreasing=TRUE)) +
  labs(title="Task 2: Gap on spectral embedding U",
       subtitle=paste0("d_th=", unique(tab$d_threshold),
                       " | m=", unique(tab$m_embed),
                       " | K.max=", unique(tab$Kmax),
                       " | selector=", unique(tab$selector),
                       " | B=", unique(tab$B)),
       x="Max radius", y="Estimated k (k_hat)") +
  theme_minimal(base_size=12)
ggsave(opt$fig, p, width=10, height=6, dpi=300)

cat("Task2-onU summary ->", opt$tab, "\n")
cat("Task2-onU failure ->", opt$fail_tab, "\n")
cat("Task2-onU figure  ->", opt$fig, "\n")
