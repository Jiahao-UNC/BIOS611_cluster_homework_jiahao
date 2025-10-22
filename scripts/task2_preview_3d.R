#!/usr/bin/env Rscript
# Generate ONE sample of concentric 3D shells and make an interactive 3D scatter (plotly)
# Also save a lightweight 2D snapshot (PCA projection) as PNG for README/report.

suppressPackageStartupMessages({
  library(optparse); library(dplyr); library(ggplot2)
  library(plotly);   library(htmlwidgets)
})

option_list <- list(
  make_option("--out_rds",  type="character", default="data/processed/t2_preview_3d.rds",
              help="Output RDS for the preview dataset"),
  make_option("--out_html", type="character", default="results/figs/t2_preview_3d.html",
              help="Output HTML (plotly widget)"),
  make_option("--out_png",  type="character", default="results/figs/t2_preview_3d_snapshot.png",
              help="Output PNG snapshot (2D PCA projection)"),
  make_option("--n_shells", type="integer", default=4,   help="Number of shells"),
  make_option("--k_per",    type="integer", default=500, help="Points per shell"),
  make_option("--max_r",    type="double",  default=10,  help="Max (outer) radius"),
  make_option("--noise_sd", type="double",  default=0.05,help="Radial noise sd")
)
opt <- parse_args(OptionParser(option_list=option_list))

dir.create(dirname(opt$out_rds),  showWarnings=FALSE, recursive=TRUE)
dir.create(dirname(opt$out_html), showWarnings=FALSE, recursive=TRUE)
dir.create(dirname(opt$out_png),  showWarnings=FALSE, recursive=TRUE)

# ---- helper: random unit vectors in R^3 ----
rand_unit_vec3 <- function(n){
  M <- matrix(rnorm(3*n), ncol=3)
  M <- sweep(M, 1, sqrt(rowSums(M^2)), "/")
  M
}

# ---- build one 3D sample of concentric spherical shells ----
set.seed(611)
n_shells <- opt$n_shells
k_per    <- opt$k_per
max_r    <- opt$max_r
noise_sd <- opt$noise_sd

r_base <- max_r / n_shells
N      <- n_shells * k_per
XYZ    <- matrix(NA_real_, nrow=N, ncol=3)
lab    <- integer(N)

idx <- 1L
for (s in seq_len(n_shells)) {
  r_s <- s * r_base
  U   <- rand_unit_vec3(k_per)
  rr  <- r_s + rnorm(k_per, 0, noise_sd)
  pts <- sweep(U, 1, rr, "*")
  XYZ[idx:(idx+k_per-1),] <- pts
  lab[idx:(idx+k_per-1)]  <- s
  idx <- idx + k_per
}

df <- as.data.frame(XYZ); names(df) <- c("x","y","z")
df$shell <- factor(lab)

saveRDS(df, opt$out_rds)
cat("Preview dataset saved ->", opt$out_rds, "\n")

# ---- interactive 3D (plotly) ----
plt3d <- plot_ly(
  df, x = ~x, y = ~y, z = ~z,
  color = ~shell, colors = RColorBrewer::brewer.pal(n_shells, "Set1"),
  type = "scatter3d", mode = "markers",
  marker = list(size = 2, opacity = 0.7)
) %>%
  layout(
    scene = list(
      xaxis = list(title="x"), yaxis = list(title="y"), zaxis = list(title="z"),
      aspectmode = "data"
    ),
    legend = list(title=list(text="shell"))
  )

saveWidget(plt3d, opt$out_html, selfcontained = FALSE, libdir = "t2_preview_3d_files")
cat("Interactive 3D plot saved ->", opt$out_html, "\n")

# ---- static PNG: PCA projection (df -> PC1/PC2) ----
pr <- prcomp(df[,c("x","y","z")], scale.=FALSE)
p2 <- data.frame(pc1 = pr$x[,1], pc2 = pr$x[,2], shell = df$shell)
g <- ggplot(p2, aes(pc1, pc2, color = shell)) +
  geom_point(alpha=0.6, size=0.6) +
  coord_equal() + theme_minimal(base_size=12) +
  labs(title = "Concentric 3D shells (PCA projection)", x="PC1", y="PC2", color="shell")
ggsave(filename = opt$out_png, plot = g, width = 7, height = 5, dpi = 300)
cat("Static snapshot saved ->", opt$out_png, "\n")
