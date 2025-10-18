# Spectral clustering wrapper for clusGap:
# FUNcluster must be function(x, k) returning an object with $cluster

suppressPackageStartupMessages({
  library(Matrix); library(RSpectra)
})

spectral_k_wrapper <- function(d_threshold = 1.0, normalized = TRUE) {
  force(d_threshold); force(normalized)
  function(x, k) {
    x <- as.matrix(x)
    # Build epsilon-neighborhood graph
    D <- as.matrix(dist(x))
    A <- (D < d_threshold) * 1
    diag(A) <- 0

    # Degree and Laplacian
    deg <- rowSums(A)
    L <- diag(deg) - A
    if (normalized) {
      invsqrtD <- diag(1 / sqrt(pmax(deg, 1e-8)))
      L <- invsqrtD %*% L %*% invsqrtD
    }

    # Smallest k eigenvectors
    k_eff <- min(k, nrow(L))
    eig <- RSpectra::eigs(L, k = k_eff, which = "SR")
    U <- Re(eig$vectors)

    # Row normalization, then kmeans
    U <- sweep(U, 1, sqrt(rowSums(U^2)) + 1e-8, "/")
    km <- kmeans(U, centers = k, nstart = 20, iter.max = 50)
    km
  }
}
