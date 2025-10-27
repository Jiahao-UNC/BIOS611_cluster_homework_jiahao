# Spectral embedding & clustering helpers
# - spectral_embed(): build epsilon-graph -> normalized Laplacian -> smallest m eigenvectors
# - robust fallbacks so clusGap won't crash on degenerate graphs

suppressPackageStartupMessages({
  library(Matrix)
  library(RSpectra)
})

# ---- Spectral Embedding (for Gap on U) ----
spectral_embed <- function(x, m, d_threshold = 1.0, normalized = TRUE) {
  # x: numeric matrix (n x p)
  # m: embedding dimension (we'll set m >= K.max)
  x <- as.matrix(x)
  n <- nrow(x)
  if (n < 2) return(matrix(0, nrow = n, ncol = max(1L, m)))

  # Build epsilon-neighborhood graph
  D <- as.matrix(dist(x))
  A <- (D < d_threshold) * 1
  diag(A) <- 0

  # If there are no edges, return original x (so clusGap can still run)
  if (sum(A) == 0) return(x[, seq_len(min(ncol(x), m)), drop = FALSE])

  # Normalized Laplacian
  deg <- rowSums(A)
  L <- diag(deg) - A
  if (normalized) {
    invsqrtD <- diag(1 / sqrt(pmax(deg, 1e-12)))
    L <- invsqrtD %*% L %*% invsqrtD
  }

  m_eff <- min(max(2L, m), max(2L, n - 1L))

  # Prefer symmetric eigensolver on Laplacian
  eig <- tryCatch(
    RSpectra::eigs_sym(L, k = m_eff, which = "SM"),
    error = function(e) NULL,
    warning = function(w) NULL
  )

  U <- if (is.null(eig) || any(!is.finite(eig$values))) {
    ev <- tryCatch(eigen(L, symmetric = TRUE), error = function(e) NULL)
    if (is.null(ev)) {
      return(x[, seq_len(min(ncol(x), m_eff)), drop = FALSE])
    }
    ev$vectors[, seq_len(m_eff), drop = FALSE]
  } else {
    Re(eig$vectors)
  }

  # Row-normalize + tiny jitter
  rn <- sqrt(rowSums(U^2))
  U <- sweep(U, 1, rn + 1e-12, "/")
  U[!is.finite(U)] <- 0
  U <- U + matrix(rnorm(n * ncol(U), 0, 1e-8), nrow = n)

  U
}

# (Optional) The old "direct spectral clustering" wrapper is kept for reference/testing.
spectral_k_wrapper <- function(d_threshold = 1.0, normalized = TRUE) {
  force(d_threshold); force(normalized)
  function(x, k) {
    # embed into k dims then kmeans
    U <- spectral_embed(x, m = max(2L, k), d_threshold = d_threshold, normalized = normalized)
    kmeans(U, centers = k, nstart = 20, iter.max = 50)
  }
}
