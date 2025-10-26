# BIOS611 Clustering Homework

This repository contains my solutions for Task 1 and Task 2 of the BIOS611 clustering assignment. Code is in R, and the project is reproducible via a Makefile.

## Task 1: Hypercube clusters + Gap Statistic
Goal: simulate n-dimensional hypercube datasets and estimate the number of clusters with the Gap Statistic.

Method:
1) Generate data where each vertex of an n-dimensional hypercube is a true cluster center (K_true = 2^n).
2) Sample 100 points per center with Gaussian noise.
3) Run K-means (nstart=20, iter.max=50) and compute the Gap Statistic (cluster::clusGap).
4) Sweep side length L from large to small and record k_hat versus L.

Outcome:
- A plot of k_hat vs L with a dashed reference line for K_true.
- A table and a failure-point summary showing where underestimation begins.

Main outputs:
- results/tables/t1_gap_summary.csv
- results/tables/t1_failure_points.csv
- results/figs/t1_gap_curve.png

## Task 2: Concentric shells + Spectral clustering + Gap
Goal: detect concentric shells using spectral clustering and locate the failure radius as shells collapse.

Method:
1) Build an epsilon-neighborhood graph (A_ij=1 if distance(i,j) < d_threshold).
2) Compute the symmetric normalized Laplacian L_sym = D^{-1/2}(D-A)D^{-1/2}.
3) Take the smallest eigenvectors to form a spectral embedding U (fixed dimension near the true K).
4) Run K-means on U. Use clusGap to estimate k_hat while scanning max_radius from 10 down to 0.
5) Compare several thresholds (d_threshold = 0.8, 1.0, 1.2) and report the first radius R_fail where k_hat < 4.

Outcome:
- For large radii, k_hat is about 4. As max_radius decreases, k_hat drops, revealing a failure point.
- Larger thresholds make the graph denser and cause earlier failure; smaller thresholds delay failure.
- Includes a small 3D preview (plotly) to visually confirm the shell structure before simulations.

Main outputs (examples):
- results/tables/t2_gap_summary_*.csv
- results/tables/t2_failure_points_*.csv
- results/figs/t2_gap_curve_*.png
- results/figs/t2_preview_3d.html
- results/figs/t2_preview_3d_snapshot.png
## Project structure
- scripts/: data generation, spectral embedding, Gap evaluation, plotting
- data/: raw and processed data
- results/: tables (.csv) and figures (.png/.html)
- report/: final Quarto report
- Makefile: one-click workflow

### Task 2 — Important caveat for 3D shells (fixed epsilon, fixed points per shell)

In the 3D concentric-shell experiment we keep the **distance threshold** \(d_{\text{th}}\) fixed while sweeping the **maximum radius**, and we also keep the **same number of points per shell**. The surface area of a sphere grows like \(4\pi r^2\), so as the radius increases, the **point density on a shell decreases** and the **typical within-shell neighbor spacing increases** (roughly \(\propto r/\sqrt{N}\) for \(N\) points). With a single global \(d_{\text{th}}\), the \(\varepsilon\)-neighborhood graph can therefore **become too sparse within a shell at large radii**, while at very small radii it can **become over-connected and form cross-shell bridges**. Either degeneracy makes the spectral embedding less informative and can drive the **Gap statistic to return \( \hat{k} = 1 \)** across the sweep.
