#!/usr/bin/env Rscript
# Install minimal packages used in the project (add more later)

pkgs <- c(
  "optparse","readr","dplyr","ggplot2","tibble",
  "cluster","Matrix","igraph","RSpectra","scales"
)

# Install missing packages
need <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
if (length(need)) {
  install.packages(need, repos="https://cloud.r-project.org")
}

message("Packages ready: ", paste(pkgs, collapse=", "))
