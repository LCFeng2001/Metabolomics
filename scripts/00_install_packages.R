## =========================================================
## 00_install_packages.R
## v3: add ComplexHeatmap / circlize for publication heatmaps
## =========================================================

cran_pkgs <- c(
  "tidyverse",
  "readxl",
  "readr",
  "openxlsx",
  "ggplot2",
  "ggrepel",
  "ggpubr",
  "svglite",
  "RColorBrewer",
  "scales",
  "patchwork"
)

bioc_pkgs <- c(
  "limma",
  "ropls",
  "ComplexHeatmap",
  "circlize",
  "clusterProfiler",
  "KEGGREST"
)

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

message("Package check finished.")
