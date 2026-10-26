# Install every R package used by the analysis scripts.
# Run once with:  Rscript install_r_packages.R

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

cran_packages <- c(
  "data.table", "dplyr", "tidyr", "tibble", "stringr", "readr",
  "ggplot2", "ggrepel", "ggpubr", "ggvenn", "ggraph", "igraph",
  "pheatmap", "corrplot", "RColorBrewer", "gplots", "UpSetR",
  "rstatix", "rio", "scales"
)

bioc_packages <- c(
  "GEOquery", "limma", "WGCNA", "GSVA",
  "clusterProfiler", "org.Hs.eg.db", "STRINGdb"
)

for (p in cran_packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    message("Installing CRAN package: ", p)
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}

for (p in bioc_packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    message("Installing Bioconductor package: ", p)
    BiocManager::install(p, ask = FALSE, update = FALSE)
  }
}

message("All packages installed.")
sessionInfo()
