pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr", "igraph",
          "Seurat", "anndata", "stats", "pheatmap", "pracma", "SeuratDisk"
)
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}



