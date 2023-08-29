pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", "optparse",
          "jhtools", "glue", "openxlsx", "ggsci", "tidyverse", "dplyr", "Seurat")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}


# what? no need to write scripts now




