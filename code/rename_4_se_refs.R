pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "tidyverse", "dplyr")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}

# set directory
se_dir <- "/cluster/home/jhuang/reference/scell/celldex/human"
file.rename(glue::glue("{se_dir}/ref_se_aml_30827681.rds"), 
            glue::glue("{se_dir}/aml_pmid30827681.rds"))
file.rename(glue::glue("{se_dir}/ref_se_fetal_bm_34588693.rds"), 
            glue::glue("{se_dir}/fetal_bm_pmid34588693.rds"))
file.rename(glue::glue("{se_dir}/ref_se_adult_bm_ERP122984.rds"), 
            glue::glue("{se_dir}/adult_bm_ERP122984.rds"))
file.rename(glue::glue("{se_dir}/ref_se_cord_ERP122984.rds"), 
            glue::glue("{se_dir}/cordblood_ERP122984.rds"))





