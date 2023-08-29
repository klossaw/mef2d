pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr",
          "SingleR", "SummarizedExperiment", "Matrix"
)
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}

# first, prepare a count matrix (dgcmatrix) in which rows represent genes
# second, prepare corresponding annotation data which row number equals to 
# column number (cells) in above matrix
# third, normalize the counts for singler input

ref_se_gen <- function(ref_mat, anno_data){
  se <- SummarizedExperiment(assays = list(counts = ref_mat), colData = anno_data)
  se_norm <- scuttle::logNormCounts(se)
  se_norm
}

# (data generation: see 0_ref_data_preprocess.R)
# ======= generate SummarizedExperiment object for SingleR ========
data_dir <- "/cluster/home/yjliu_jh/projects/mef2d/data"

# ======= the leukemia dataset collected from GSE116256 ========
#leu1_mat_s <- readr::read_rds(glue::glue("{data_dir}/public/rds/leu1_assay_sub.rds"))
leu1_mat <- readr::read_rds(glue::glue("{data_dir}/public/rds/ref_cells_lsc_mat.rds"))
leu1_anno <- readr::read_rds(glue::glue("{data_dir}/ref_cells_lsc_anno.rds"))
#leu1_ref_small <- ref_se_gen(leu1_mat_s, leu1_anno)
leu1_ref <- ref_se_gen(leu1_mat, leu1_anno)
#readr::write_rds(leu1_ref_small, glue::glue("{data_dir}/ref_se_aml_30827681_small.rds"))
readr::write_rds(leu1_ref, glue::glue("{data_dir}/ref_se_aml_30827681.rds"))
# in the official SingleR example, the reference contained 19000+ genes
# one may test the results based on the matrix designed for other tools (filtered)


# ======= fetal bone marrow data ========
fbm <- readr::read_rds(glue::glue("{data_dir}/h5ad1_mat_dgc.rds"))
fbm_c <- t(fbm)
pdata_fbm <- readr::read_rds(glue::glue("{data_dir}/h5ad1_obs.rds"))
sce_fbm_norm <-  ref_se_gen(fbm_c, pdata_fbm)
readr::write_rds(sce_fbm_norm, glue::glue("{data_dir}/ref_se_fetal_bm_34588693.rds"))

# ======= adult bone marrow data ========
abm <- readr::read_rds(glue::glue("{data_dir}/public/rds/adultbm_cells_matrix_sub.rds"))
abm_c <- t(abm)
pdata_abm <- readr::read_rds(glue::glue("{data_dir}/public/rds/adultbm_cells_metadata.rds"))
sce_abm_norm <-  ref_se_gen(abm_c, pdata_abm)
readr::write_rds(sce_abm_norm, glue::glue("{data_dir}/ref_se_adult_bm_ERP122984.rds"))

# ======= cordblood data ========
cord <- readr::read_rds(glue::glue("{data_dir}/public/rds/cordblood_cells_matrix_sub.rds"))
cord_c <- t(cord)
pdata_cord <- readr::read_rds(glue::glue("{data_dir}/public/rds/cordblood_cells_metadata.rds"))
sce_cord_norm <-  ref_se_gen(cord_c, pdata_cord)
readr::write_rds(sce_cord_norm, glue::glue("{data_dir}/ref_se_cord_ERP122984.rds"))


