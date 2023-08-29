pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr",
          "loomR", "SummarizedExperiment", "Matrix"
)
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}


# select best n for constructing ref matrix 

k_select_fc <- function(mat, mat_anno, cluster_col = "cluster", cell_group_col){
  i <- 1
  mat_con_num <- 1000000
  # select n from 20 to 150, calculation may take a while
  for (n in c(seq(20, 60, by = 2), seq(70, 150, by = 5))){
    mat_top <- mat %>% arrange(desc(absFC)) %>% group_by(.[[cluster_col]]) %>% dplyr::slice(1:n)
    mat_s <- t(mat_assay[unique(mat_top$gene), ])
    mat_s0 <- aggregate(mat_s, by = list(mat_anno[, id_col]), mean)
    mat0 <- as.matrix(mat0[, -1])
    rownames(mat0) <- mat_s0$Group.1
    i <- i + 1
    mat_con_num[i] <- kappa(mat0)
    if (diff(mat_con_num)[i - 1] > 0) { break }
  }
  return(list(ref_mat = mat0, ref_assay = mat_s))
}


# ====== construct reference matrices for adult bone marrow data =======
abm_coldata <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/data/public/rds/adultbm_cells_metadata.rds")
abm_assay <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/data/public/rds/adultbm_cells_matrix_sub.rds")
all_markers_abm <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/data/public/rds/all_markers_abm.rds")
all_markers_abm$absFC <- abs(all_markers_abm$avg_log2FC)
am_abm_sub <- all_markers_abm[all_markers_abm$p_val_adj < 0.05, ]
am_abm_sub <- am_abm_sub[am_abm_sub$gene %in% colnames(abm_assay), ]
abm_exp_perc <- colSums(abm_assay == 0)
exp_genes_abm <- names(abm_exp_perc[abm_exp_perc < (0.8 * nrow(abm_assay))])
am_abm_sub2 <- am_abm_sub[am_abm_sub$gene %in% exp_genes_abm, ]
table(am_abm_sub2$cluster)


res_abm <- k_select_fc(am_abm_sub2, abm_coldata, "cluster", "cell_identity")
abm_sub <- res_abm[["ref_assay"]]
abm_sub <- abm_sub[-332, ]  ## remove duplicated MALAT1 (due to source data?)
readr::write_rds(res_abm[["ref_mat"]], "/cluster/home/yjliu_jh/projects/mef2d/data/public/rds/ref_mat_abm.rds")
readr::write_rds(abm_sub, "/cluster/home/yjliu_jh/projects/mef2d/data/public/rds/abm_assay_sub.rds")


# ==== construct reference matrices for fetal bone marrow data =======


# ==== construct reference matrices for cordblood data =======


# ==== construct reference matrices for aml data from PMID:30827681 =======
# merits all the same as above, fill up later


