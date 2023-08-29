






# ==== find markers for leukemia dataset from fbm =======

ref_mat <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/data/h5ad1_mat_dgc.rds")
ref_anno <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/data/h5ad1_obs.rds")
ref_anno <- ref_anno[rownames(ref_mat), ]
ref_mat2 <- t(ref_mat)
leu_ref <- CreateSeuratObject(counts = ref_mat2, project = "fbm_anno", meta.data = ref_anno,
                              min.cells = 3, min.features = 200)
Idents(leu_ref) <- leu_ref@meta.data$cell.labels
leu_ref <- NormalizeData(leu_ref)
future::plan("multicore", workers = 20)
all_markers_leu <- FindAllMarkers(leu_ref)
future::plan("sequential")
readr::write_rds(all_markers_leu, "/cluster/home/yjliu_jh/projects/mef2d/data/public/rds/all_markers_fbm.rds")




# ====== construct reference matrices for fetal bone marrow data =======
fbm_coldata <- ref_anno  
fbm_assay <- ref_mat  ## row:cell x column:gene

future::plan("multicore", workers = 20)
all_markers_fbm <- FindAllMarkers(leu_ref)
future::plan("sequential")
readr::write_rds(all_markers_fbm, "/cluster/home/yjliu_jh/projects/mef2d/data/public/rds/all_markers_fbm.rds")
all_markers_fbm$absFC <- abs(all_markers_fbm$avg_log2FC)
am_fbm_sub <- all_markers_fbm[all_markers_fbm$p_val_adj < 0.05, ]
am_fbm_sub <- am_fbm_sub[am_fbm_sub$gene %in% colnames(fbm_assay), ]
fbm_exp_perc <- colSums(fbm_assay == 0)
exp_genes_fbm <- names(fbm_exp_perc[fbm_exp_perc < (0.8 * nrow(fbm_assay))])
am_fbm_sub2 <- am_fbm_sub[am_fbm_sub$gene %in% exp_genes_fbm, ]
table(am_fbm_sub2$cluster)
identical()




k_select_fc <- function(mat, mat_assay, mat_anno, cluster_col = "cluster", cell_group_col){
  i <- 1
  mat_con_num <- 1000000
  # select n from 20 to 150, calculation may take a while
  for (n in c(seq(20, 60, by = 2), seq(70, 150, by = 5))){
    mat_top <- mat %>% arrange(desc(absFC)) %>% group_by(.[[cluster_col]]) %>% dplyr::slice(1:n)
    mat_s <- mat_assay[, unique(mat_top$gene)]
    mat_s0 <- aggregate(mat_s, by = list(mat_anno[, cell_group_col]), mean)
    mat0 <- as.matrix(mat_s0[, -1])
    rownames(mat0) <- mat_s0$Group.1
    i <- i + 1
    mat_con_num[i] <- kappa(mat0)
    if (diff(mat_con_num)[i - 1] > 0) { break }
  }
  return(list(ref_mat = mat0, ref_assay = mat_s))
}




res_fbm <- k_select_fc(am_fbm_sub2, fbm_assay, fbm_coldata, "cluster", "cell.labels")
readr::write_rds(res_fbm, "/cluster/home/yjliu_jh/projects/mef2d/data/public/rds/res_fbm_temp.rds")







