# pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", "optparse",
#           "vroom", "jhtools", "glue", "openxlsx", "ggsci", "patchwork", "cowplot",
#           "tidyverse", "dplyr", "Seurat", "ggridges", "viridis", "circlize", "ComplexHeatmap",
#           "matrixStats")
# for (pkg in pkgs){
#   suppressPackageStartupMessages(library(pkg, character.only = T))
# }
# 
# hcc_samples <- c("hcc5T", "hcc3T", "hcc4T", "hcc2T", "hcc1T", "hcc4C", "hcc1C", "hcc3C")
# hcc_samples_two <- paste0(hcc_samples, "_two")
# hcc_tumor <- hcc_samples[1:5]
# hcc_normal <- hcc_samples[6:8]
# 
# 
# 
# for (i in 1:length(hcc_samples)) {
#   temp <- readr::read_rds(glue::glue("/cluster/home/ztao_jh/projects/liver/analysis/qzhang/human/meta/checkfilter/{hcc_samples[i]}_obj_res2_filter2_two_matrix.rds"))
#   assign(hcc_samples_two[i], temp)
# }
# 
# 
# 
# 
# top_10_dist <- function(ind){
#   stopifnot("make sure input are two-column indices" = (dim(ind)[2] == 2 ))
#   dist <- RANN::nn2(ind, ind, k = 10)$nn.dists
#   dist2 <- mean(rowMeans(dist))
#   dist2
# }
# 
# top_10_in_20_dist <- function(ind){
#   stopifnot("make sure input are two-column indices" = (dim(ind)[2] == 2 ))
#   dist <- RANN::nn2(ind, ind, k = 20)$nn.dists
#   dist_order <- order(dist)
#   dist2 <- mean(sort(rowMeans(dist))[1:10])
#   dist2
# }
# 
# 

# 
# for (j in 1:length(hcc_samples)) {
#   hcc_temp <- get(hcc_samples_two[j])
#   spexp <- as.matrix(hcc_temp@assays$SCT@data)
#   spexp <- spexp[rowQuantiles(spexp, probs = (ncol(spexp) - 30) / ncol(spexp)) > 0, ]
#   exp_cl <- as.data.frame(matrix(0, nrow(spexp), 4))
#   colnames(exp_cl) <- c("gene", "dist10", "dist20", "count")
#   for (i in 1:nrow(spexp)){
#     inten <- spexp[i, ]
#     coor <- hcc_temp@images$slice1@coordinates[, 2:3]
#     ind10 <- coor[names(inten[order(inten, decreasing = T)])[1:10], ]
#     ind20 <- coor[names(inten[order(inten, decreasing = T)])[1:20], ]
#     exp_cl[i, 1] <- rownames(spexp)[i]
#     exp_cl[i, 2] <- top_10_dist(ind10)
#     exp_cl[i, 3] <- top_10_in_20_dist(ind20)
    exp_cl[i, 4] <- sum(inten > 0.5 * max(inten))
#   }
#   readr::write_rds(exp_cl, glue::glue("{out_dir}/gene_cl_mat_{hcc_samples[j]}_new.rds"))
# }

# temp <- list()
# for (i in 1:length(hcc_samples)) {
#   temp[[hcc_samples[i]]] <- readr::read_rds(glue::glue("{out_dir}/gene_cl_mat_{hcc_samples[i]}_new.rds"))
# }
# gene_cl_new <- bind_rows(temp, .id = "samples")
# 
# gene_cl_new <- gene_cl_new[order(gene_cl_new$dist10), ]
# 
# 



# cl_gene_top50000 <- readr::read_rds("/cluster/home/yjliu_jh/projects/hcc/output/cl_gene_50k.rds")
# cl_gene_top50000$m <- cl_gene_top50000$diff * cl_gene_top50000$cons 
# cl_gene_top <- cl_gene_top50000[cl_gene_top50000$m > 0.01, ]
# 
# 


# s_samples <- unique(cl_gene_top$samples)
# # get filtered metabolite list and use to calculate gene overlap
# for (i in 2:length(s_samples)){
#   res_list <- list()
#   gene_part <- cl_gene_top[cl_gene_top$samples %in% s_samples[i], ]
#   genes <- gene_part$gene
   counts <- gene_part$count
#   hcc_temp <- get(paste0(s_samples[i], "_two"))
#   spmeta <- hcc_temp@assays$metoblism@data
#   spexp <- as.matrix(hcc_temp@assays$SCT@data)
#   spexp <- spexp[, intersect(colnames(spmeta), colnames(spexp))]
#   spexp <- spexp[genes, ]
#   spot_coordinates <- hcc_temp@images$slice1@coordinates[, 2:3]
#   for(m in 1:length(genes)){
#     if (length(genes) == 1) {inten <- spexp
#     } else {
#       inten <- spexp[m, ]}
    coor_gene <- spot_coordinates[names(sort(inten, decreasing = T)[1:counts[m]]), ]
#     spmeta_fil <- spmeta[rowQuantiles(spmeta, probs = 1 - (counts[m] / length(inten))) > 0, ]
#     metabolites <- rownames(spmeta_fil)
#     metab_spots_mat <- matrix(0, counts[m], length(metabolites))
#     colnames(metab_spots_mat) <- metabolites
#     metab_spots_mat <- as.data.frame(metab_spots_mat)
#     for (n in 1:length(metabolites)) {
#       metab_spots_mat[, n] <- order(spmeta_fil[metabolites[n], ], decreasing = T)[1:counts[m]]
#     }
#     res_list[[genes[m]]] <- sapply(metab_spots_mat, 
#                                    \(x) {group_distance6(spot_coordinates[x, ], coor_gene)})
#   }
#   readr::write_rds(res_list, glue::glue("{out_dir}/gene_based_dist_{s_samples[i]}_1.rds"))
#   res_list <- list()
# }











