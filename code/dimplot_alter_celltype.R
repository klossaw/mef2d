library("Seurat")
se <- read_rds("/cluster/home/jhuang/projects/mef2d/analysis/zwn/human/tenx/seurat/merged/mef2d.rds")
#cell_id <- colnames(se)
#new_cell_id <- str_sub(cell_id, 1 + str_locate(cell_id, ":")[, 1] / 2)
#se <- RenameCells(se, new.names = new_cell_id)

# set colors
all_cell_color <- c("HSC_MPP" = "indianred1", "CLP" = "indianred3", "pre_proB" = "darkseagreen1",
                    "proB" = "lightgreen", "preB" = "darkolivegreen1",  "preB_I" = "yellowgreen",
                    "preB_II" = "chartreuse4", "imatureB" = "#3CB371", "matureB" = "darkgreen",
                    "erythroid_cell" = "lightpink2", "NK_T" = "#9575aa", "myeloid" = "deeppink4", "tumor1" = "grey",
                    "tumor_cell1" = "yellow4", "tumor_cell2" = "aquamarine2", "tumor_cell3" = "lightseagreen")

# # first plot
# psub1 <- DimPlot(se, group.by = "cell_types", reduction='umap', label = T) + NoLegend() + 
#   scale_color_manual(values = all_cell_color)
# 
# psub2 <- DimPlot(se, group.by = "orig.ident", reduction='umap', label = T) + NoLegend()
# 
# psub_all <- psub1 + psub2 + plot_layout(nrow = 1)
# 
# #ggsave(psub_all, filename = "/cluster/home/yjliu_jh/projects/temp/111_umap_allsamplestest.pdf", width = 9, height = 4)
# 
# 
# # second plot
# psub1test <- DimPlot(se, group.by = "SCT_snn_res.0.6", reduction='umap', label = T) + NoLegend() 
# 
# psub2test <- DimPlot(se, group.by = "seurat_clusters", reduction='umap', label = T) + NoLegend() 
# 
# psub_alltest <- psub1test + psub2test + plot_layout(nrow = 1)
# 
# #ggsave(psub_alltest, filename = "/cluster/home/yjliu_jh/projects/temp/111_umap_allsamplestest.pdf", width = 9, height = 4)
# 
# 
# # match 3 types of cells to seurat_clusters
# # (can be done by code, but quicker by eyes)
# 
# # myeloid: 21 39     NK-T: 7 16      erythroid_cell  0 22 25 28 33 40 41
# 
se_meta <- se@meta.data
# se_meta$cell_types[se_meta$seurat_clusters %in% c(21, 39)] <- "myeloid"
# se_meta$cell_types[se_meta$seurat_clusters %in% c(7, 16)] <- "NK_T"
# se_meta$cell_types[se_meta$seurat_clusters %in% c(0, 22, 25, 28, 33, 40, 41)] <- "erythroid_cell"
# se@meta.data <- se_meta
# #ggsave(psub_all, filename = "/cluster/home/yjliu_jh/projects/temp/111_umap_allsamples.pdf", width = 9, height = 4)
# 
# 


# set regions based on observation and existing data
umap_loc <- Embeddings(se, "umap")
# observation:
in_ery_ind <- umap_loc[, 1] < -5.8 & umap_loc[, 2] < 5
se_meta$cell_types[in_ery_ind] <- "erythroid_cell"
se_meta$cell_types[!in_ery_ind & se_meta$cell_types %in% "erythroid_cell"] <- NA

in_mye_ind <- umap_loc[, 1] > -5 & umap_loc[, 2] > 13.9
se_meta$cell_types[in_mye_ind] <- "myeloid"
se_meta$cell_types[!in_mye_ind & se_meta$cell_types %in% "myeloid"] <- NA

in_nkt_ind <- umap_loc[, 1] > -2.955 & umap_loc[, 1] < 4.5 &
              umap_loc[, 2] > 9 & umap_loc[, 2] < 13.9
se_meta$cell_types[in_nkt_ind] <- "NK_T"
se_meta$cell_types[!in_nkt_ind & se_meta$cell_types %in% "NK_T"] <- NA
in_nkt_ind <- umap_loc[, 1] > -1.8 & umap_loc[, 1] < 0 &
              umap_loc[, 2] > 8.4 & umap_loc[, 2] < 13.9
se_meta$cell_types[in_nkt_ind] <- "NK_T"


se@meta.data <- se_meta

# plot
psub1 <- DimPlot(se, group.by = "cell_types", reduction='umap', label = T) + NoLegend() + 
  scale_color_manual(values = all_cell_color)
ggsave(psub1, filename = "/cluster/home/yjliu_jh/projects/temp/111_umap_all_altered.pdf", width = 5, height = 4)



readr::write_rds(se, "/cluster/home/jhuang/projects/mef2d/analysis/zwn/human/tenx/seurat/merged/mef2d.rds")




# data: not working!!!

# mye_sub_loc <- umap_loc[se_meta$cell_types %in% "myeloid", ]
# in_mye_ind <- umap_loc[, 1] < quantile(mye_sub_loc[, 1], 0.99) &
#   umap_loc[, 1] > quantile(mye_sub_loc[, 1], 0.01) &
#   umap_loc[, 2] < quantile(mye_sub_loc[, 2], 0.99) & 
#   umap_loc[, 2] > quantile(mye_sub_loc[, 2], 0.01)
# se_meta$cell_types[in_mye_ind] <- "myeloid"









