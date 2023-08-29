



library(phylogram)
library(dendextend)

tree_fn <- "/cluster/home/ylxie_jh/projects/leukemia/analysis/weinazhang/human/lyj_0321/inferCNV/M/result/infercnv.observations_dendrogram.txt"
#ic_dend <- read.dendrogram(tree_fn) ## doesn't work anymore because multiple trees
if_trees <- phytools::read.newick(tree_fn)
ic_dend_2 <- if_trees[[2]]
ic_dend_3 <- if_trees[[3]]
ic_dend_4 <- if_trees[[4]]

# Cut tree 
ic_labels_pre <- cutree(ic_dend_2, k = 7)
ic_labels_pre_II <- cutree(ic_dend_3, k = 2)

pre_cells <- data.frame(cell_id_new = names(ic_labels_pre), subcluster = ic_labels_pre)

sce_tumor <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/data/merged_4_tumors_new.rds")

tm_meta <- sce_tumor@meta.data
tm_meta$cell_id <- rownames(tm_meta)
tm_meta$cell_id_new <- sub(":", "-", tm_meta$cell_id)
pre_cells <- left_join(pre_cells, tm_meta[, c("cell_id_new", "cell_types_all", "cell_types_broad")])
pre_cells$sample <- substr(pre_cells$cell_id, 1, 2)
table(pre_cells[, c("subcluster", "cell_types_all")])
table(pre_cells[, c("subcluster", "sample")])

pre2_cells <- data.frame(cell_id_new = names(ic_labels_pre_II), subcluster = ic_labels_pre_II)
pre2_cells <- left_join(pre2_cells, tm_meta[, c("cell_id_new", "cell_types_all", "cell_types_broad")])
pre2_cells$sample <- substr(pre2_cells$cell_id, 1, 2)
table(pre2_cells[, c("subcluster", "cell_types_all")])

what <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/data/singleR_anno_fbm.rds")
what$cell_id_new <- sub(":", "-", rownames(what))
what <- what[, 3:5]
what <- as.data.frame(what)


pre_cells <- left_join(pre_cells, what)
pre_cells$sample <- substr(pre_cells$cell_id, 1, 2)
temp$cell_id_new <- sub(":", "-", rownames(temp))
pre_cells <- left_join(pre_cells, temp[, c("cell_id_new", "cell_types_broad", "cell_types_all")])




table(ic_labels)
# Color labels
the_bars <- as.data.frame(tableau_color_pal("Tableau 20")(20)[ic_labels])
colnames(the_bars) <- "inferCNV_tree"
the_bars$inferCNV_tree <- as.character(the_bars$inferCNV_tree)

ic_dend %>% set("labels",rep("", nobs(ic_dend)) ) %>% plot(main = "inferCNV dendrogram") %>%
  colored_bars(colors = as.data.frame(the_bars), dend = ic_dend,
               sort_by_labels_order = FALSE, add = T, y_scale = 100, y_shift = 0)

# compare with the main plot: 3 and 4 has strong patterns


ocells <- names(ic_labels)[ic_labels %in% 3:4]
ocells_origin <- substr(ocells, 1, 2)

xcells <- names(ic_labels)[ic_labels %in% 5:6]
xcells_origin <- substr(xcells, 1, 2)

readr::write_rds(ocells, "/cluster/home/yjliu_jh/projects/mef2d/output/cells_deletion_pattern_infercnv.rds")
#ocells <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/output/cells_deletion_pattern_infercnv.rds")
mef2d_meta <- sce@meta.data
rownames(mef2d_meta) <- sub(":", "-", rownames(mef2d_meta))
meta_sub_cnv <- mef2d_meta[rownames(mef2d_meta) %in% ocells, ]

sce_single_ball <- readr::read_rds("/cluster/home/ylxie_jh/projects/leukemia/analysis/weinazhang/human/split_hd/sce_single_ball_ctype_new.rds")
mef2d_meta_new <- sce_single_ball@meta.data
rownames(mef2d_meta_new) <- sub(":", "-", rownames(mef2d_meta_new))
meta_sub_cnv2 <- mef2d_meta_new[rownames(mef2d_meta_new) %in% ocells, ]
table(meta_sub_cnv2$cell_types_new)

















# new anno




Idents(merge) <- "cell_types_broad"
DotPlot(object = merge, features = c(gene03, "MME", "IGLL1", "ATXN1", "KIT", "CHCHD10", "HDAC9", "MEF2C", "MEF2D")) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
Idents(merge) <- "orig.ident"


Idents(sce_hd) <- "cell_types"
DotPlot(object = sce_hd, features = c(gene03, "MME", "IGLL1", "ATXN1", "KIT", "CHCHD10")) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
Idents(sce_hd) <- "orig.ident"





Idents(merge_sub) <- "temp_anno"
DotPlot(object = merge_sub, features = gene05) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
Idents(merge_sub) <- "orig.ident"



new_anno <- readr::read_csv("/cluster/home/ylxie_jh/projects/leukemia/analysis/weinazhang/human/anno_cor/anno_all_xyl.csv")
new_anno$broad_xyl <- new_anno$cell_types_xyl
new_anno$broad_xyl[grepl("preB_II", new_anno$broad_xyl)] <- "pre_B_II"
new_anno$broad_xyl[grepl("preB", new_anno$broad_xyl)] <- "preB_I"
new_anno$broad_xyl[grepl("proB", new_anno$broad_xyl)] <- "proB"
new_anno$broad_xyl[new_anno$cell_types_all %in% "preproB-like"] <- "pre_proB"
new_anno$broad_xyl[new_anno$broad_xyl %in% c("pre_B_II")] <- "preB_II"
new_anno$broad_xyl[new_anno$broad_xyl %in% c("NK", "T")] <- "NK_T"
new_anno$broad_xyl[new_anno$broad_xyl %in% c("MLP")] <- "CLP"

#table(new_anno$broad_xyl)
#identical(colnames(sce_tumor), new_anno$barcode)
sce_tumor@meta.data$cell_types_all <- new_anno$cell_types_xyl
sce_tumor@meta.data$cell_types_broad <- new_anno$broad_xyl
sce_tumor@meta.data$cell_types_new <- NULL
sce_tumor@meta.data$cell_types_broad_new <- NULL
sce_tumor@meta.data$cell_types <- sce_tumor@meta.data$cell_types_broad
sce_tumor@meta.data$cell_types_all <- NULL
sce_tumor@meta.data$cell_types_broad <- NULL


readr::write_rds(sce_tumor, "/cluster/home/yjliu_jh/projects/mef2d/data/sce_tumor_0329.rds")







Idents(sce_tumor) <- "cell_types"
sce_tumor@active.ident <- factor(sce_tumor@active.ident, 
                                 levels = c("HSC_MPP", "CLP", "pre_proB", "proB", "preB_I", 
                                            "preB_II", "immatureB", "matureB", "myeloid", "NK_T",
                                            "erythroid_cell"))
DotPlot(object = sce_tumor, features = gened2) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  scale_color_distiller(palette = "RdYlBu")
Idents(sce_tumor) <- "orig.ident"


sce_hd <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/analysis/zwn/human/tenx/seurat/merge/sce_hd.rds")

Idents(sce_hd) <- "cell_types"
sce_hd@active.ident <- factor(sce_hd@active.ident, 
                                 levels = c("HSC_MPP", "CLP", "pre_proB", "proB", "preB_I", 
                                            "preB_II", "immatureB", "matureB", "myeloid", "NK_T",
                                            "erythroid_cell"))
DotPlot(object = sce_hd, assay = "SCT", features = gened2[1:5]) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  scale_color_distiller(palette = "RdYlBu")
Idents(sce_hd) <- "orig.ident"







sce_hd_type <- readr::read_csv("/cluster/home/ylxie_jh/projects/leukemia/analysis/weinazhang/human/lyj_0321/projection_hd_old/sce_hd_ctype.csv")
sce_hd <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/data/healthy_donors.rds")
sce_hd_fil <- sce_hd[, sce_hd_type$barcode]


sce_tumor_old <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/data/merged_4_tumors_new.rds")
sce_tumor_old@meta.data$cell_types_broad[sce_tumor_old$cell_types_all %in% "preproB-like"] <- "pre_proB"

Idents(sce_tumor_old) <- "cell_types_all"

DotPlot(object = sce_tumor_old, assay = "RNA", features = c(gened2, "PAX5", "VPREB1")) + 
theme(axis.text.x = element_text(angle = 45, hjust=1)) +
scale_color_distiller(palette = "RdYlBu")


Idents(sce_tumor_old) <- "cell_types_broad"
DotPlot(object = sce_tumor, assay = "RNA", features = c(gened2, "PAX5", "VPREB1")) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  scale_color_distiller(palette = "RdYlBu")


FeaturePlot(sce_tumor, features = c("PAX5", "VPREB1", "DNTT", "SRGN", "CD72"))




sce_tumor_old@active.ident <- factor(sce_tumor_old@active.ident, 
                                     levels = c("HSC_MPP", "CLP", "pre_proB", "proB", "preB_I", 
                                                "preB_II", "immatureB", "matureB", "myeloid", "NK_T",
                                                "erythroid_cell"))


meta0 <- sce_tumor@meta.data
bar_data <- meta0[, c("orig.ident", "cell_types")]
colnames(bar_data)[1] <- "sample_id"
bar_data <- bar_data %>% group_by(sample_id) %>% mutate(ct = n()) %>%
  group_by(sample_id, cell_types) %>% mutate(ct2 = n()) %>% reframe(percentage = ct2 / ct) %>% unique()
bar_data$cell_types <- factor(bar_data$cell_types, 
                                    levels = c("HSC_MPP", "CLP", "pre_proB", "proB", "preB_I", 
                                               "preB_II", "immatureB", "matureB", "myeloid", "NK_T",
                                               "erythroid_cell"))


gene_dotplot <- c("MEIS1", "MLLT3", "ATP8B4", "CD99", "DNTT", "IGLL1", "VPREB1", "MME",
                  "LEF1", "EBF1", "CD24", "STMN1", "CD19", "HLA-DRA", "CD79A", "CD74", "IGHM", 
                  "IGKC", "MS4A1", "CST3", "SPI1", "CD68", "S100A9", "LTB", "NKG7", 
                  "CD3G", "CD3D", "HBA1", "HBA2", "CA1", "ALAS2", "SPINK2", "CD83", "SRGN")

gened2 <- c(gene_dotplot, setdiff(gene02, c(gene03, "MME", "IGLL1", "ATXN1", "KIT", "CHCHD10", "HDAC9", "MEF2C", "MEF2D")))










cell_color <- c("HSC_MPP" = "indianred1", "CLP" = "indianred3", "pre_proB" = "darkseagreen1",
                "proB" = "darkolivegreen1", "preB" = "lightgreen",  "preB_I" = "yellowgreen",
                "preB_II" = "#3CB371", "immatureB" = "chartreuse4", "matureB" = "darkgreen",
                "erythroid_cell" = "lightpink2", "NK_T" = "#9575aa", "myeloid" = "deeppink4")

# 1a  temp, just use a standard ggplot currently
ggplot(bar_data, aes(fill = cell_types, y = percentage, x = sample_id)) + 
  geom_bar(position = "fill", stat = "identity") + scale_fill_manual(values = cell_color)



readr::write_rds(new_anno, "/cluster/home/yjliu_jh/projects/mef2d/data/")

























tree_fn <- "/cluster/home/ylxie_jh/projects/leukemia/analysis/weinazhang/human/res0329/inferCNV/E/result/infercnv.observations_dendrogram.txt"
if_trees <- phytools::read.newick(tree_fn)   
ic_dend_5 <- if_trees[[5]]

ic_dend_4 <- if_trees[[4]]
mature_cells <- dendextend::cutree(ic_dend_4, k = 2)
mature_cells <- data.frame(cell_id_new = names(mature_cells), subcluster = mature_cells)
mature_cells <- left_join(mature_cells, tm_meta[, c("cell_id_new", "cell_id",  "cell_types")])

# Cut tree 
ic_labels_pro <- dendextend::cutree(ic_dend_5, k = 2)

pro_cells <- data.frame(cell_id_new = names(ic_labels_pro), subcluster = ic_labels_pro)

sce_tumor <- readr::read_rds("/cluster/home/ylxie_jh/projects/leukemia/analysis/weinazhang/human/tenx/seurat/merge/merge.rds")

tm_meta <- sce_tumor@meta.data
tm_meta$cell_id <- rownames(tm_meta)
tm_meta$cell_id_new <- sub(":", "-", tm_meta$cell_id)
pro_cells <- left_join(pro_cells, tm_meta[, c("cell_id_new", "cell_id",  "cell_types")])


procells2 <- pro_cells$cell_id[pro_cells$subcluster == 2]





pro_cells$sample <- substr(pro_cells$cell_id, 1, 2)
table(pre_cells[, c("subcluster", "cell_types_all")])
table(pre_cells[, c("subcluster", "sample")])











