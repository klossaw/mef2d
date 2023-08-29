pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr", "igraph",
          "Seurat", "anndata", "stats", "pheatmap", "pracma", "SeuratDisk",
          "phylogram", "dendextend", "clusterProfiler", "org.Hs.eg.db"
)
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}

mdata <- anndata::read_h5ad('/cluster/home/yjliu_jh/projects/mef2d/data/metacells.h5ad')

mdata_mat <- as(as.matrix(mdata$X), "dgCMatrix")
umis <- as.matrix(mdata$X)
fractions <- umis / rowSums(umis)
log_fractions <- log2(1e-5 + fractions)

genes_per_metacell <- 2
minimal_max_log_fraction_of_interesting_genes <- -10
minimal_relative_log_fraction_of_candidate_genes <- 2

median_log_fractions_of_genes <- apply(log_fractions, 2, median)
relative_log_fractions <- sweep(log_fractions, 2, median_log_fractions_of_genes)
max_log_fractions_of_genes <- apply(log_fractions, 2, max)


genes_hspc <- c("KIT", "IL7R", "ATXN1", "FLT3")
genes_prob <- c("IL7R", "CD45", "CD43", "CD19", "RAG1", "RAG2", "EBF1", "PAX5", "SOX4", "LEF1")
genes_preb <- c("IL7R", "CD45", "CD43", "CD19", "RAG1", "RAG2", "SOX4", "LEF1", "EIF4EBP1")
genes_imb <- c("CD19", "EBF1", "CD24", "EIF4EBP1", "TNFRSF13C")



"VPREB1"





sce_new1 <- readr::read_rds("/cluster/home/ylxie_jh/projects/leukemia/analysis/weinazhang/human/split_hd/sce_single_ball_ctypenew1.rds")
sce_new1@meta.data$cell_id_new <- sub(":", "-", rownames(sce_new1@meta.data))
sce_new1@meta.data <- left_join(sce_new1@meta.data, ctomc[, c("cell_id", "cell_id_new", "cluster", "cn_cluster")])

table(sce_new1@meta.data$SCT_snn_res.0.6[sce_new1@meta.data$cell_id_new %in% ocells])





# ==== others    map metacell clusters to current umap =====


# read seurat
sce <- readr::read_rds("~/projects/mef2d/analysis/zwn/human/tenx/seurat/merged/mef2d.rds")
sce@meta.data$group <- ifelse(sce@meta.data$orig.ident %in% c("B069_Dx_CD19", "M2", "M3", "M4"), "MEF2D", "others")
Idents(object = sce) <- "group"

# some random differential tests
sce_tumor <- subset(sce, `orig.ident` %in% c("B069_Dx_CD19", "M2", "M3", "M4", "E1"))
sce_preanno <- readr::read_rds("/cluster/home/ylxie_jh/projects/leukemia/analysis/weinazhang/human/public/final/mef2d_tumor_preanno.rds")
# seems this file contains all cells except erythoid, myeloid and nk_t
temp_join <- data.frame(cell_id = rownames(sce_preanno@meta.data), pre_anno = sce_preanno@meta.data$pre_anno)
temp <- sce_tumor@meta.data
temp <- left_join(temp, temp_join)
temp <- temp %>% mutate(pre_anno2 = coalesce(pre_anno, cell_types))
cnv_scores <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/data/cnvscores.rds")
temp <- left_join(temp, cnv_scores)


sce_tumor@meta.data$pre_anno <- temp$pre_anno
sce_tumor@meta.data$pre_anno2 <- temp$pre_anno2
sce_tumor@meta.data$cnv1 <- sqrt(temp$cnv1)

tm_meta <- sce_tumor@meta.data
tm_meta$cell_id_new <- sub(":", "-", rownames(tm_meta))
tm_meta$cn <- ifelse(tm_meta$cell_id_new %in% ocells, "altered", "not_altered")
#table(tm_meta$seurat_clusters[tm_meta$cell_id_new %in% ocells])
sce_tumor@meta.data$cn <- tm_meta$cn
sce_tumor@meta.data$cell_id <- rownames(sce_tumor@meta.data)

a3 <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/output/mc_tumor_clust.rds")
sce_tumor@meta.data <- left_join(sce_tumor@meta.data, a3)
rownames(sce_tumor@meta.data) <- sce_tumor@meta.data$cell_id


sce_tumor <- FindNeighbors(sce_tumor, dims = 1:30, verbose = F) 
sce_tumor <- FindClusters(sce_tumor, resolution = 0.5, verbose = F)
sce_tumor <- RunUMAP(sce_tumor, dims = 1:30, verbose = F)


sce_tumor <- RunUMAP(sce_tumor, dims = 1:30, verbose = F, min.dist = 0.5, spread = 1)


sce <- FindNeighbors(sce, dims = 1:30, verbose = F) 
sce <- FindClusters(sce, resolution = 0.5, verbose = F)
sce <- RunUMAP(sce, dims = 1:30, verbose = F, min.dist = 0.05, spread = 4)
sce@meta.data$cell_id_new <- sub(":", "-", rownames(sce@meta.data))
sce@meta.data$cn <- ifelse(sce@meta.data$cell_id_new %in% ocells, "altered", "not_altered")



table(sce_tumor@meta.data$SCT_snn_res.1)
sce_tumor@meta.data$cell_id_new <- sub(":", "-", rownames(sce_tumor@meta.data))
table(sce_tumor@meta.data$SCT_snn_res.2[sce_tumor@meta.data$cell_id_new %in% ocells])
sce_tumor@meta.data <- left_join(sce_tumor@meta.data, ctomc[, c("cell_id", "cell_id_new", "cluster", "cn_cluster")])
sce_tumor@meta.data$cn <- ifelse(sce_tumor@meta.data$cell_id_new %in% ocells, "altered", "not_altered")
rownames(sce_tumor@meta.data) <- names(sce_tumor$cell_id)

# then visualization

facet(DimPlot(sce_tumor, group.by = 'clust', reduction = 'umap'), facet.by = "clust")
facet(DimPlot(sce_tumor, group.by = 'pre_anno2', reduction = 'umap') , facet.by = "pre_anno2")
facet(DimPlot(sce_tumor, group.by = 'pre_anno', reduction = 'umap') , facet.by = "pre_anno")
facet(DimPlot(sce_tumor, group.by = 'seurat_clusters', reduction = 'umap') , facet.by = "seurat_clusters")
facet(DimPlot(sce_tumor_fil, group.by = 'cluster_anno', reduction = 'umap') , facet.by = "cluster_anno")
DimPlot(sce_tumor2, group.by = 'cn_cluster', reduction = 'umap') 




FeaturePlot(sce_tumor, features = c("cnv1", "MKI67", "ATXN1", "DNTT", "CD19", "MS4A1", "CD22",
                                    "KIT", "FLT3", "CD34", "EIF4EBP1", "IGKC", "IGLC2"))

plotht <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/output/plotht.rds")



pht5data <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/output/pht5data.rds")
gene01 <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/output/gene01.rds")

Idents(sce_tumor_fil) <- "cluster_anno"
DotPlot(object = sce_tumor_fil, features = gene01 #, split.by = 'cluster_final_anno',  cols = fa_colors
        ) + theme(axis.text.x = element_text(angle = 45, hjust=1))

mcn_order <- pht5_data$tree_col$labels[pht5$tree_col$order]

DotPlot(object = sce_tumor_fil, features = gene02) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))

gene02 <- c("DNTT", "CD24", "VPREB1", "CD38", "LEF1", "FLT3", "BLK",
            "EBF1", "PAX5", "IGHM", "TNFRSF13C", "IGKC", "MS4A1", 
            "CD72", "CD19", "CD79A", "HLA-DRA", "CD74", "SPI1", "CD68",
            "S100A9", "LILRB2", "MEIS1", "CST3", "LTB", "CHCHD10",
            "PDE4B", "IGHG1", "NKG7", "KLRB1", "CD27", "CD3D", "CD3G", 
            "MLLT3", "CD2", "LAT", "SOX4", "CD22", "C1QTNF4", "SMIM24", 
            "BCL11B", "IL7R", "EIF4EBP1", "STMN1", "SPINK2", "HBA1", "HBA2",
            "MKI67", "CD83", "ATP8B4", "CD99", "IGHD", "ATXN1", "SRGN")



fa_colors <- c("#ee5d5a", "#34956c", "#6b6498", "#71a3a2", "#c88978", "#cc9b32", "#53096a", 
                              RColorBrewer::brewer.pal(7, "Set2"), RColorBrewer::brewer.pal(6, "Pastel1"))



sce_tumor2 <- subset(sce, `orig.ident` %in% c("B069_Dx_CD19", "M2", "M3", "M4", "E1"))
sce_tumor2@meta.data <- sce_tumor@meta.data
DimPlot(sce_tumor2, group.by = 'cluster_new', reduction = 'umap')


umap_pa <- as.data.frame(sce_preanno@reductions$umap@cell.embeddings)
umap_pa$cell_id <- rownames(umap_pa)
umap_pa <- left_join(umap_pa, ctomc[, c("cell_id", "cell_id_new", "cluster")])

ggscatter(umap_pa, "UMAP_1", "UMAP_2", color = "cluster", size = 0.2)


umap_tm <- as.data.frame(sce_tumor@reductions$umap@cell.embeddings)
umap_tm$cell_id <- rownames(umap_tm)
umap_tm <- left_join(umap_tm, ctomc[, c("cell_id", "cell_id_new", "cluster", "cn_cluster")])
umap_tm$cn_cluster <- as.character(umap_tm$cn_cluster)

facet(ggscatter(umap_tm, "UMAP_1", "UMAP_2", color = "cn_cluster", size = 0.2), facet.by = "cn_cluster")



umap_tm2 <- as.data.frame(sce_tumor@reductions$umap@cell.embeddings)
umap_tm2$cell_id <- rownames(umap_tm2)
umap_tm2 <- left_join(umap_tm2, ctomc[, c("cell_id", "cell_id_new", "cluster", "cn_cluster")])
umap_tm2$cn_cluster <- as.character(umap_tm2$cn_cluster)

#facet(ggscatter(umap_tm2, "UMAP_1", "UMAP_2", color = "cn_cluster", size = 0.2), facet.by = "cn_cluster")
# seems pretty valid if use meta-cell clusters   but the cn subgroup is separated 


umap_tm2$cn <- ifelse(umap_tm2$cell_id_new %in% ocells, "altered", "not_altered")
ggscatter(umap_tm2, "UMAP_1", "UMAP_2", color = "cn", size = 0.2)





# read metacell info
umis <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/output/metacells_umis.rds")
ctomc <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/output/cells_to_metacells.rds")


ocells <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/output/cells_deletion_pattern_infercnv.rds")




# 
ddr_test <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/output/ddrtree_object_testxx.rds")
ddr_test@phenoData@data$clust[is.na(ddr_test@phenoData@data$clust)] <- "NA"
ddr_test2 <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/output/ddrtree_object_1.rds")
ddr_test2@phenoData@data$cell_types2 <- ddr_test@phenoData@data$cluster_new
# ddr_test@phenoData@varMetadata

samples <- unique(ddr_test@phenoData@data$clust)
sample_color <- c("#ee5d5a", "#9f2f6a", "#34956c", "#6b6498","#b29bc9", "#71a3a2", 
                  "#c88978",  "#a15217", "#ce662f", "#cc9b32", "#53096a", "#6569c2", "#66676c",
                  RColorBrewer::brewer.pal(7, "Set2"))[1:length(samples)]
names(sample_color) <- samples


monocle::plot_cell_trajectory(ddr_test2, color_by = "cell_types2", alpha = 0.3, size = 1, show_backbone = TRUE) +
  theme(legend.title = element_blank()) + scale_color_discrete(sample_color) + facet_wrap("~cell_types2", nrow = 3)








# ===== cn =======

ocells <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/output/cells_deletion_pattern_infercnv.rds")
ic_dend <- read.dendrogram("/cluster/home/ylxie_jh/projects/leukemia/analysis/weinazhang/human/split_hd/infercnv/result/infercnv.observations_dendrogram.txt")
ic_labels <- cutree(ic_dend, k = 10)
sce_tumor_cn <- sce_tumor[, sce_tumor$cell_id_new %in% names(ic_labels)]
sce_tumor_cn_fil <- sce_tumor_cn[, sce_tumor_cn$cell_id %in% sce_tumor_fil$cell_id]
Idents(sce_tumor_cn_fil) <- "cn"

diff_markers_cnv <- FindMarkers(sce_tumor_cn_fil, ident.1 = "altered", ident.2 = "not_altered", method = "DESeq2")
diff_markers_cnv$perc_diff <- diff_markers_cnv$pct.1 - diff_markers_cnv$pct.2
diff_res <- diff_markers_cnv[order(diff_markers_cnv$perc_diff), ]  ## order: low expressed in ocells
diff_cn_genes <- rownames(diff_res)

check_ego_gene <- function(data){
  ego <- enrichGO(gene = data, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                  ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
  ego <- ego[order(ego$p.adjust), ]
  ego
}

check_ego_gene(head(diff_cn_genes, 50))  ## high in other cells
# mitochondrial respirasome, ATP synthesis, translation initiation... and percent differs significantly
check_ego_gene(tail(diff_cn_genes, 50))  
# 





# ====== markers from yali ======
MEP = c("FCER1A", "GATA1", "KLF1", "PBX1")
Stromal = c("MCAM", "NGFR", "VCAM1", "THY1")
CMP = c("CD34", "CD38", "FLT3")
MLP = c("C1QTNF4", "SPINK2", "SMIM24", "CSF3R")
NK_T = c("CD3D","CD3G","NKG7", "KLRC1")
T_B = c("MS4A1")
Early_ERP = c("GATA1", "CNRIP1", "SLC40A1", "PKIG")
Erythroid_cell = c("HBD", "CA1")
Macrophage = c("CD68")
Monocyte = c("CD14","S100A9","CST3")
MDP = c("MPO", "ELANE", "LYZ")
cDC = c("CD1C", "FCER1A")
pDC = c("IL3RA", "IRF8")
HPC = c("CD34", "CD38", "HMGB3", "SRGN")
HSC_MPP = c("CD34","AVP", "SPINK2", "MLLT3")
LMPP = c("CD34", "ACY3", "RUNX2")
CLP = c("CYGB", "LTB", "IL7R", "CD99")
#Pre_pro_B = c("DNTT","IGLL1")
Pro_B = c("DNTT","CD38")
Pre_B = c("IGLL1","MME")
Immature_B = c("CD19")
Mature_B = c("MS4A1", "CD79A")
Plasma_B = c("CD27")

genes_hspc <- c("KIT", "IL7R", "ATXN1", "FLT3")
genes_prepro_prob <- c("IL7R", "CD45", "CD43", "CD19", "RAG1", "RAG2", "EBF1", "PAX5", "SOX4", "LEF1", "DNTT", "VPREB1", "IGLL1")
genes_preb <- c("IL7R", "CD45", "CD43", "CD19", "RAG1", "RAG2", "SOX4", "LEF1", "EIF4EBP1")
genes_imb <- c("CD19", "EBF1", "CD24", "EIF4EBP1", "TNFRSF13C")
genes_matureb <- c("CD20", "CD79A", "CD19", "EIF4EBP1", "TNFRSF13C")
genes_plasma <- c("IGHM", "IGHD", "TNFRSF13C", "IGKC", "IGLC2", "CD20", "CD79A")
genes_other <- c("CD22", "CHCHD10", "CD74", "MKI67", "STMN1", "MS4A1")


add_genes <- unique(c(MEP, Stromal, CMP, MLP, NK_T, T_B,
                      Early_ERP, Erythroid_cell, Macrophage, Monocyte, MDP, cDC, pDC, HPC, HSC_MPP,
                      LMPP, CLP, Pro_B, Pre_B, Immature_B, Mature_B, Plasma_B, genes_hspc, genes_prepro_prob,
                      genes_preb, genes_imb, genes_matureb, genes_plasma, genes_other))


