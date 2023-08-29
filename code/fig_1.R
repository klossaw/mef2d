pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr", 
          "Seurat", "ComplexHeatmap")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}

# read in data
sce_tumor <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/data/merged_4_tumors.rds")
sce_hd <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/data/healthy_donors.rds")
meta0 <- sce_tumor@meta.data

# set colors
cell_color <- c("HSC_MPP" = "indianred1", "CLP" = "indianred3", "pre_proB" = "darkseagreen1",
                "proB" = "darkolivegreen1", "preB" = "lightgreen",  "preB_I" = "yellowgreen",
                "preB_II" = "#3CB371", "immatureB" = "chartreuse4", "matureB" = "darkgreen",
                "erythroid_cell" = "lightpink2", "NK_T" = "#9575aa", "myeloid" = "deeppink4")

# 1a data
bar_data <- meta0[, c("orig.ident", "cell_types")]
colnames(bar_data)[1] <- "sample_id"
bar_data <- bar_data %>% group_by(sample_id) %>% mutate(ct = n()) %>%
  group_by(sample_id, cell_types) %>% mutate(ct2 = n()) %>% reframe(percentage = ct2 / ct) %>% unique()
bar_data$cell_types <- factor(bar_data$cell_types, 
                              levels = c("HSC_MPP", "CLP", "pre_proB", "proB", "preB_I", 
                                         "preB_II", "immatureB", "matureB", "myeloid", "NK_T",
                                         "erythroid_cell"))

# 1a plot (temp, just use a standard ggplot currently)
pdf("/cluster/home/yjliu_jh/projects/mef2d/output/fig1atemp2.pdf", height = 4, width = 4)
ggplot(bar_data, aes(fill = cell_types, y = percentage, x = sample_id)) + 
  geom_bar(position = "fill", stat = "identity") + scale_fill_manual(values = cell_color) +
  cowplot::theme_cowplot()
dev.off()





# 1d dotplot ver.

gene_dotplot <- c("MEIS1", "MLLT3", "ATP8B4", "CD99", "DNTT", "IGLL1", "VPREB1", "MME",
                  "LEF1", "EBF1", "CD24", "STMN1", "CD19", "HLA-DRA", "CD79A", "CD74", "IGHM", 
                  "IGKC", "MS4A1", "CST3", "SPI1", "CD68", "S100A9", "LTB", "NKG7", 
                  "CD3G", "CD3D", "HBA1", "HBA2", "CA1", "ALAS2", "SPINK2", "CD83", "SRGN")

Idents(sce_hd) <- "cell_types"
sce_hd@active.ident <- factor(sce_hd@active.ident, 
                              levels = c("HSC_MPP", "CLP", "pre_proB", "proB", "preB_I", "preB_II", 
                                         "immatureB", "matureB", "myeloid", "NK_T", "erythroid_cell"))
DotPlot(object = sce_hd, assay = "SCT", features = gene_dotplot) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  scale_color_distiller(palette = "RdYlBu")
Idents(sce_hd) <- "orig.ident"




Idents(sce_tumor) <- "cell_types"
sce_tumor@active.ident <- factor(sce_tumor@active.ident, 
                              levels = c("CLP", "pre_proB", "proB", "preB_I", "preB_II", 
                                         "immatureB", "matureB", "myeloid", "NK_T", "erythroid_cell"))
DotPlot(object = sce_tumor, assay = "SCT", features = unique(c(gene_dotplot, gene02))) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  scale_color_distiller(palette = "RdYlBu")
Idents(sce_tumor) <- "orig.ident"








# 1b temp 
DimPlot(sce_tumor, group.by = 'orig.ident', reduction = 'umap')
DimPlot(sce_tumor, group.by = 'cell_types_broad', reduction = 'umap')


# 1c temp   may not so significant but just this now
Idents(sce_tumor) <- "cell_types_broad"
DotPlot(object = sce_tumor, features = gene03) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))

Idents(sce_tumor) <- "cell_types_all"
DotPlot(object = sce_tumor, features = gene03) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))



temp_cnv <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/temp_cnvtype.rds")
meta0 <- sce_tumor@meta.data
meta0$barcode <- sub(":", "-", rownames(meta0))
meta0 <- left_join(meta0, temp_cnv)
sce_tumor$ctype_cnv <- NULL
sce_tumor$ctype_cnv <- meta0$ctype_cnv
sce_tumor$ctype_cnv[is.na(sce_tumor$ctype_cnv)] <- sce_tumor$cell_types[is.na(sce_tumor$ctype_cnv)] 


Idents(sce_tumor) <- "ctype_cnv"
DotPlot(object = sce_tumor, assay = "SCT", features = unique(c(gene_dotplot, gene02))) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  scale_color_distiller(palette = "RdYlBu")
Idents(sce_tumor) <- "orig.ident"



Idents(sce_hd) <- "cell_types"
DotPlot(object = sce_hd, assay = "SCT", features = unique(c(gene_dotplot, gene02))) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  scale_color_distiller(palette = "RdYlBu")
Idents(sce_hd) <- "orig.ident"







  
  




sce_tumor <- NormalizeData(sce_tumor)

future::plan("multicore", workers = 20) ## do parallel
all_markers_tumor <- FindAllMarkers(sce_tumor)
future::plan("sequential")



cell_color <- c("HSC_MPP" = "indianred1", "CLP" = "indianred3", "pre_proB" = "darkseagreen1",
                "proB" = "darkolivegreen1", "preB" = "lightgreen",  "preB_I" = "yellowgreen",
                "preB_II" = "#3CB371", "immatureB" = "chartreuse4", "matureB" = "darkgreen",
                "erythroid_cell" = "lightpink2", "NK_T" = "#9575aa", "myeloid" = "deeppink4")

pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr", "igraph",
          "Seurat", "CytoTRACE", "reticulate", "monocle", "numDeriv"#, "monocle3"
)
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}

rna_assay <- GetAssayData(sce_tumor, assay = "RNA", slot = "data")
meta_mef2d <- sce_tumor@meta.data
# subset to remove cells
Idents(sce_tumor) <- "cell_types"
mef2d_sub <- subset(sce_tumor, `cell_types` %notin% c("NK_T", "erythroid_cell", "myeloid"))
cyt_res_mef2d_sub <- CytoTRACE(as.matrix(mef2d_sub[["RNA"]]@counts), ncores = 30)

meta_sub <- mef2d_sub@meta.data
mef2d_sub_data <- GetAssayData(mef2d_sub, assay = "RNA", slot = 'counts')
fd <- data.frame(gene_short_name = row.names(mef2d_sub_data), row.names = row.names(mef2d_sub_data))
pd <- new('AnnotatedDataFrame', data = meta_sub)
fd1 <- new('AnnotatedDataFrame', data = fd)

cds0 <- monocle::newCellDataSet(mef2d_sub_data,
                                phenoData = pd,
                                featureData = fd1,
                                lowerDetectionLimit = 0.5,
                                expressionFamily = VGAM::negbinomial.size())


cds0 <- cds0 %>% estimateSizeFactors() %>% estimateDispersions()
cds1 <- monocle::detectGenes(cds0, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(cds1), num_cells_expressed >= 0.1 * dim(mef2d_sub_data)[2])) 

# get genes from 

all_markers_fil <- all_markers_tumor[all_markers_tumor$p_val_adj < 0.05, ]
all_markers_fil$absFC <- abs(all_markers_fil$avg_log2FC)
all_markers_fil$absdiff <- abs(all_markers_fil$pct.1 - all_markers_fil$pct.2)
top_markers <- all_markers_fil %>% arrange(desc(absFC)) %>% group_by(cluster) %>% dplyr::slice(1:100)
top_markers2 <- all_markers_fil %>% arrange(desc(absdiff)) %>% group_by(cluster) %>% dplyr::slice(1:100)
top_markers3 <- all_markers_fil %>% arrange(desc(absFC)) %>% group_by(cluster) %>% dplyr::slice(1:150)
markers <- unique(c(top_markers3$gene, top_markers2$gene))
ordering_genes <- intersect(markers, expressed_genes)


cds2 <- monocle::setOrderingFilter(cds1[expressed_genes, ], ordering_genes)
cds_test2 <- monocle::reduceDimension(cds2, max_components = 2, method = 'DDRTree')
readr::write_rds(cds_test2, "/cluster/home/yjliu_jh/projects/mef2d/analysis/cds_ordered_m2y.rds")

#cds_test2 <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/analysis/cds_ordered_m2x.rds")
#source("/cluster/home/yjliu_jh/projects/mef2d/code/order_cells.R")

cds_test3 <- orderCells(cds_test2)
readr::write_rds(cds_test3, "/cluster/home/yjliu_jh/projects/mef2d/analysis/cds_ordered_m2yy.rds")

cds_test2 <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/analysis/cds_ordered_m20y.rds")

pm1 <- monocle::plot_cell_trajectory(cds_test2, color_by = "Pseudotime", size = 1, show_backbone = TRUE)
pm2 <- monocle::plot_cell_trajectory(cds_test2, color_by = "cell_types", size = 1, show_backbone = TRUE) +
  scale_color_discrete(cell_color)
pm_all <- pm1 + pm2 + plot_layout(nrow = 1)

pm2 + facet_wrap("cell_types")
pm2 + facet_wrap("orig.ident")

table(pData(cds_test3)[, c("cell_types", "State")])



ordering_genes <- intersect(top_markers3$gene, expressed_genes)
cds2 <- monocle::setOrderingFilter(cds1[expressed_genes, ], ordering_genes)
cds_test2 <- monocle::reduceDimension(cds2, max_components = 2, method = 'DDRTree')
cds_test2 <- orderCells(cds_test2)
readr::write_rds(cds_test2, "/cluster/home/yjliu_jh/projects/mef2d/analysis/cds_ordered_m20y.rds")







  
