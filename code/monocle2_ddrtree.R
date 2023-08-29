
pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr", "igraph",
          "Seurat", "CytoTRACE", "reticulate", "monocle", "numDeriv"#, "monocle3"
)
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}


mef2d <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/data/sce_tumor_final.rds")
meta_mef2d <- mef2d@meta.data

rna_assay <- GetAssayData(mef2d, assay = "RNA", slot = "data")
cell_color <- c("CLP" = "indianred3", "pre_proB" = "darkseagreen1",
                "proB" = "lightgreen", "preB_I" = "yellowgreen",
                "preB_II" = "chartreuse4", "immatureB" = "#3CB371", "matureB" = "darkgreen",
                "erythroid_cell" = "lightpink2", "NK_T" = "#9575aa", "myeloid" = "deeppink4",
                "tumor_cell3" = "lightseagreen")
samples <- unique(meta_mef2d$orig.ident)
sample_color <- c("#ee5d5a", "#9f2f6a", "#34956c", "#6b6498","#b29bc9", "#71a3a2", 
                           "#c88978",  "#a15217", "#ce662f", "#cc9b32", "#53096a", "#6569c2", "#66676c")[1:length(samples)]
names(sample_color) <- samples

mef2d_sub <- subset(mef2d, `cell_types_final_broad` %notin% c("NK_T", "erythroid_cell", "myeloid", "CLP"))




cyt_res_mef2d_sub <- CytoTRACE(as.matrix(mef2d_sub[["RNA"]]@counts), ncores = 12)

cell_type <- as.character(mef2d_sub@meta.data$cell_types_final_broad)
orig_ident <- as.character(mef2d_sub@meta.data$orig.ident)
names(cell_type) <-rownames(mef2d_sub@meta.data)
names(orig_ident) <- rownames(mef2d_sub@meta.data)

meta_sub <- mef2d_sub@meta.data
emb <- as.matrix(mef2d_sub@reductions$umap@cell.embeddings)
mat <- t(cyt_res_mef2d_sub$exprMatrix)[rownames(emb), ]
cyto <- cyt_res_mef2d_sub$CytoTRACE[rownames(emb)]

cyto_plot_datfm <- data.frame(emb, CytoTRACE_score = cyto, HDAC9 = mat[, "HDAC9"], 
                              MEF2D = mat[, "MEF2D"], meta_sub)

# CytoTRACE plot by ggplot2
temp_color <- RColorBrewer::brewer.pal(11, "Spectral")
temp_color[6] <- "gold"
  rbPal <- colorRampPalette(temp_color)
  
  p_cyt <- ggplot(cyto_plot_datfm, aes(x = UMAP_1, y = UMAP_2, color = CytoTRACE_score)) + geom_point(size = .5) + 
    ggpubr::theme_pubr() + scale_colour_gradientn(name = "CytoTRACE score", colours = rev(rbPal(50)),
                                                  guide = ggplot2::guide_colourbar(ticks.colour = "black", ticks.linewidth = 1, frame.colour = "black"),
                                                  breaks = seq(0, 1, 0.2), labels = c("0.0 (More diff.)", 0.2, 0.4, 0.6, 0.8, "1.0 (Less diff.)")) +
    theme(legend.text = ggplot2::element_text(size = 6), legend.title = ggplot2::element_text(size = 7), 
          plot.title = ggplot2::element_text(size = 10, hjust = 0.5), axis.title.x = ggplot2::element_text(size = 8), 
          axis.title.y = ggplot2::element_text(size = 7), axis.text = ggplot2::element_text(size = 6), 
          legend.position = "right", plot.margin = ggplot2::unit(c(0.5, 1, 0.5, 1), "cm"))
  
  

  
  
  "reducedDimW<-" <- function (cds, value)
  {
    stopifnot(is(cds, "CellDataSet"))
    cds@reducedDimW <- value
    validObject(cds)
    cds
  }
  
  
  "reducedDimK<-" <- function (cds, value)
  {
    stopifnot(is(cds, "CellDataSet"))
    cds@reducedDimK <- value
    validObject(cds)
    cds
  }
  
# --- monocle2 
  
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
diff_test_res <- monocle::differentialGeneTest(cds1[expressed_genes, ],
                                               fullModelFormulaStr = "~cell_types_final_broad", cores = 1)
deg <- subset(diff_test_res, qval < 0.001)
deg <- deg[order(deg$qval, decreasing = F), ]
ordering_genes <- row.names(deg[1:1000,])
cds2 <- monocle::setOrderingFilter(cds1[expressed_genes, ], ordering_genes)
cds_test2 <- monocle::reduceDimension(cds2, max_components = 2, method = 'DDRTree')
cds_test2 <- orderCells(cds_test2)
readr::write_rds(cds_test2, "/cluster/home/yjliu_jh/projects/mef2d/analysis/cds_mef_final.rds")

pdf("/cluster/home/yjliu_jh/projects/mef2d/trajectory_ddrtree_monocle2.pdf", width = 5, height = 5)
monocle::plot_cell_trajectory(cds_test2, color_by = "Pseudotime", size = 1, show_backbone = TRUE)
dev.off()


pm1 <- monocle::plot_cell_trajectory(cds_test2, color_by = "Pseudotime", size = 1, show_backbone = TRUE)
pm2 <- monocle::plot_cell_trajectory(cds_test2, color_by = "cell_types_final_broad", size = 1, show_backbone = TRUE) +
  scale_color_discrete(cell_color)
pm_all <- pm1 + pm2 + plot_layout(nrow = 1)

pdf("/cluster/home/yjliu_jh/projects/mef2d/trajectory_ddrtree_xxx.pdf", width = 10, height = 5)
pm_all
dev.off()



cds_test2 <- read_rds("/cluster/home/yjliu_jh/projects/mef2d/analysis/cds_mef_final.rds")

pm3 <- monocle::plot_cell_trajectory(cds_test2, color_by = "cell_types_final_broad", size = 1, show_backbone = TRUE) +
  scale_color_discrete(cell_color) + facet_wrap("cell_types_final_broad")


pm3





