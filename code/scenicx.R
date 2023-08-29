pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse",
          "SCopeLoomR", "AUCell", "SCENIC", "KernSmooth", "BiocParallel", "Seurat",
          "plotly", "ComplexHeatmap", "data.table", "grid")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}


scenicLoomPath <- "/cluster/home/yjliu_jh/projects/mef2d/analysis/zwn/home/SCENIC/sample_SCENIC_all.loom"

loom <- open_loom(scenicLoomPath)
# Read information from loom file:
exprMat <- get_dgem(loom)
exprMat_log <- log2(exprMat+1) # Better if it is logged/normalized
regulons_incidMat <- get_regulons(loom, column.attr.name = "Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name = "RegulonsAUC")
regulonAucThresholds <- get_regulon_thresholds(loom)
embeddings <- get_embeddings(loom)
cellClusters <- get_clusterings(loom)
close_loom(loom)

# length(regulons);  head(names(regulons))
# regulonAUC
# length(regulonsAucThresholds)
plot(embeddings$`SCENIC AUC UMAP`)




cells <- colnames(exprMat)
cellClusters <- data.frame(cell_type = sce_tumor$cell_types_all, 
                           cell_type_broad = sce_tumor$cell_types_broad)

rss <- calcRSS(AUC = getAUC(regulonAUC), cellAnnotation = cellClusters[colnames(regulonAUC), "cell_type"])
# cellClusters: one-column data.frame(rownames = cell_ids, colname = 0?)
# we can alter and modify this step as we need
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)
plotRSS_oneSet(rss, setName = "myeloid")

# sub_regulonAUC

table(cellClusters$cell_type)


# regulon separated
cellsPerCluster <- split(rownames(cellClusters), cellClusters[, "cell_type"]) 
regulonAUCx <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)), ]
# Calculate average expression:
regulonActivity_byCellType <- sapply(cellsPerCluster,
                                     function(cells) rowMeans(getAUC(regulonAUCx)[, cells]))
# Scale expression:
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale = T))
# plot:
hm <- draw(ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled[70:120, ], name = "Regulon activity",
                                   row_names_gp = grid::gpar(fontsize=6))) # row font size
regulonOrder <- rownames(regulonActivity_byCellType_Scaled)[row_order(hm)] # to save the clustered regulons for later


# check top regulators
topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators$CellType <- factor(as.character(topRegulators$CellType))
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0), ]
dim(topRegulators)




# ====== separate mef2d and tcf3 samples ======

mef_id <- colnames(sce_tumor)[sce_tumor$orig.ident %in% c("M2", "M3", "M4")]
regulonAUC_mef <- regulonAUC[, match(mef_id, colnames(regulonAUC))]
expr_mef <- exprMat[, mef_id]

cells_mef <- colnames(expr_mef)
cellClusters_mef <- data.frame(cell_type = sce_tumor$cell_types_all[match(mef_id, colnames(sce_tumor))], 
                               cell_type_broad = sce_tumor$cell_types_broad[match(mef_id, colnames(sce_tumor))])
rss <- calcRSS(AUC = getAUC(regulonAUC_mef), cellAnnotation = cellClusters_mef[colnames(regulonAUC_mef), "cell_type_broad"])
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)

# just like this...


tcf_id <- colnames(sce_tumor)[sce_tumor$orig.ident %in% c("E1")]
regulonAUC_tcf <- regulonAUC[, match(tcf_id, colnames(regulonAUC))]
expr_tcf <- exprMat[, tcf_id]

cells_tcf <- colnames(expr_tcf)
cellClusters_tcf <- data.frame(cell_type = sce_tumor$cell_types_all[match(tcf_id, colnames(sce_tumor))], 
                               cell_type_broad = sce_tumor$cell_types_broad[match(tcf_id, colnames(sce_tumor))])



# regulon separated
cellsPerCluster <- split(rownames(cellClusters), cellClusters[, "cell_type"]) 
regulonAUCx <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)), ]
# Calculate average expression:
regulonActivity_byCellType <- sapply(cellsPerCluster,
                                     function(cells) rowMeans(getAUC(regulonAUCx)[, cells]))
# Scale expression:
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale = T))
# plot:
hm <- draw(ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled[1:60, ], name = "Regulon activity",
                                   row_names_gp = grid::gpar(fontsize=6))) # row font size
regulonOrder <- rownames(regulonActivity_byCellType_Scaled)[row_order(hm)] # to save the clustered regulons for later


# check top regulators
topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators$CellType <- factor(as.character(topRegulators$CellType))
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0), ]
dim(topRegulators)


