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
















# yzy's code
# plot --------------------------------------------------------------------


work.dir <- "/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/tenx/pyscenic/"
mkdir <- function(dir.list){
  lapply(dir.list,function(x){
    if (!dir.exists(x)) {
      dir.create(x)
    }
  })
}
work.dir <- '/cluster/home/yzy_jh/projects/mef2d/analysis/zwn/human/tenx/pyscenic/'
setwd(work.dir)

od <- file.path(work.dir,"plot/")
checkdir(od)
projectname <- "mef2d"

# function

sce <- read_rds('~/projects/mef2d/analysis/zwn/human/tenx/annotation/sce.rds')
cellInfo <- sce@meta.data
rownames(cellInfo)<-gsub('-','.',rownames(cellInfo))#将sce和auc两个文件基因名格式统一

scenicLoomPath="~/projects/mef2d/analysis/zwn/human/tenx/pyscenic/sample_SCENIC.loom"
loom <- open_loom(scenicLoomPath)
# Read information from loom file: 
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name="RegulonsAUC") 


if(!is.null(sce@meta.data$seurat_clusters)){
  regulonActivity_bycluster <- sapply(split(rownames(cellInfo), cellInfo$seurat_clusters),
                                      function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
  colnames(regulonActivity_bycluster)<-factor(colnames(regulonActivity_bycluster))
  

  p.all <- pheatmap::pheatmap(regulonActivity_bycluster,
                              scale ='row',
                              cluster_cols = F,
                              show_rownames = T, fontsize_row=5, 
                              color=colorRampPalette(c("blue","white","red"))(100),
                              breaks=seq(-3, 3, length.out = 100),
                              treeheight_row=10, treeheight_col=10,
                              border_color=NA)
  p.50 <- pheatmap::pheatmap(regulonActivity_bycluster[sample(1:length(row.names(regulonActivity_bycluster)),50,replace = F),],
                             scale ='row',
                             cluster_cols = F,
                             show_rownames = T, #fontsize_row=3, 
                             color=colorRampPalette(c("blue","white","red"))(100),
                             breaks=seq(-3, 3, length.out = 100),
                             treeheight_row=10, treeheight_col=10,
                             border_color=NA)
  p.100 <- pheatmap::pheatmap(regulonActivity_bycluster[sample(1:length(row.names(regulonActivity_bycluster)),100,replace = F),],
                              scale ='row',
                              cluster_cols = F,
                              show_rownames = T, #fontsize_row=3, 
                              color=colorRampPalette(c("blue","white","red"))(100),
                              breaks=seq(-3, 3, length.out = 100),
                              treeheight_row=10, treeheight_col=10,
                              border_color=NA)
  # save image 
  ggsave(p.50,filename = glue("{od}/cluster_TF50.pdf"),width =5 ,height =10)
  ggsave(p.100,filename = glue("{od}/cluster_TF100.pdf"),width =5 ,height =20 )
  ggsave(p.all,filename = glue("{od}/cluster_TFall.pdf"),width =5 ,height =40 )
  
}

if(!is.null(sce@meta.data$cell_types)){
  regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$cell_types),
                                       function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
  colnames(regulonActivity_byCellType)<-factor(colnames(regulonActivity_byCellType))
  # plot heatmap 注意转录因子数量一般情况是200到300，如果少于100则修改参数
  p.all <- pheatmap::pheatmap(regulonActivity_byCellType,
                              scale ='row',
                              cluster_cols = F,
                              show_rownames = T, fontsize_row=5, 
                              color=colorRampPalette(c("blue","white","red"))(100),
                              breaks=seq(-3, 3, length.out = 100),
                              treeheight_row=10, treeheight_col=10,
                              border_color=NA)
  p.50 <- pheatmap::pheatmap(regulonActivity_byCellType[sample(1:length(row.names(regulonActivity_byCellType)),50,replace = F),],
                             scale ='row',
                             cluster_cols = F,
                             show_rownames = T, #fontsize_row=3, 
                             color=colorRampPalette(c("blue","white","red"))(100),
                             breaks=seq(-3, 3, length.out = 100),
                             treeheight_row=10, treeheight_col=10,
                             border_color=NA)
  p.100 <- pheatmap::pheatmap(regulonActivity_byCellType[sample(1:length(row.names(regulonActivity_byCellType)),100,replace = F),],
                              scale ='row',
                              cluster_cols = F,
                              show_rownames = T, #fontsize_row=3, 
                              color=colorRampPalette(c("blue","white","red"))(100),
                              breaks=seq(-3, 3, length.out = 100),
                              treeheight_row=10, treeheight_col=10,
                              border_color=NA)
  # save image
  ggsave(p.50,filename = glue("{od}/celltype_TF50.pdf"),width =5 ,height =10)
  ggsave(p.100,filename = glue("{od}/celltype_TF100.pdf"),width =5 ,height =20 )
  ggsave(p.all,filename = glue("{od}/celltype_TFall.pdf"),width =5 ,height =40 )
}























# ====================================================================
#mef2d_sub <- read_rds('/cluster/home/jhuang/projects/mef2d/analysis/zwn/human/tenx/seurat/merge/harmony/merged_anno_reg.rds')
merge <- read_rds("/cluster/home/jhuang/projects/mef2d/analysis/zwn/human/tenx/seurat/merge/merge.rds")
sce_tumor2 <- read_rds("/cluster/home/yjliu_jh/projects/mef2d/data/sce_tumor_final.rds")

DimPlot(mef2d_sub, reduction = "umap", group.by = "cell_types", split.by = "cell_types")

filtered_cells <- sce_tumor2$cell_id
new_cell_id <- str_sub(colnames(merge), 1 + str_locate(colnames(merge), ":")[, 1] / 2)
merge <- RenameCells(merge, new.names = new_cell_id)
merge_sub <- merge[, filtered_cells]
DimPlot(merge_sub, reduction = "umap", group.by = "cell_types_final", split.by = "cell_types_final")
DimPlot(merge_sub, reduction = "umap", group.by = "cell_types_final")
DimPlot(merge_sub, reduction = "umap", group.by = "cell_types")

# subset works, keep the umap structure

# 1 check current annotations

# see spatial clustering


facet(DimPlot(merge, reduction = "umap", group.by = "seurat_clusters", label = T), facet.by = "seurat_clusters")

Idents(merge_sub) <- "seurat_clusters"
DotPlot(object = merge_sub, features = c(gene03, "MME", "IGLL1")) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))

DotPlot(object = merge, features = c(gene03, "MME", "IGLL1")) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))



# 2 re-calculate the seurat clusters     resolution
# 3 

# get results from ddrtree and annotate cells with pesudotime

fbm_anno <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/data/singleR_anno_fbm.rds")
fbm_anno <- as.data.frame(fbm_anno[, 2:4])
fbm_anno$cell_id <- rownames(fbm_anno)
temp <- sce_tumor2@meta.data
temp2 <- data.frame(cell_id = rownames(merge_sub@meta.data), seu = merge_sub$seurat_clusters)
temp2 <- left_join(temp2, temp)
temp2$cell_types_old <- merge_sub@meta.data$cell_types
temp2 <- left_join(temp2, fbm_anno[, c("cell_id", "labels")])

table(temp2[, c("temp_anno", "labels")])
table(temp2[, c("temp_anno", "cell_types_final")])


temp2$temp_anno <- as.character(temp2$seu)
temp2$temp_anno[temp2$temp_anno %in% "9"] <- "T"
temp2$temp_anno[temp2$temp_anno %in% "16"] <- "erythroid_cell"
temp2$temp_anno[temp2$temp_anno %in% "22"] <- "myeloid"
temp2$temp_anno[temp2$temp_anno %in% "23"] <- "MLP"
temp2$temp_anno[temp2$temp_anno %in% "14"] <- "NK"
temp2$temp_anno[temp2$temp_anno %in% "NK" & temp2$cell_types_final %in% "T"] <- "NK_T"
temp2$temp_anno[temp2$temp_anno %in% c("4", "5", "25")] <- "preB_MKI67p" #
temp2$temp_anno[temp2$temp_anno %in% c("0", "6", "11")] <- "proB_IGHDp" #
temp2$temp_anno[temp2$temp_anno %in% c("2", "8", "24")] <- "pre_proB"
temp2$temp_anno[temp2$temp_anno %in% c("1", "17")] <- "preB_II"
temp2$temp_anno[temp2$temp_anno %in% c("12")] <- "contains_matureB"
temp2$temp_anno[temp2$temp_anno %in% c("13", "10", "15")] <- "preB_y"
temp2$temp_anno[temp2$temp_anno %in% c("19")] <- "preB_IGKCp"
temp2$temp_anno[temp2$temp_anno %in% c("7")] <- "preB_BLKp"
temp2$temp_anno[temp2$temp_anno %in% c("18", "20", "21")] <- "immatureB"
temp2$temp_anno[temp2$temp_anno %in% c("3")] <- "preB_x"



temp2$temp_anno <- as.character(temp2$seu)
temp2$temp_anno[temp2$temp_anno %in% "9"] <- "NK_T"
temp2$temp_anno[temp2$temp_anno %in% "16"] <- "erythroid_cell"
temp2$temp_anno[temp2$temp_anno %in% "22"] <- "myeloid"
temp2$temp_anno[temp2$temp_anno %in% "23"] <- "MLP"
temp2$temp_anno[temp2$temp_anno %in% "14"] <- "NK_T"
temp2$temp_anno[temp2$temp_anno %in% "NK" & temp2$cell_types_final %in% "T"] <- "NK_T"
temp2$temp_anno[temp2$temp_anno %in% c("4", "5", "25")] <- "preB_I" #
temp2$temp_anno[temp2$temp_anno %in% c("0", "6", "11")] <- "proB" #
temp2$temp_anno[temp2$temp_anno %in% c("2", "8", "24")] <- "pre_proB"
temp2$temp_anno[temp2$temp_anno %in% c("1", "17")] <- "preB_II"
temp2$temp_anno[temp2$temp_anno %in% c("12")] <- "matureB"
temp2$temp_anno[temp2$temp_anno %in% c("13", "10", "15")] <- "preB_I"
temp2$temp_anno[temp2$temp_anno %in% c("19")] <- "preB_I"
temp2$temp_anno[temp2$temp_anno %in% c("7")] <- "preB_I"
temp2$temp_anno[temp2$temp_anno %in% c("18", "20", "21")] <- "immatureB"
temp2$temp_anno[temp2$temp_anno %in% c("3")] <- "preB_I"






temp3 <- left_join(temp3, fbm_anno)


merge$temp_anno[temp3$labels %in% "ELP"] <- "ELP"

DotPlot(object = merge, features = c(gene03, "MME", "IGLL1")) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))





merge_sub$cn_cluster <- temp2$cn_cluster



DimPlot(merge_sub, reduction = "umap", group.by = "temp_anno", label = T)


gene04 <- setdiff(unique(c(gene03, gene02)), gene03[19:36])
gene05 <- c(gene04[-c(22:30)], "MME", "IGLL1", "KIT")

merge_sub$temp_anno <- temp2$temp_anno
Idents(merge_sub) <- "temp_anno"
DotPlot(object = merge_sub, features = c(gene03, "MME", "IGLL1", "ATXN1", "KIT", "CHCHD10")) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
Idents(merge_sub) <- "orig.ident"


merge_sub$temp_anno <- temp2$temp_anno
Idents(merge_sub) <- "temp_anno"
DotPlot(object = merge_sub, features = gene05) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
Idents(merge_sub) <- "orig.ident"






merge_sub$cell_types_final <- temp2$cell_types_final
merge_sub$cn_cluster <- temp2$cn_cluster
merge_sub$cell_types_broad <- temp2$cell_types_final_broad

meta_time <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/output/ddrtree_anno.rds")
  
DimPlot(merge_sub, reduction = "umap", group.by = "orig.ident", split.by = "orig.ident")



T_id <- temp2$cell_id[temp2$cell_types_final %in% "T"]



temp3 <- merge@meta.data
temp3$cell_id <- rownames(temp3)
temp3$temp_anno <- as.character(temp3$seurat_clusters)
temp3$temp_anno[temp3$temp_anno %in% "9"] <- "NK_T"
temp3$temp_anno[temp3$temp_anno %in% "16"] <- "erythroid_cell"
temp3$temp_anno[temp3$temp_anno %in% "22"] <- "myeloid"
temp3$temp_anno[temp3$temp_anno %in% "23"] <- "CLP"
temp3$temp_anno[temp3$temp_anno %in% "14"] <- "NK_T"
temp3$temp_anno[temp3$temp_anno %in% "NK" & temp3$cell_id %in% T_id] <- "NK_T"
temp3$temp_anno[temp3$temp_anno %in% c("4", "5", "25")] <- "preB_I" #
temp3$temp_anno[temp3$temp_anno %in% c("0", "6", "11")] <- "proB" #
temp3$temp_anno[temp3$temp_anno %in% c("2", "8", "24")] <- "proB"
temp3$temp_anno[temp3$temp_anno %in% c("1", "17")] <- "preB_II"
temp3$temp_anno[temp3$temp_anno %in% c("12")] <- "matureB"
temp3$temp_anno[temp3$temp_anno %in% c("13", "10", "15")] <- "preB_I"
temp3$temp_anno[temp3$temp_anno %in% c("19")] <- "preB_I"
temp3$temp_anno[temp3$temp_anno %in% c("7")] <- "preB_I"
temp3$temp_anno[temp3$temp_anno %in% c("18", "20", "21")] <- "immatureB"
temp3$temp_anno[temp3$temp_anno %in% c("3")] <- "preB_I"


merge$temp_anno <- temp3$temp_anno

temp3 <- left_join(temp3, fbm_anno)
merge$temp_anno[temp3$labels %in% "ELP"] <- "ELP-like"
merge$temp_anno[merge$temp_anno %in% "proB" & temp3$labels %in% "pre B progenitor"] <- "preB_I"
merge$temp_anno[merge$temp_anno %in% "NK_T" & temp3$labels %in% "naive B cell"] <- "immatureB"

merge$cell_types_broad <- merge$temp_anno



temp3 <- merge@meta.data
temp3$cell_id <- rownames(temp3)
temp3$temp_anno <- as.character(temp3$seurat_clusters)
temp3$temp_anno[temp3$temp_anno %in% "9"] <- "T"
temp3$temp_anno[temp3$temp_anno %in% "16"] <- "erythroid_cell"
temp3$temp_anno[temp3$temp_anno %in% "22"] <- "myeloid"
temp3$temp_anno[temp3$temp_anno %in% "23"] <- "MLP"
temp3$temp_anno[temp3$temp_anno %in% "14"] <- "NK"
temp3$temp_anno[temp3$temp_anno %in% "NK" & temp3$cell_id %in% T_id] <- "NK_T"
temp3$temp_anno[temp3$temp_anno %in% c("4", "5", "25")] <- "preB_MKI67p" #
temp3$temp_anno[temp3$temp_anno %in% c("0", "6", "11")] <- "proB_IGHDp" #
temp3$temp_anno[temp3$temp_anno %in% c("2", "8", "24")] <- "preproB-like"
temp3$temp_anno[temp3$temp_anno %in% c("1", "17")] <- "preB_II"
temp3$temp_anno[temp3$temp_anno %in% c("12")] <- "contains_matureB"
temp3$temp_anno[temp3$temp_anno %in% c("13", "10", "15")] <- "preB_y"
temp3$temp_anno[temp3$temp_anno %in% c("19")] <- "preB_IGKCp"
temp3$temp_anno[temp3$temp_anno %in% c("7")] <- "preB_BLKp"
temp3$temp_anno[temp3$temp_anno %in% c("18", "20", "21")] <- "immatureB"
temp3$temp_anno[temp3$temp_anno %in% c("3")] <- "preB_x"

merge$temp_anno <- temp3$temp_anno

temp3 <- left_join(temp3, fbm_anno)
merge$temp_anno[temp3$labels %in% "ELP"] <- "ELP-like"
merge$temp_anno[grepl("proB", merge$temp_anno) & temp3$labels %in% "pre B progenitor"] <- "preB_z"
merge$temp_anno[merge$temp_anno %in% "NK_T" & temp3$labels %in% "naive B cell"] <- "immatureB"

merge$cell_types_all <- merge$temp_anno


DimPlot(merge, reduction = "umap", group.by = "cell_types_all", split.by = "cell_types_all")


merge$temp_anno <- NULL

readr::write_rds(merge, "/cluster/home/yjliu_jh/projects/mef2d/data/merged_4_tumors.rds")
readr::write_rds(sce_hd, "/cluster/home/yjliu_jh/projects/mef2d/data/healthy_donors.rds")

sce_hd$cell_types[sce_hd$cell_types %in% "imatureB"] <- "immatureB"












Idents(merge) <- "temp_anno"
DotPlot(object = merge, features = c(gene05, "CD34", "CD27", "MPO", "CSF1R", "CD7", "CD244", "PROM1")) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
Idents(merge) <- "orig.ident"




Idents(sce_hd) <- "cell_types"
sce_hd <- subset(sce_hd, subset = cell_types %in% unique(names(table(sce_hd$cell_types))))


Idents(sce_hd) <- "cell_types"
DotPlot(object = sce_hd, features = unique(c(gene02, "CD34", "CD27", "MPO", "CSF1R", "CD7", "CD244", "PROM1"))) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
Idents(sce_hd) <- "orig.ident"






DimPlot(merge, reduction = "umap", group.by = "temp_anno", split.by = "temp_anno")


readr::write_rds(merge, "/cluster/home/yjliu_jh/projects/mef2d/data/sce0307.rds")
readr::write_rds(temp2, "/cluster/home/yjliu_jh/projects/mef2d/data/temp2.rds")
  

merged <- read_rds("/cluster/home/jhuang/projects/mef2d/analysis/zwn/human/tenx/seurat/merged/mef2d.rds")
Idents(merged) <- "orig.ident"
sce_hd <- subset(merged, subset = orig.ident %in% c("H1_CD19", "H1_CD34pCD19n", "H2_CD19", "H2_CD34pCD19n"))




t3_plot <- temp3 %>% group_by(orig.ident, temp_anno) %>% summarise(count = n())
ggplot(t3_plot, aes(fill = temp_anno, y = count, x = orig.ident)) + 
  geom_bar(position = "fill", stat = "identity")







sce_tumor <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/data/merged_4_tumors.rds")
merge2 <- merge %>% RunPCA() %>% FindNeighbors() %>%
  FindClusters(resolution = 2)
merge2 <- merge %>% FindClusters(resolution = 1.7)

meta2 <- merge2@meta.data
meta2$cell_id <- rownames(meta2)
meta3 <- left_join(meta2, temp2[, c("cell_id", "cell_types_final")])
meta3$temp <- temp3$temp_anno
meta3$temp[meta3$temp %in% "12" & meta3$SCT_snn_res.2 %in% 21] <- "matureB"
meta3$temp[meta3$temp %in% "12" & meta3$SCT_snn_res.2 %in% c(6, 10)] <- "preB_x"
meta3$temp[meta3$temp %in% "12" & meta3$SCT_snn_res.2 %in% c(0, 25, 27)] <- "proB_IGHDp"
meta3$temp[meta3$temp %in% "12" & meta3$SCT_snn_res.2 %in% c(1, 2, 5, 12, 26, 31)] <- "preB_MKI67p"
meta3$temp[meta3$temp %in% "12" & meta3$SCT_snn_res.2 %in% c(19)] <- "preB_II"
meta3$temp[meta3$temp %in% "12" & meta3$SCT_snn_res.2 %in% c(34)] <- "T"
meta3$temp[meta3$temp %in% "12" & meta3$SCT_snn_res.2 %in% c(11, 13, 16, 20)] <- "preB_y"
meta3$temp[meta3$temp %in% "12" & meta3$SCT_snn_res.2 %in% c(28)] <- "preB_BLKp"

meta3$cellxx <- sce_tumor$cell_types_all
meta3$temp[meta3$cellxx %in% "preB_z"] <- "preB_z"





merge@meta.data$cell_types_all <- meta3$temp
merge@meta.data$temp_anno <- NULL


tempu <- merge@reductions$umap@cell.embeddings
T_cells <- rownames(merge@meta.data)[merge$cell_types_all %in% "T"]
ex_id <- rownames(tempu[T_cells, ][tempu[T_cells, "UMAP_1"] < 9, ])
meta3$temp[colnames(merge) %in% ex_id] <- "preB_y"
merge@meta.data$cell_types_all <- meta3$temp


meta3$temp_broad <- meta3$temp
meta3$temp_broad[meta3$temp %in% c("NK", "NK_T", "T")] <- "NK_T"
meta3$temp_broad[meta3$temp %in% c("MLP")] <- "CLP"
meta3$temp_broad[meta3$temp %in% c("preB_BLKp", "preB_IGKCp", "preB_MKI67p", "preB_x", "preB_y", "preB_z")] <- "preB_I"
meta3$temp_broad[meta3$temp %in% c("preproB-like", "proB_IGHDp")] <- "proB"
merge@meta.data$cell_types_broad <- meta3$temp_broad


DimPlot(merge, reduction = "umap", group.by = "cell_types_broad")
DimPlot(merge, reduction = "umap", group.by = "cell_types")

readr::write_rds(merge, "/cluster/home/yjliu_jh/projects/mef2d/data/merged_4_tumors_new.rds")


mature_cells <- sce_tumor2$cell_id[sce_tumor2$cell_types_final %in% "matureB"]
sub_mature <- sce_tumor@meta.data[rownames(sce_tumor@meta.data) %in% mature_cells, ]
DimPlot(sce_tumor, reduction = "umap", group.by = "SCT_snn_res.1.6")




gene04 <- setdiff(unique(c(gene03, gene02)), gene03[19:36])
gene05 <- c(gene04[-c(22:30)], "MME", "IGLL1", "KIT")


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


