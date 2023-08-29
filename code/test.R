if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SingleCellExperiment")
a
# install devtools if necessary
install.packages('devtools')

# install the MuSiC package
devtools::install_github('xuranw/MuSiC')


BiocManager::install("Seurat")
n

xx <- read_h5ad("F:\\codeMEF2D\\mef2d\\data\\metacells.h5ad", backed = NULL)

mdata_mat <- as(as.matrix(xxxx), "dgCMatrix")

sce <- readr::read_rds("F:\\codeMEF2D\\sce_tumor_fil.rds")
meta <- sce@meta.data
meta0 <- meta %>% group_by(metacell) %>% mutate(perc = max(table(orig.ident)) / sum(table(orig.ident)))
meta00 <- meta0 %>% group_by(metacell) %>% summarise(id = names(table(orig.ident)[which.max(table(orig.ident))]))
mcxxx <- unique(meta0$metacell[meta0$perc < 0.5])
meta00$id[meta00$metacell %in% mcxxx] <- NA
meta00 <- meta00[order(meta00$metacell), ]

ctomc <- ctomc %>% group_by(metacell) 


rna_mat <- sce[["RNA"]]
sct_mat <- sce[["SCT"]]

rna <- as.matrix(rna_mat)

rnat<- t(as.matrix(rna_mat[]))
sctt<- t(sct_mat[])
identical(rownames(rnat), rownames(meta))

D <- data.frame(metacell = meta$metacell, sctt)
#small_rna <- aggregate(. ~ metacell, D, sum)
rna <- D %>% group_by(metacell) %>%
  summarise_at(vars(-group_cols()), sum)  ## 5 min in laptop
sct <- D %>% group_by(metacell) %>%
  summarise_at(vars(-group_cols()), sum)

rm(list = c("rna_mat", "sct_mat", "D"))
gc()

umis_rna <- rna %>% as.data.frame() %>% column_to_rownames(var = "metacell") %>% t() 
umis_sct <- sct %>% as.data.frame() %>% column_to_rownames(var = "metacell") %>% t() 

#fractions <- umis / rowSums(umis)

rownames(umis_rna) <- sub("\\.", "-", rownames(umis_rna))
rownames(umis_sct) <- sub("\\.", "-", rownames(umis_sct))

#log_fractions <- log2(1e-5 + fractions)
log_rna <- log2(1 + umis_rna)
log_rna <- t(log_rna)

log_sct <- log2(0.7 + umis_sct)
log_sct <- t(log_sct)




median_log_rna <- apply(log_rna, 2, median)
relative_log_rna <- sweep(log_rna, 2, median_log_rna)
median_log_sct <- apply(log_sct, 2, median)
relative_log_sct <- sweep(log_sct, 2, median_log_sct)




#pht6 <- readr::read_rds("F:\\codeMEF2D\\mef2d\\data\\plotht.rds")




mcn_order <- pht6x$tree_col$labels[pht6x$tree_col$order]
mcn_colclust <- cutree(pht6x$tree_col, k = 16)
mcn_cco <- mcn_colclust[mcn_order]

annotation_col3 <- unique(meta[, c("metacell", "cell_types_all", "pre_anno_pred",
                                   "labels_leu1_pred", "cell_types_pred")])
annotation_col3 <- annotation_col3[order(annotation_col3$metacell), ]
annotation_col3$cluster <- mcn_colclust
annotation_col3$id <- meta00$id

cluster_colors2 <- c("#ee5d5a", "#34956c", "#6b6498", "#71a3a2", "#c88978", "#cc9b32", "#53096a", "#03bb26", 
                              RColorBrewer::brewer.pal(8, "Set2"), RColorBrewer::brewer.pal(8, "Pastel1"), 
                              RColorBrewer::brewer.pal(8, "Set1"))[1:16]
names(cluster_colors2) <- unique(annotation_col3$cluster)
leu1_colors <- c("seagreen", "yellow3", "red", "brown3", "ivory", "lightgreen", "lightblue", "orange")
names(leu1_colors) <- sort(na.omit(unique(annotation_col3$labels_leu1_pred)))
pre_anno_colors <- c("magenta", "red", "green2", "green", "orange", "chartreuse", "limegreen", "lawngreen", "lightgreen")
names(pre_anno_colors) <- sort(na.omit(unique(annotation_col3$pre_anno_pred)))
cell_type_colors <- c("brown3", "red", "green", "purple", "orange", "limegreen", RColorBrewer::brewer.pal(4, "Paired"))
names(cell_type_colors) <- sort(na.omit(unique(annotation_col3$cell_types_pred)))
cluster_colors_all <- c("#ee5d5a", "#34956c", "#6b6498", "#71a3a2", "#c88978", "#cc9b32", "#53096a", "#03bb26", 
                              RColorBrewer::brewer.pal(8, "Set2"), RColorBrewer::brewer.pal(3, "Pastel1"))
names(cluster_colors_all) <- unique(annotation_col3$cell_types_all)
id_colors <- c("#c88978", "#cc9b32", "gray", RColorBrewer::brewer.pal(3, "Pastel1"))
names(id_colors) <- unique(annotation_col3$id)

annotation_colors4 <- list(cluster = cluster_colors2,
                           labels_leu1_pred = leu1_colors,
                           pre_anno_pred = pre_anno_colors,
                           cell_types_pred = cell_type_colors,
                           cell_types_all = cluster_colors_all,
                           id = id_colors)
annotation_col4 <- annotation_col3 %>% remove_rownames() %>% column_to_rownames(var = "metacell")


pht6 <- readr::read_rds("F:\\codeMEF2D\\mef2d\\data\\plotht.rds")
genex <- pht6x$tree_row$labels[pht6x$tree_row$order]

gene02 <- genex[-c(1:3, 5, 6 , 8, 11, 14, 26:28, 43:60, 77:82, 84:85)]
gene03 <- setdiff(gene02, c("MS4A4A", "IRF8", "BCL11B", "IL7R", "S100A9", "IGHG1", "MLLT3", "CD3G",
                            "SPI1", "HBD", "CD68"))
pht5data <-  t(relative_log_sct[, gene03])
pht5data <-  t(relative_log_rna[, gene01])



pht6data <- round(pht5data, 4)




breaks <- pracma::interp1(0:7, c(-3, -2, -1, 0, 1, 2, 3, 4), 0:140/20)
colors <- colorRampPalette(c('darkred', 'red', 'white', 'white', 'lightblue', 'blue', 'darkblue'))(141)
options(repr.plot.width = 25, repr.plot.height = 13)





pht6x <- pheatmap::pheatmap(
  pht6data,
  # cluster_rows = F,
  clustering_method = "mcquitty",
  treeheight_col= 80,
  treeheight_row=0,
  cellwidth=1,
  cellheight=10,
  show_rownames=TRUE,
  show_colnames=FALSE,  
  main='Metacell Clusters - BALL anno Genes',
  color=colors,
  breaks=breaks,
  #annotation_row=annotation_row3,
  annotation_col=annotation_col4,
  annotation_colors=annotation_colors4,
  annotation_legend=TRUE,
  legend=TRUE
)

pdf("D:\\test_heatmap666.pdf", width = 15, height = 20)
pht6x
dev.off()


pht6y <- pheatmap::pheatmap(
  pht6data,
  # cluster_rows = F,
  clustering_method = "ward.D2",
  treeheight_col=0,
  treeheight_row=0,
  cellwidth=1,
  cellheight=10,
  show_rownames=TRUE,
  show_colnames=FALSE,  
  main='Metacell Clusters - BALL anno Genes',
  color=colors,
  breaks=breaks,
  #annotation_row=annotation_row3,
  annotation_col=annotation_col4,
  annotation_colors=annotation_colors4,
  annotation_legend=TRUE,
  legend=TRUE
)

pdf("D:\\test_heatmap6666.pdf", width = 15, height = 20)
pht6y
dev.off()







readr::write_rds(annotation_col4, "D:\\ancol4.rds")
table(annotation_col4[, c("cluster", "cell_types_all")])


quick_anno <- function(sce, anno){
  temp <- sce@meta.data
  temp <- left_join(temp, data.frame(metacell = 0:823, anno = as.character(anno)))
  sce@meta.data$anno <- temp$anno
  sce
}
sce$anno <- NULL
sce <- quick_anno(sce, annotation_col4$cluster)
dimcolors <- cluster_colors2
names(dimcolors) <- as.character(names(dimcolors))



ggpubr::facet((DimPlot(sce, group.by = "anno", cols = dimcolors)), facet.by = "anno")
ggpubr::facet(DimPlot(sce, group.by = "cell_types_all"), facet.by = "cell_types_all")
DimPlot(sce, group.by = "orig.ident")
FeaturePlot(sce, features = c("cnv1", "DNTT", "LEF1", "EBF1", "STMN1", "CD19"))



sce <- sce %>% RunUMAP(min.dist = 0.3, dims = 1:10#, umap.method = umap-learn
                       )

sce$



