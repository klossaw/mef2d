pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr", "igraph",
          "Seurat", "anndata", "stats", "pheatmap", "pracma", "SeuratDisk"
)
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}

mdata <- anndata::read_h5ad('/cluster/home/yjliu_jh/projects/mef2d/data/metacell/mef2d_tumor_metacells.h5ad')

mdata_mat <- as(as.matrix(mdata$X), "dgCMatrix")
umis <- as.matrix(mdata$X)
fractions <- umis / rowSums(umis)     ## for each metacell, sum = 1
log_fractions <- log2(1e-5 + fractions)

genes_per_metacell <- 2
minimal_max_log_fraction_of_interesting_genes <- -10
minimal_relative_log_fraction_of_candidate_genes <- 2

median_log_fractions_of_genes <- apply(log_fractions, 2, median)
relative_log_fractions <- sweep(log_fractions, 2, median_log_fractions_of_genes) ## indicates fraction order among metacells
max_log_fractions_of_genes <- apply(log_fractions, 2, max)

interesting_genes_mask <- (max_log_fractions_of_genes
                           >= minimal_max_log_fraction_of_interesting_genes)
marker_genes <- unique(
  as.character(
    unlist(
      apply(
        relative_log_fractions[,interesting_genes_mask],
        1,
        function(relative_log_fraction_of_metacell) {
          candidate_genes_mask <- (relative_log_fraction_of_metacell
                                   >= minimal_relative_log_fraction_of_candidate_genes)
          names(
            head(
              sort(-relative_log_fraction_of_metacell[candidate_genes_mask]),
              n=genes_per_metacell
            )
          )
        }
      )
    )
  )
)



# some of these (e.g. JUN) are suspicious genes! 
forbidden_gene <- mdata$var$forbidden_gene
names(forbidden_gene) <- as.character(mdata$var_names)
forbidden_marker_gene <- unlist(
  lapply(
    forbidden_gene[marker_genes],
    function(value) { if (value) { 'Forbidden' } else { 'Allowed' } }
  )
)
annotation_row <- data.frame(forbidden=forbidden_marker_gene)
rownames(annotation_row) <- marker_genes
forbidden_colors <- c('gray', 'black')
names(forbidden_colors) <- c('Allowed', 'Forbidden')
annotation_colors <- list(forbidden=forbidden_colors)


breaks <- pracma::interp1(0:7, c(-3, -2, -1, 0, 1, 2, 3, 4), 0:140/20)
colors <- colorRampPalette(c('darkblue', 'blue', 'lightblue', 'white', '#ffcccb', 'red', 'darkred', 'darkred'))(141)
options(repr.plot.width = 25, repr.plot.height = 13)


pht <- pheatmap::pheatmap(
  t(relative_log_fractions[,marker_genes]),
  treeheight_col=0,
  treeheight_row=0,
  cellwidth=1,
  cellheight=10,
  show_rownames=TRUE,
  show_colnames=FALSE,
  main='Metacell Marker Genes',
  color=colors,
  breaks=breaks,
  annotation_row=annotation_row,
  
  annotation_colors=annotation_colors,
  annotation_legend=TRUE,
  legend=TRUE
)

pdf("/cluster/home/yjliu_jh/projects/mef2d/output/test_heatmap_tumor.pdf", width = 15, height = 12)
pht
dev.off()




# ===== use BALL genes to draw heatmap ======

mdata$obs$cluster <- as.character(cutree(pht$tree_col, k = 7))

gg <- ggpubr::ggscatter(mdata$obs, "umap_x", "umap_y", color = "cluster", size = 0.4)
pdf("/cluster/home/yjliu_jh/projects/mef2d/output/test_umap_tumor.pdf", width = 8, height = 8)
gg
dev.off()



ball_genes <- readxl::read_excel("/cluster/home/yjliu_jh/projects/mef2d/data/public/ball_clu_genes_32470390.xlsx",
                                 sheet = 2)
ball_genes_fil <- ball_genes[ball_genes$gene %in% colnames(relative_log_fractions), ]
ball_genes_fil <- ball_genes_fil[!(ball_genes_fil$gene %in% "FLT3" & ball_genes_fil$celltype %in% "HSPC"), ]



annotation_row2 <- data.frame(ball=ball_genes_fil$celltype)
rownames(annotation_row2) <- ball_genes_fil$gene
ball_colors <- RColorBrewer::brewer.pal(6, "Set2")
names(ball_colors) <- unique(ball_genes_fil$celltype)


annotation_col <- data.frame(cluster = mdata$obs$cluster)
rownames(annotation_col) <- rownames(mdata$obs)
cluster_colors <- c("#ee5d5a", "#34956c", "#6b6498", "#71a3a2", "#c88978", "#cc9b32", "#53096a")
names(cluster_colors) <- unique(mdata$obs$cluster)

annotation_colors2 <- list(cluster = cluster_colors, ball = ball_colors)

pht2data <- t(relative_log_fractions[, ball_genes_fil$gene])

pht2 <- pheatmap::pheatmap(
  pht2data,
  treeheight_col=0,
  treeheight_row=0,
  cellwidth=1,
  cellheight=10,
  show_rownames=TRUE,
  show_colnames=FALSE,
  main='Metacell Marker Genes',
  color=colors,
  breaks=breaks,
  annotation_row=annotation_row2,
  annotation_col=annotation_col,
  annotation_colors=annotation_colors2,
  annotation_legend=TRUE,
  legend=TRUE
)

pdf("/cluster/home/yjliu_jh/projects/mef2d/output/test_heatmap2_tumor.pdf", width = 15, height = 18)
pht2
dev.off()


# remove genes with no variation
rowclust = hclust(dist(pht2data))
reordered = pht2data[rowclust$order, ]

tempa <- apply(pht2data, 1, \(x) max(x) - min(x))
fil_genes <- names(tempa[tempa > 4])

# ======= 3rd round heatmap ======

ball_genes_fil2 <- ball_genes_fil[ball_genes_fil$gene %in% fil_genes, ]

annotation_row3 <- data.frame(ball=ball_genes_fil2$celltype)
rownames(annotation_row3) <- ball_genes_fil2$gene

pht3data <- t(relative_log_fractions[, ball_genes_fil2$gene])

pht3 <- pheatmap::pheatmap(
  pht3data,
  treeheight_col = 2,
  treeheight_row = 0,
  cellwidth = 1,
  cellheight = 10,
  show_rownames=TRUE,
  show_colnames=FALSE,
  main='Metacell Clusters - BALL anno Genes',
  color=colors,
  breaks=breaks,
  annotation_row=annotation_row3,
  annotation_col=annotation_col,
  annotation_colors=annotation_colors2,
  annotation_legend=TRUE,
  legend=TRUE
)

pdf("/cluster/home/yjliu_jh/projects/mef2d/output/test_heatmap3_tumor.pdf", width = 15, height = 8)
pht3
dev.off()


mc_order <- pht3$tree_col$labels[pht3$tree_col$order]
mc_colclust <- cutree(pht3$tree_col, k = 18)
mc_cco <- mc_colclust[mc_order]
# check subgroup for T markers




# === some other tests before the final heatmap


cells_temp <- anndata::read_h5ad('/cluster/home/yjliu_jh/projects/mef2d/data/metacell/mef2d_tumor_cells.h5ad')
ctomc2 <- cells_temp$obs
ctomc2$cell_id <- rownames(ctomc2)
mta <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/output/mta.rds")
ctomc2 <- left_join(ctomc2, mta[, -c(2, 6)])

calc_lhood <- function(df, a = "metacell", b){
  ## note that .data is not equal to .
  df <- df %>% group_by(.data[[a]]) %>% 
    mutate(pred_b = ifelse(is.null(names(which.max(table(.data[[b]])))), NA, names(which.max(table(.data[[b]])))),
           perc_b = max(table(.data[[b]])) / sum(table(.data[[b]])))
  df_sub <- df[df$perc_b > 0.5, ]
  df_s <- df_sub %>% group_by(.data[[a]]) %>% summarise(pred = unique(pred_b))
  colnames(df_s)[2] <- paste0(b, "_pred")
  df_s
}

annocol <- data.frame(metacell = rownames(annotation_col), cluster = annotation_col$cluster)
annocol$metacell <- as.numeric(annocol$metacell)

construct_anno <- function(df, columns){
  for (i in 1:length(columns)){
    temp_join <- calc_lhood(ctomc2, b = columns[i])
    df <- left_join(df, temp_join)
  }
  df <- df %>% column_to_rownames(var = "metacell")
  df
}

temp <- readr::read_rds("/cluster/home/ylxie_jh/projects/leukemia/analysis/weinazhang/human/public/final/mef2d_tumor_preanno.rds")
ctomc2 <- left_join(ctomc2, data.frame(cell_id = rownames(temp@meta.data), pre_anno = temp@meta.data$pre_anno))
readr::write_rds(ctomc2, "/cluster/home/yjliu_jh/projects/mef2d/output/cells_to_mc_tumor.rds")

annotation_col3 <- construct_anno(annocol, c("cell_types", "labels_fbm", ## warnings if NA
                                  "labels_leu1", "labels_leu1s", "pre_anno"))




tempb <- names(mc_colclust[mc_colclust %in% c(10, 18, 14, 9)])
table(ctomc2$cell_types[ctomc2$metacell %in% tempb])  ## almost all nk_t

to_join2 <- data.frame(metacell = as.numeric(names(mc_colclust)), clu18 = mc_colclust)
ctomc2 <- as.data.frame(left_join(ctomc2, to_join2))

#table(ctomc2$clu18[ctomc2$cn %in% "altered"])
# 6 and 11, over half
# pax5? ATG7?



sce <- readr::read_rds("~/projects/mef2d/analysis/zwn/human/tenx/seurat/merged/mef2d.rds")
sce_tumor <- subset(sce, `orig.ident` %in% c("B069_Dx_CD19", "M2", "M3", "M4", "E1"))
ocells <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/output/cells_deletion_pattern_infercnv.rds")
summary(tumor_exp["PAX5", colnames(tumor_exp) %in% ocells])


# ========

annotation_col3$cluster <- mc_colclust
cluster_colors2 <- c("#ee5d5a", "#34956c", "#6b6498", "#71a3a2", "#c88978", "#cc9b32", "#53096a", 
                              RColorBrewer::brewer.pal(7, "Set2"), RColorBrewer::brewer.pal(4, "Pastel1"))
names(cluster_colors2) <- unique(annotation_col3$cluster)
leu1_colors <- c("seagreen", "yellow3", "red", "brown3", "ivory", "lightgreen", "lightblue", "orange")
names(leu1_colors) <- sort(na.omit(unique(annotation_col3$labels_leu1_pred)))
pre_anno_colors <- c("magenta", "red", "green2", "green", "orange", "chartreuse", "limegreen", "lawngreen", "lightgreen")
names(pre_anno_colors) <- sort(na.omit(unique(annotation_col3$pre_anno_pred)))
cell_type_colors <- c("brown3", "red", "green", "purple", "orange", "limegreen", RColorBrewer::brewer.pal(4, "Paired"))
names(cell_type_colors) <- sort(na.omit(unique(annotation_col3$cell_types_pred)))
fbm_colors <- c("orange", "orangered", "gold", "red2", "brown3", "blue", "green2", "limegreen", "chartreuse", "lightgreen")
names(fbm_colors) <- sort(na.omit(unique(annotation_col3$labels_fbm_pred)))

annotation_colors3 <- list(cluster = cluster_colors2,
                           ball = ball_colors,
                           labels_leu1_pred = leu1_colors,
                           labels_leu1s_pred = leu1_colors,
                           pre_anno_pred = pre_anno_colors,
                           cell_types_pred = cell_type_colors,
                           labels_fbm_pred = fbm_colors)




pht5 <- pheatmap::pheatmap(
  pht3data,
  # cluster_rows = F,
  #clustering_method = "ward.D2",
  treeheight_col=0,
  treeheight_row=0,
  cellwidth=1,
  cellheight=10,
  show_rownames=TRUE,
  show_colnames=FALSE,  
  main='Metacell Clusters - BALL anno Genes',
  color=colors,
  breaks=breaks,
  annotation_row=annotation_row3,
  annotation_col=annotation_col3,
  annotation_colors=annotation_colors3,
  annotation_legend=TRUE,
  legend=TRUE
)

pdf("/cluster/home/yjliu_jh/projects/mef2d/output/test_heatmap5y_tumor.pdf", width = 18, height = 22)
pht5
dev.off()


mcn_order <- pht5$tree_col$labels[pht5$tree_col$order]
mcn_colclust <- cutree(pht5$tree_col, k = 18)
mcn_cco <- mcn_colclust[mcn_order]

to_join3 <- data.frame(metacell = as.numeric(names(mcn_colclust)), clust = mcn_colclust)
a3 <- cbind(annotation_col3, to_join3)
a3 <- left_join(ctomc2[, c("cell_id", "metacell")], a3)
colnames(a3)[3] <- "cluster_new"
readr::write_rds(a3, "/cluster/home/yjliu_jh/projects/mef2d/output/mc_tumor_clust.rds")



# ======


annotation_col3$cluster <- mcn_colclust
cluster_colors2 <- c("#ee5d5a", "#34956c", "#6b6498", "#71a3a2", "#c88978", "#cc9b32", "#53096a", 
                              RColorBrewer::brewer.pal(7, "Set2"), RColorBrewer::brewer.pal(4, "Pastel1"))
names(cluster_colors2) <- unique(annotation_col3$cluster)
leu1_colors <- c("seagreen", "yellow3", "red", "brown3", "ivory", "lightgreen", "lightblue", "orange")
names(leu1_colors) <- sort(na.omit(unique(annotation_col3$labels_leu1_pred)))
pre_anno_colors <- c("magenta", "red", "green2", "green", "orange", "chartreuse", "limegreen", "lawngreen", "lightgreen")
names(pre_anno_colors) <- sort(na.omit(unique(annotation_col3$pre_anno_pred)))
cell_type_colors <- c("brown3", "red", "green", "purple", "orange", "limegreen", RColorBrewer::brewer.pal(4, "Paired"))
names(cell_type_colors) <- sort(na.omit(unique(annotation_col3$cell_types_pred)))
fbm_colors <- c("orange", "orangered", "gold", "red2", "brown3", "blue", "green2", "limegreen", "chartreuse", "lightgreen")
names(fbm_colors) <- sort(na.omit(unique(annotation_col3$labels_fbm_pred)))

annotation_colors3 <- list(cluster = cluster_colors2,
                           ball = ball_colors,
                           labels_leu1_pred = leu1_colors,
                           labels_leu1s_pred = leu1_colors,
                           pre_anno_pred = pre_anno_colors,
                           cell_types_pred = cell_type_colors,
                           labels_fbm_pred = fbm_colors)



pht4data <- t(relative_log_fractions[, colnames(relative_log_fractions) %in% fil_genes2])

pht5 <- pheatmap::pheatmap(
  pht4data,
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
  annotation_row=annotation_row3,
  annotation_col=annotation_col3,
  annotation_colors=annotation_colors3,
  annotation_legend=TRUE,
  legend=TRUE
)

pdf("/cluster/home/yjliu_jh/projects/mef2d/output/test_heatmap5z_tumor.pdf", width = 18, height = 22)
pht5
dev.off()




# unique(mcn_cco)
# [1] 13 16 17 15 10  9  6 18  5  4  2  1 14 11  8  3  7 12

co1 <- c(13, 16, 17, 15, 10,  9,  6, 18,  5,  4,  2, 1, 14, 11,  8,  3,  7, 12)
co2 <- c("ATXN1_high_BCL11Bn_preB", "preB_IGKCp", "IGKCn_SPINK2p_proB", 
         "preB_x", "T", "NK", "ery", "ery-like", "naiveB1", "naiveB2", "myeloid", "matureB1",
         "matureB", "proB-like_HLADRn", "proB_x", "CLP", "prob_y", "prob_z")

# b1 pre b  b2 pro b



to_join_cl <- data.frame(cluster = co1, cluster_anno = co2)
temp <- annotation_col3
temp <- left_join(temp, to_join_cl)
annotation_col3$cluster_anno = temp$cluster_anno



# ====== BIG plot ======

temp_gene <- unique(c(add_genes, names(tempa)))
pht5data <- t(relative_log_fractions[, colnames(relative_log_fractions) %in% temp_gene])
tempc <- apply(pht5data, 1, \(x) max(x) - min(x))
fil_genes3 <- names(tempc[tempc > 4])

pht5data <- t(relative_log_fractions[, colnames(relative_log_fractions) %in% fil_genes3])

cluster_colors3 <- cluster_colors2
cluster_colors4 <- cluster_colors2
cluster_colors5 <- cluster_colors2
names(cluster_colors3) <- unique(annotation_col3$cluster_anno)
names(cluster_colors4) <- unique(annotation_col3$clustx)
names(cluster_colors5) <- unique(annotation_col3$clustx1)

annotation_colors4 <- list(cluster_anno = cluster_colors3,
                           ball = ball_colors,
                           labels_leu1_pred = leu1_colors,
                           labels_leu1s_pred = leu1_colors,
                           pre_anno_pred = pre_anno_colors,
                           cell_types_pred = cell_type_colors,
                           labels_fbm_pred = fbm_colors#,
                           #clustx = cluster_colors4,
                           #clustx1 = cluster_colors5
                           )

pht6 <- pheatmap::pheatmap(
  pht5data,
  #cluster_rows = F,
  #clustering_method = "ward.D2",
  treeheight_col=3,
  treeheight_row=0,
  cellwidth=1,
  cellheight=10,
  show_rownames=TRUE,
  show_colnames=FALSE,  
  main='Metacell Clusters - BALL anno Genes',
  color=colors,
  breaks=breaks,
  annotation_row=annotation_row3,
  annotation_col=annotation_col3,
  annotation_colors=annotation_colors4,
  annotation_legend=TRUE,
  legend=TRUE
)

pdf("/cluster/home/yjliu_jh/projects/mef2d/output/test_heatmap6ae_tumor.pdf", width = 18, height = 22)
pht6
dev.off()

readr::write_rds(pht6, "/cluster/home/yjliu_jh/projects/mef2d/output/plotht.rds")





nc_order <- pht6$tree_col$labels[pht6$tree_col$order]
nc_colclust <- cutree(pht6$tree_col, k = 18)
nc_cco <- nc_colclust[nc_order]

to_join4 <- data.frame(metacell = as.numeric(names(nc_colclust)), clustx1 = nc_colclust)
a4 <- cbind(annotation_col3, to_join4)
a4 <- left_join(ctomc2[, c("cell_id", "metacell")], a4)
colnames(a4)[3] <- "cluster_final"



a4 <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/output/anno_metacell.rds")
temp <- sce_tumor@meta.data
temp <- left_join(temp, a4[, c("cell_id", "clustx1")])
sce_tumor@meta.data$cluster_final <- temp$clustx1
readr::write_rds(sce_tumor, "/cluster/home/yjliu_jh/projects/mef2d/output/sce_tumor.rds")



mdata$obs$clust <- as.character(mcn_colclust)
gg <- ggpubr::ggscatter(mdata$obs, "umap_x", "umap_y", color = "clust", size = 0.4)

pdf("/cluster/home/yjliu_jh/projects/mef2d/output/test_umap_1.pdf", width = 8, height = 8)
gg
dev.off()

#ggsave("/cluster/home/yjliu_jh/projects/mef2d/output/test_umap.png", gg,
#       width = 800, height = 800, units = "px")




cells_temp <- anndata::read_h5ad("/cluster/home/yjliu_jh/projects/mef2d/data/metacell/mef2d_tumor_cells.h5ad")
ctomc2 <- cells_temp$obs
ctomc2$cell_id <- rownames(ctomc2)
ctomc2 <- ctomc2 %>% group_by(metacell) %>% mutate(perc = max(table(pre_anno)) / sum(table(pre_anno)))

# annotation should be in high concordance
# now try more annotations

mta <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/output/mta.rds")

ctomc <- left_join(ctomc, mta[, -2])

ctomc <- ctomc %>% group_by(metacell) %>% mutate(perc_leu1 = max(table(labels_leu1)) / sum(table(labels_leu1)))

ctomc <- ctomc %>% group_by(metacell) %>% mutate(perc_fbm = max(table(labels_fbm)) / sum(table(labels_fbm)))

to_join1 <- data.frame(metacell = as.numeric(rownames(annotation_col)), cluster = annotation_col$cluster)
ctomc <- left_join(ctomc, to_join1)

readr::write_rds(as.data.frame(ctomc), "/cluster/home/yjliu_jh/projects/mef2d/output/cells_to_metacells.rds")
readr::write_rds(umis, "/cluster/home/yjliu_jh/projects/mef2d/output/metacells_umis.rds")


# maybe use umap to backwards visualize the current annotation







# first check CD19 expression for metacells









