pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr", "igraph",
          "Seurat", "anndata", "stats", "pheatmap", "pracma", "SeuratDisk"
)
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}

mdata <- anndata::read_h5ad('/cluster/home/yjliu_jh/projects/mef2d/data/metacells.h5ad')


temp <- read_h5ad("/cluster/home/yjliu_jh/projects/mef2d/data/public/fig1b_fbm_scaled_gex_updated_dr_20210104.h5ad")

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
colors <- colorRampPalette(c('darkred', 'red', 'white', 'white', 'lightblue', 'blue', 'darkblue'))(141)
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

pdf("/cluster/home/yjliu_jh/projects/mef2d/output/test_heatmap.pdf", width = 15, height = 12)
pht
dev.off()




# ===== use BALL genes to draw heatmap ======

mdata$obs$cluster <- as.character(cutree(pht$tree_col, k = 7))
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
cluster_colors <- c("#ee5d5a", "#34956c", "#6b6498", "#71a3a2",
                             "#c88978", "#cc9b32", "#53096a")
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

pdf("/cluster/home/yjliu_jh/projects/mef2d/output/test_heatmap2.pdf", width = 15, height = 18)
pht2
dev.off()


# remove genes with no variation
rowclust = hclust(dist(pht2data))
reordered = pht2data[rowclust$order, ]
fil_genes <- rownames(reordered)[c(1:25, 79:94)]


# ======= 3rd round heatmap ======

ball_genes_fil2 <- ball_genes_fil[ball_genes_fil$gene %in% fil_genes, ]

annotation_row3 <- data.frame(ball=ball_genes_fil2$celltype)
rownames(annotation_row3) <- ball_genes_fil2$gene

pht3data <- t(relative_log_fractions[, ball_genes_fil2$gene])

pht3 <- pheatmap::pheatmap(
  pht3data,
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
  annotation_col=annotation_col,
  annotation_colors=annotation_colors2,
  annotation_legend=TRUE,
  legend=TRUE
)

pdf("/cluster/home/yjliu_jh/projects/mef2d/output/test_heatmap3.pdf", width = 15, height = 8)
pht3
dev.off()



# ========

pht4 <- pheatmap::pheatmap(
  pht3data,
  cluster_rows = F,
  treeheight_col=0,
  treeheight_row=0,
  cellwidth=1,
  cellheight=10,
  show_rownames=TRUE,
  show_colnames=FALSE,     t
  main='Metacell Clusters - BALL anno Genes',
  color=colors,
  breaks=breaks,
  annotation_row=annotation_row3,
  annotation_col=annotation_col,
  annotation_colors=annotation_colors2,
  annotation_legend=TRUE,
  legend=TRUE
)

pdf("/cluster/home/yjliu_jh/projects/mef2d/output/test_heatmap4.pdf", width = 15, height = 8)
pht4
dev.off()















summary(mdata$obs$grouped)
table(mdata$obs$pile)# pile indicates parallel calculation, not related to groups





gg <- ggpubr::ggscatter(mdata$obs, "umap_x", "umap_y", color = "cluster", size = 0.4)

pdf("/cluster/home/yjliu_jh/projects/mef2d/output/test_umap.pdf", width = 8, height = 8)
gg
dev.off()

#ggsave("/cluster/home/yjliu_jh/projects/mef2d/output/test_umap.png", gg,
#       width = 800, height = 800, units = "px")



readr::write_rds(fbm_mat, "/cluster/home/yjliu_jh/projects/mef2d/data/h5ad1_mat_dgc.rds")
readr::write_rds(h5ad1$var, "/cluster/home/yjliu_jh/projects/mef2d/data/h5ad1_var.rds")
readr::write_rds(h5ad1$obs, "/cluster/home/yjliu_jh/projects/mef2d/data/h5ad1_obs.rds")



cells_temp <- anndata::read_h5ad('/cluster/home/yjliu_jh/projects/mef2d/data/cells.h5ad')
ctomc <- cells_temp$obs
ctomc$cell_id <- rownames(ctomc)
ctomc <- ctomc %>% group_by(metacell) %>% mutate(perc = max(table(pre_anno)) / sum(table(pre_anno)))

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








