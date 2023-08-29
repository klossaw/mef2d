pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "cowplot",
          "Seurat", "vroom", "loomR", "anndata", "harmony", "SingleR", "scuttle")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}

project <- "project"
dataset <- "dataset"
species <- "human"
workdir <- glue("~/projects/{project}/analysis/{dataset}/{species}")
setwd(workdir)


sce <- readr::read_rds("~/projects/mef2d/analysis/zwn/human/tenx/seurat/merged/mef2d.rds")
sce@meta.data$group <- ifelse(sce@meta.data$orig.ident %in% c("B069_Dx_CD19", "M2", "M3", "M4"), "MEF2D", "others")
Idents(object = sce) <- "group"

# some random differential tests


sce_tumor <- subset(sce, `orig.ident` %in% c("B069_Dx_CD19", "M2", "M3", "M4", "E1"))


# first locate differentially expressed genes between the 2 groups
diff_markers <- FindMarkers(sce, ident.1 = "MEF2D", ident.2 = "others", method = "DESeq2")
diff_markers$perc_diff <- diff_markers$pct.1 - diff_markers$pct.2
head(diff_markers[order(diff_markers$perc_diff), ], 10)

# some TCF4 MAN1A1 MEF2C HLA-DQA1 ACSM3 ...
# CD37 is low (not so high) expressed in MEF2D samples!


# then, for the marker genes, find the co-expression pattern





# TCF4 might play in concert with CREB, MEF2 and other transcription factors 
# to modulate BDNF levels following neuronal activity










# ===== suddenly: =======
# read infercnv result files
library(phylogram)
library(dendextend)
ic_dend <- read.dendrogram("/cluster/home/ylxie_jh/projects/leukemia/analysis/weinazhang/human/split_hd/infercnv/result/infercnv.observations_dendrogram.txt")

# Cut tree 
ic_labels <- cutree(ic_dend, k = 10)
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


# ========== find and interpret differences  ======

to_classify <- sce[, sce$`orig.ident` %in% c("E1", "B069_Dx_CD19", "M2", "M3", "M4")]
meta_tumor <- to_classify@meta.data
rownames(to_classify@meta.data) <- sub(":", "-", rownames(to_classify@meta.data))
to_classify@meta.data$group <- ifelse(rownames(to_classify@meta.data) %in% ocells, "CN_altered", "others")
Idents(object = to_classify) <- "group"

# some random differential tests




# check gene expression difference between the two groups
diff_markers_cnv <- FindMarkers(to_classify, ident.1 = "CN_altered", ident.2 = "others", method = "DESeq2")
diff_markers_cnv$perc_diff <- diff_markers_cnv$pct.1 - diff_markers_cnv$pct.2
diff_cn_genes <- rownames(diff_markers_cnv[order(diff_markers_cnv$perc_diff), ])

library(clusterProfiler)
library(org.Hs.eg.db)

check_ego_gene <- function(data){
  ego <- enrichGO(gene = data, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                  ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
  ego <- ego[order(ego$p.adjust), ]
  ego
}

check_ego_gene(diff_cn_genes[1:50])


# check differences within cancer group
to_classify@meta.data$group2 <- ifelse(to_classify@meta.data$orig.ident %in% "E1", "others", "MEF2D")
Idents(object = to_classify) <- "group2"
diff_markers_it <- FindMarkers(to_classify, ident.1 = "MEF2D", ident.2 = "others", method = "DESeq2")
diff_markers_it$perc_diff <- diff_markers_it$pct.1 - diff_markers_it$pct.2
dm_ordered_it <- diff_markers_it[order(diff_markers_it$perc_diff), ]

# ataxin-1 binds specifically to histone deacetylase-4 (HDAC4) and MEF2


# CD44 low expression? 






















# ======= cell annotations =======
# temp scripts
new_anno <- mef2d_meta_new
new_anno$singler3 <- singler_leu1$pruned.labels
readr::write_rds(new_anno, "/cluster/home/yjliu_jh/projects/mef2d/output/new_cell_anno.rds")
cor_mat3 <- table(new_anno[, c("cell_types", "singler3")])

# seems alright, construct a plot function haoleya

scaled_corr_map <- function(mat){
  mat_row_perc_scale <- mat / rowSums(mat)
  ht <- ComplexHeatmap::Heatmap(mat_row_perc_scale, cluster_columns = F, cluster_rows = F)
  ht
}


temp2 <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/output/new_cell_anno.rds")

scaled_corr_map(table(temp2[, c("cell_types", "singler3")]))  ## color needs adjustment 










