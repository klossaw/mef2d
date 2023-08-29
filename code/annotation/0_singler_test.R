pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "cowplot",
          "Seurat", "anndata", "harmony", "SingleR", "scuttle")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}

sce <- readr::read_rds("/cluster/home/jhuang/projects/mef2d/analysis/zwn/human/tenx/seurat/merged/mef2d.rds")
to_classify <- sce[, sce$`orig.ident` %in% c("E1", "B069_Dx_CD19", "M2", "M3", "M4")]
data_dir <- "/cluster/home/yjliu_jh/projects/mef2d/data"

# read 
abm <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/data/public/rds/adultbm_cells_matrix_sub.rds")
abm_c <- t(abm)
pdata_abm <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/data/public/rds/adultbm_cells_metadata.rds")

cord <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/data/public/rds/cordblood_cells_matrix_sub.rds")
cord_c <- t(cord)
pdata_cord <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/data/public/rds/cordblood_cells_metadata.rds")





singler_abm <- SingleR(test = sce_clsf_norm, ref = sce_abm_norm, 
                       labels = sce_abm_norm$cell_label, de.method = "wilcox")


singler_cord <- SingleR(test = sce_clsf_norm, ref = sce_fbm_norm, 
                       labels = sce_fbm_norm$cell.labels, de.method = "wilcox")




# ======== new

leu_s <- readr::read_rds(glue::glue("{data_dir}/ref_se_aml_30827681_small.rds"))


comm_genes_leu_s <- intersect(rownames(to_classify[["RNA"]]@counts), rownames(leu_s@assays@data@listData$counts))
to_classify <- to_classify[comm_genes_leu_s, ]
sce_clsf <- SummarizedExperiment(assays = list(counts = to_classify[["RNA"]]@counts))
sce_clsf_norm <- scuttle::logNormCounts(sce_clsf)

leu_s <- leu_s[comm_genes_leu_s, ]
leu_s_norm <- scuttle::logNormCounts(leu_s)


singler_leu1_s <- SingleR(test = sce_clsf_norm, ref = leu_s_norm, 
                       labels = leu_s@colData@listData$CellType, de.method = "wilcox")

singler_leu1 <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/data/singleR_anno_leu1.rds")
comb <- data.frame(leu1 = singler_leu1$labels, leu1_s = singler_leu1_s$labels)
table(comb)





# ====== read existing singler annotations ======


singler_fbm <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/data/singleR_anno_fbm.rds")
singler_abm <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/data/singleR_anno_abm.rds")
singler_leu1 <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/data/singleR_anno_leu1.rds")
meta_tumor <- to_classify@meta.data
meta_tumor$cell_id <- rownames(meta_tumor)
singler_fbm_sub <- as.data.frame(singler_fbm[, c("delta.next", "pruned.labels")])
singler_abm_sub <- as.data.frame(singler_abm[, c("delta.next", "pruned.labels")])
singler_leu1_sub <- as.data.frame(singler_leu1[, c("delta.next", "pruned.labels")])
colnames(singler_fbm_sub) <- c("delta_fbm", "labels_fbm")
colnames(singler_abm_sub) <- c("delta_abm", "labels_abm")
colnames(singler_leu1_sub) <- c("delta_leu1", "labels_leu1")
singler_fbm_sub$cell_id <- rownames(singler_fbm_sub)
singler_abm_sub$cell_id <- rownames(singler_abm_sub)
singler_leu1_sub$cell_id <- rownames(singler_leu1_sub)

meta_tumor_all <- meta_tumor %>% left_join(singler_fbm_sub) %>% 
  left_join(singler_abm_sub) %>% left_join(singler_leu1_sub)
meta_tumor_all <- meta_tumor_all[, c("cell_id", "cell_types", "labels_fbm", "labels_abm", "labels_leu1")]



# test singler usage on small data
singler_leu1_s <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/data/singleR_anno_leu1_s.rds")
singler_leu1_s_sub <- as.data.frame(singler_leu1_s[, c("delta.next", "pruned.labels")])
colnames(singler_leu1_s_sub) <- c("delta_leu1s", "labels_leu1s")
singler_leu1_s_sub$cell_id <- rownames(singler_leu1_s_sub)
singler_leu1_s_sub <- singler_leu1_s_sub[, c("cell_id", "labels_leu1s")]
mta <- left_join(meta_tumor_all, singler_leu1_s_sub)



# ====== based on the fbm annotation ======

fbm_summary <- meta_tumor_all[, c("labels_fbm", "cell_types")]
fbm_table <- table(fbm_summary)
fbm_short <- c("other", "other", "other", "other", "other", "other", "other", 
              "other", "other", "other", "ELP", "other", "other", "other", 
              "other", "immature B cell", "other", "other", "other", "other", "other", 
              "other", "other", "monocytoid macrophage", "other", "other", "other", 
              "other", "naive B cell", "other", "other", "other", "other", 
              "pre B progenitor", "pre pro B progenitor", "pro B progenitor",
              "other", "other", "other", "other", "other")
fbm_join <- data.frame(labels_fbm = rownames(fbm_table), group = fbm_short)
fbm_table <- as.data.frame(fbm_table[, c("tumor_cell1", "tumor_cell2", "tumor_cell3", "tumor1")])
fbm_table <- left_join(fbm_table, fbm_join)
fbm_table_s <- fbm_table %>% group_by(group, cell_types) %>% summarise(Freq = sum(Freq)) %>%
  group_by(cell_types) %>% mutate(prop =  100 * Freq / sum(Freq)) %>% mutate(ypos = 100 - cumsum(prop) + 0.5 * prop)

readr::write_rds(fbm_table_s, "/cluster/home/yjliu_jh/projects/mef2d/output/fbm_anno_freq.rds")
fbm_table_s <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/output/fbm_anno_freq.rds")


tb_tmcell1 <- fbm_table_s[fbm_table_s$cell_types %in% "tumor_cell1", ]
tb_tmcell2 <- fbm_table_s[fbm_table_s$cell_types %in% "tumor_cell2", ]
tb_tmcell3 <- fbm_table_s[fbm_table_s$cell_types %in% "tumor_cell3", ]

# quick function with fixed parameters and column names
quick_pie <- function(data){
  p <- ggplot(data, aes(x = "", y = prop, fill = group)) +
    geom_bar(stat = "identity", width = 10, size = 2, color = "white") +
    coord_polar("y", start = 0) +
    theme_void() + 
    theme(legend.position = "none") +
    geom_text(aes(y = ypos, label = paste(group, round(prop, 2), sep = "\n")), 
              color = "black", size = 3, nudge_x = 3) +
    scale_fill_brewer(palette = "Set2")
  p
}

p_tmcell1 <- quick_pie(tb_tmcell1)
p_tmcell2 <- quick_pie(tb_tmcell2)
p_tmcell3 <- quick_pie(tb_tmcell3)

p_all <- p_tmcell1 + p_tmcell2 + p_tmcell3

# different composition
# echoed the KMT2A publication (ELP)




# ========== based on the leu1 annotation ======






