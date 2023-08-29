pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr", "igraph",
          "Seurat", "anndata", "stats", "pheatmap", "pracma", "SeuratDisk"
)
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}



infercnv_obj <- readr::read_rds("/cluster/home/ylxie_jh/projects/leukemia/analysis/weinazhang/human/split_hd/infercnv/result/run.final.infercnv_obj")
cnv_mat <- infercnv_obj@expr.data
cnv_flat <- as.vector(cnv_mat)



test_func1 <- function(x){sum((x - median(x))^2)}
cnv1 <- apply(cnv_mat, 2, test_func1)

test_func2 <- function(x){sum((x - 1)^2)}
cnv2 <- apply(cnv_mat, 2, test_func2)

# scale to [-1, 1], zero  and calculate sum of squares sth.?
# row scale? 

test_func3 <- function(x){
  x <- x - median(x)  ## to avoid zero as denominator
  x[x < 0] <- x / min(x)
  x[x > 0] <- x / max(x)
  x
}
cnv_mat_scaled1 <- apply(cnv_mat, 1, test_func3)

cnv3 <- apply(cnv_mat_scaled1, 1, \(x) sum(x^2))

test_func4 <- function(x){
  x <- x - 1  ## to avoid zero as denominator
  x[x < 0] <- x / min(x)
  x[x > 0] <- x / max(x)
  x
}
cnv_mat_scaled2 <- apply(cnv_mat, 1, test_func4)

# these are not right



cnv4 <- apply(cnv_mat, 2, test_func4)

temp = data.frame(cell_id = names(cnv1), cnv1, cnv2)
readr::write_rds(temp, "/cluster/home/yjliu_jh/projects/mef2d/data/cnvscores.rds")



temp <- sce_tumor@meta.data
temp <- left_join(temp, data.frame(metacell = as.numeric(rownames(annotation_col3)), clustx = annotation_col3$clustx))


annoc <- c("matureB1", "myeloid", "CLP", "naive_proB", "naive_proB2", "erythoid", 
  "proB_1", "proB_2", "NK", "T", "proB_like", "NK-T", 
  "ATXN1_high_BCL11Bn_preB", "matureB2", "proB_3", "preB", "preB_IGKCp", "ery-like")
annoc <- data.frame(clustx = 1:18, cluster_anno = annoc)
temp <- left_join(temp, annoc)
sce_tumor@meta.data$cluster_anno <- temp$cluster_anno

annoc2 <- c("matureB", "myeloid", "CLP", "naive_proB", "proB1", "naive_proB2", 
            "proB2", "erythoid", "proB3", "proB4", "NK-T", "proB5", 
            "proB6", "preB1", "preB2", "preB-like", "preB", "ery-like")
sannoc2 <- data.frame(cluster_final = 1:18, cluster_final_anno = annoc2)
temp <- sce_tumor@meta.data
temp <- left_join(temp, sannoc2)
sce_tumor@meta.data$cluster_final_anno <- temp$cluster_final_anno

sce_tumor@meta.data[sce_tumor@meta.data$cluster_final_anno %in% "T", "cluster_final_anno"] <- sce_tumor@meta.data[sce_tumor@meta.data$cluster_final_anno %in% "T", "cluster_anno"]





