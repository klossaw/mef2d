pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr", "igraph",
          "Seurat", "Asgard", "celldex", "cmapR", "pracma", "SeuratDisk",
          "cowplot", "edgeR"
)
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}


sce_tumor <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/data/merged_4_tumors.rds")
sce_hd <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/data/healthy_donors.rds")
sce_all <- readr::read_rds("/cluster/home/jhuang/projects/mef2d/analysis/zwn/human/tenx/seurat/merged/all_merged.rds")
new_cell_id <- str_sub(colnames(sce_all), 1 + str_locate(colnames(sce_all), ":")[, 1] / 2)
sce_all <- RenameCells(sce_all, new.names = new_cell_id)

tm_id <- colnames(sce_tumor)
hd_id <- colnames(sce_hd)
sce_all_fil <- sce_all[, c(tm_id, hd_id)]
sce_all_fil@meta.data$celltype <- "other"
case = c("E1", "M2", "M3", "M4")
control = c("H1_CD19", "H1_CD34pCD19n", "H2_CD19", "H2_CD34pCD19n")
sce_all_fil$celltype[sce_all_fil$orig.ident %in% case] <- sce_tumor$cell_types_broad[sce_tumor$orig.ident %in% case]
sce_all_fil$orig.ident[sce_all_fil$orig.ident %in% "Healthy"] <- sce_hd$orig.ident
sce_all_fil$celltype[sce_all_fil$orig.ident %in% control] <- sce_hd$cell_types[sce_hd$orig.ident %in% control]


sce_all_fil <- sce_all_fil %>% RunPCA() %>% FindNeighbors() %>% RunUMAP(dims = 1:15)

readr::write_rds(sce_all_fil, "/cluster/home/yjliu_jh/share/mef_sce_all.rds")
sce_all_fil <- readr::read_rds("/cluster/home/yjliu_jh/share/mef_sce_all.rds")

sce_all_fil@meta.data$type <- ifelse(sce_all_fil@meta.data$orig.ident %in% c("E1", "M2", "M3", "M4"), "tumor", "normal")

sce_all_fil <- sce_all_fil %>% RunUMAP(dims = 1:15)

Idents(sce_all_fil) <- "celltype"
DimPlot(sce_all_fil, reduction = "umap", group.by = "celltype")

sce_all_fil <- sce_all_fil %>% RunHarmony("type") %>%
  RunUMAP(dims = 1:10, reduction.name = "umap10") %>%
  RunUMAP(dims = 1:20, reduction.name = "umap20") %>% #ScaleData() %>%
  RunUMAP(reduction = "harmony", dims = 1:10, reduction.name = "umap10_h") %>%
  RunUMAP(reduction = "harmony", dims = 1:20, reduction.name = "umap20_h") %>% 
 # RunUMAP(reduction = "harmony", features = gene01, reduction.name = "umap01_h") %>%
  RunUMAP(reduction = "harmony", features = gene03, reduction.name = "umap03_h") %>%
  RunUMAP(reduction = "harmony", features = gene003y, reduction.name = "umap003_h") %>% 
  RunUMAP(reduction = "harmony", features = gene02, reduction.name = "umap02_h")


DimPlot(sce_all_fil, reduction = "umap02_h", group.by = "celltype")



DimPlot(sce_tumor, reduction = "umap", group.by = "seurat_clusters", label = T)

