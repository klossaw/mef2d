# construct merged sce object
sce_tumor <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/data/merged_4_tumors_new.rds")
sce_hd <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/data/healthy_donors.rds")
sce_all <- readr::read_rds("/cluster/home/jhuang/projects/mef2d/analysis/zwn/human/tenx/seurat/merged/all_merged.rds")
new_cell_id <- str_sub(colnames(sce_all), 1 + str_locate(colnames(sce_all), ":")[, 1] / 2)
sce_all <- RenameCells(sce_all, new.names = new_cell_id)

# set basic 
case = c("E1", "M2", "M3", "M4")
control = c("H1_CD19", "H1_CD34pCD19n", "H2_CD19", "H2_CD34pCD19n")

tm_id <- colnames(sce_tumor)
hd_id <- colnames(sce_hd)
sce_all_sub <- sce_all[, c(tm_id, hd_id)]

sce <- sce_all_sub
sce@meta.data$celltype <- "other"
sce$celltype[sce$orig.ident %in% case] <- sce_tumor$cell_types_broad
sce$orig.ident[sce$orig.ident %in% "Healthy"] <- sce_hd$orig.ident
sce$celltype[sce$orig.ident %in% control] <- sce_hd$cell_types
sce@meta.data$type <- "normal"
sce@meta.data$type[sce$orig.ident %in% c("M2", "M3", "M4")] <- "MEF2D"
sce@meta.data$type[sce$orig.ident %in% "E1"] <- "TCF3"
sce@meta.data$sample <- sce@meta.data$orig.ident
sce@meta.data$cell_id <- colnames(sce)

readr::write_rds(sce, "/cluster/home/yjliu_jh/projects/mef2d/data/mef_hd_all_cells.rds")





