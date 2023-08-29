pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse",
          "CellChat", "ggplot2", "patchwork", "ggalluvial", "svglite", "Seurat")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}

work_dir <- c("/cluster/home/yjliu_jh/projects/mef2d/analysis/zwn/home/cellchat")
setwd(work_dir)

sce_tumor <- readr::read_rds("~/projects/mef2d/analysis/zwn/human/tenx/seurat/merge/merge.rds")
sce_hd <- readr::read_rds("~/projects/mef2d/analysis/zwn/human/tenx/seurat/merge/sce_hd.rds")
Idents(sce_hd) <- "cell_types"
sce_hd <- subset(x = sce_hd, subset = cell_types %notin% "remove")
dat_tumor <- sce_tumor@assays$RNA@data
meta_tumor <- data.frame(broad = sce_tumor$cell_types)
dat_hd <- sce_hd@assays$RNA@data
meta_hd <- data.frame(broad = sce_hd$cell_types)

# ====== overall cellchat object ======
cc_tumor <- createCellChat(object = dat_tumor, meta = meta_tumor, group.by = "broad")
cc_tumor@DB <- CellChatDB.human ## get DB data
cc_tumor <- subsetData(cc_tumor) ## filter for signaling genes
cc_hd <- createCellChat(object = dat_hd, meta = meta_hd, group.by = "broad")
cc_hd@DB <- CellChatDB.human ## get DB data
cc_hd <- subsetData(cc_hd) ## filter for signaling genes


# calculate cellchat
future::plan("multicore", workers = 20) ## do parallel
cc_tumor <- identifyOverExpressedGenes(cc_tumor)
cc_tumor <- identifyOverExpressedInteractions(cc_tumor)
cc_tumor <- projectData(cc_tumor, PPI.human)  ## PPI.human 4875x4875
cc_tumor <- computeCommunProb(cc_tumor)
cc_tumor <- filterCommunication(cc_tumor, min.cells = 10) ## filter groups under 10 cells

cc_hd <- identifyOverExpressedGenes(cc_hd)
cc_hd <- identifyOverExpressedInteractions(cc_hd)
cc_hd <- projectData(cc_hd, PPI.human)  ## PPI.human 4875x4875
cc_hd <- computeCommunProb(cc_hd)
cc_hd <- filterCommunication(cc_hd, min.cells = 10) ## filter groups under 10 cells

cc_tumor <- computeCommunProbPathway(cc_tumor)
cc_tumor <- aggregateNet(cc_tumor)
cc_hd <- computeCommunProbPathway(cc_hd)
cc_hd <- aggregateNet(cc_hd)

future::plan("sequential")

readr::write_rds(cc_tumor, "cellchat_tumor.rds")
readr::write_rds(cc_hd, "cellchat_healthy.rds")


# ====== separated cellchat object ======
Idents(sce_tumor) <- "orig.ident"
sce_TCF3 <- subset(x = sce_tumor, subset = orig.ident %in% "E1")
sce_MEF2D <- subset(x = sce_tumor, subset = orig.ident %notin% "E1")

dat_tcf <- sce_TCF3@assays$RNA@data
meta_tcf <- data.frame(broad = sce_TCF3$cell_types)
dat_mef <- sce_MEF2D@assays$RNA@data
meta_mef <- data.frame(broad = sce_MEF2D$cell_types)

cc_tcf <- createCellChat(object = dat_tcf, meta = meta_tcf, group.by = "broad")
cc_tcf@DB <- CellChatDB.human ## get DB data
cc_tcf <- subsetData(cc_tcf) ## filter for signaling genes
cc_mef <- createCellChat(object = dat_mef, meta = meta_mef, group.by = "broad")
cc_mef@DB <- CellChatDB.human ## get DB data
cc_mef <- subsetData(cc_mef) ## filter for signaling genes


# calculate cellchat
future::plan("multicore", workers = 20) ## do parallel
cc_tcf <- identifyOverExpressedGenes(cc_tcf)
cc_tcf <- identifyOverExpressedInteractions(cc_tcf)
cc_tcf <- projectData(cc_tcf, PPI.human)  ## PPI.human 4875x4875
cc_tcf <- computeCommunProb(cc_tcf)
cc_tcf <- filterCommunication(cc_tcf, min.cells = 10) ## filter groups under 10 cells

cc_mef <- identifyOverExpressedGenes(cc_mef)
cc_mef <- identifyOverExpressedInteractions(cc_mef)
cc_mef <- projectData(cc_mef, PPI.human)  ## PPI.human 4875x4875
cc_mef <- computeCommunProb(cc_mef)
cc_mef <- filterCommunication(cc_mef, min.cells = 10) ## filter groups under 10 cells

cc_tcf <- computeCommunProbPathway(cc_tcf)
cc_tcf <- aggregateNet(cc_tcf)
cc_mef <- computeCommunProbPathway(cc_mef)
cc_mef <- aggregateNet(cc_mef)

future::plan("sequential")

readr::write_rds(cc_tcf, "cellchat_TCF3.rds")
readr::write_rds(cc_mef, "cellchat_MEF2D.rds")



# ======= analysis part =========
# source("/cluster/home/yjliu_jh/projects/mef2d/code/common/lift.R")

# separate analysis

# === output plots functions ===
circle_out <- function(obj, dir, option = c("weight", "count")) {
  option <- match.arg(option)
  net_circle_dir <- jhtools::checkdir(dir)
  groupSize <- as.numeric(table(obj@idents))
  mat <- obj@net[[option]]
  for (i in 1:nrow(mat)) {
    pdf(glue::glue("{net_circle_dir}/{rownames(mat)[i]}_{option}_netVisual_circle.pdf"), width = 9, height = 9)
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 2.5,
                     arrow.size = 0.3, edge.weight.max = max(mat), title.name = rownames(mat)[i])
    dev.off()
    png(glue::glue("{net_circle_dir}/{rownames(mat)[i]}_{option}_netVisual_circle.png"), width = 9, height = 9, units = "in", res = 600)
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 2.5,
                     arrow.size = 0.3, edge.weight.max = max(mat), title.name = rownames(mat)[i])
    dev.off()
  }
}
# igraph 1.4 (stricter) will cause problems for this function  ## https://github.com/sqjin/CellChat/issues/547
# downgrade to igraph 1.35 or upgrade CellChat to latest


pathway_out <- function(obj, dir, option = c("chord", "circle")) {
  option <- match.arg(option)
  pathways_dir <- jhtools::checkdir(dir)
  for (i in 1:length(obj@netP$pathways)) {
    pdf(glue::glue("{pathways_dir}/{obj@netP$pathways[i]}_{option}.pdf"), width = 4.5, height = 6)
    netVisual_aggregate(obj, signaling = obj@netP$pathways[i], vertex.weight = 0.5, arrow.width = 2.5,
                        arrow.size = 0.3, layout = option)
    dev.off()
    png(glue::glue("{pathways_dir}/{obj@netP$pathways[i]}_{option}.png"), width = 4.5, height = 6, units = "in", res = 600)
    netVisual_aggregate(obj, signaling = obj@netP$pathways[i], vertex.weight = 0.5, arrow.width = 2.5,
                        arrow.size = 0.3, layout = option)
    dev.off()
  }
}


# pathway_heatmap <- function(obj, dir) {
#   pathways_dir <- jhtools::checkdir(dir)
#   for (i in 1:length(obj@netP$pathways)) {
#     pdf(glue::glue("{pathways_dir}/{obj@netP$pathways[i]}_heatmap.pdf"), width = 5, height = 5)
#     netVisual_heatmap(obj, signaling = obj@netP$pathways[i], color.heatmap=c("#2166ac", "#b2182b"))
#     dev.off()
#     png(glue::glue("{pathways_dir}/{obj@netP$pathways[i]}_heatmap.png"), width = 5, height = 5, units = "in", res = 600)
#     netVisual_heatmap(obj, signaling = obj@netP$pathways[i], color.heatmap=c("#2166ac", "#b2182b"))
#     dev.off()
#   }
# }
#   
  


# === plot ===
# update to newest cellchat object
cc_hd_up <- updateCellChat(cc_hd)
cc_mef_up <- updateCellChat(cc_mef)
cc_tcf_up <- updateCellChat(cc_tcf)


circle_out(cc_hd_up, glue::glue("{work_dir}/healthy/circle_weight"), "weight")
circle_out(cc_hd_up, glue::glue("{work_dir}/healthy/circle_count"), "count")
pathway_out(cc_hd_up, glue::glue("{work_dir}/healthy/pathway_chord"), "chord")
pathway_out(cc_hd_up, glue::glue("{work_dir}/healthy/pathway_circle"), "circle")

circle_out(cc_mef_up, glue::glue("{work_dir}/MEF2D/circle_weight"), "weight")
circle_out(cc_mef_up, glue::glue("{work_dir}/MEF2D/circle_count"), "count")
pathway_out(cc_mef_up, glue::glue("{work_dir}/MEF2D/pathway_chord"), "chord")
pathway_out(cc_mef_up, glue::glue("{work_dir}/MEF2D/pathway_circle"), "circle")


circle_out(cc_tcf_up, glue::glue("{work_dir}/TCF3/circle_weight"), "weight")
circle_out(cc_tcf_up, glue::glue("{work_dir}/TCF3/circle_count"), "count")
pathway_out(cc_tcf_up, glue::glue("{work_dir}/TCF3/pathway_chord"), "chord")
pathway_out(cc_tcf_up, glue::glue("{work_dir}/TCF3/pathway_circle"), "circle")












# plot and analyzations without loop
# overall:




# specific: 





# merge data
source("xxx.R")
cc_all <- liftCellChat_ADAPTED(cc_all, c("CLP", "erythroid_cell", "immatureB", "matureB",  
                                         "myeloid", "NK_T", "pre_proB", "preB_I", "preB_II", "proB"))



