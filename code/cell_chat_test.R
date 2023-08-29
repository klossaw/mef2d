pkgs <- c("CellChat", "ggplot2", "patchwork", "ggalluvial", "svglite", "Seurat")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}

work_dir <- c('/cluster/home/yjliu_jh/projects/mef2d/analysis/zwn/home/cellchat/')
project_dir <- file.path(work_dir,"tumor_cell_types")
setwd(project_dir)

sce_tumor2 <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/data/sce_tumor_final.rds")

# first, overall cellchat
dat_chat <- sce_tumor_fil@assays$RNA@data
meta_chat <- data.frame(group = sce_tumor_fil$cluster_final_anno)
cellchat <- createCellChat(object = dat_chat, meta = meta_chat, group.by = "group")
cellchat@DB <- CellChatDB.human ## get DB data
cellchat <- subsetData(cellchat) ## filter for signaling genes

future::plan("multicore", workers = 20) ## do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)  ## PPI.human 4875x4875
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10) ## filter groups under 10 cells


# future::plan("sequential")
#access the inferred cell-cell communications of interest
comm_net <- subsetCommunication(cellchat)

cluster_colors2 <- c("#ee5d5a", "#34956c", "#6b6498", "#71a3a2", "#c88978", "#cc9b32", "#53096a", 
                              RColorBrewer::brewer.pal(7, "Set2"), RColorBrewer::brewer.pal(6, "Pastel1"))
names(cluster_colors2) <- unique(sce_tumor_fil@meta.data$cluster_final_anno)

jhtools::df2excel(comm_net, "inferred_cell-cell_communications.xlsx")
#readr::write_rds(cellchat, "cellchat0.rds")


#Infer the cell-cell communication at a signaling pathway level
#We can calculate the aggregated cell-cell communication network 
#by counting the number of links or summarizing the communication probability.
cellchat <- computeCommunProbPathway(cellchat)
#Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)
#showing the number of interactions or the total interaction strength (weights) 
#between any two cell groups using circle plot
#readr::write_rds(cellchat, file = "aggregateNet_cellchat.rds")


#cellchat <- read_rds("aggregateNet_cellchat.rds")
if(F){
  net_circle_dir <- jhtools::checkdir(glue::glue("{project_dir}/cellTypes_net_circle"))
  groupSize <- as.numeric(table(cellchat@idents))
  interaction_types <- c("cellchat@net$count", "cellchat@net$weight")
  grDevices::cairo_pdf(glue::glue("{net_circle_dir}/all_netVisual_circle.pdf"), width = 18, height = 9)
  par(mfrow=c(1,2), xpd=TRUE
      #, mar = c(5, 4, 4, 2) +0.1
  )
  for (i in 1:2) {
    netVisual_circle(interaction_types[i], vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions", arrow.width = 2.5, arrow.size = 0.3)
  }
  dev.off()
}

#showed just one fig
groupSize <- as.numeric(table(cellchat@idents))
#par(mfrow = c(1,2), xpd=TRUE)
net_circle_dir <- jhtools::checkdir(glue::glue("{project_dir}/cellTypes_net_circle"))
pdf(glue::glue("{net_circle_dir}/all_netVisual_circle_Number.pdf"), width = 9, height = 9)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions", arrow.width = 2.5, arrow.size = 0.3)
dev.off()
pdf(glue::glue("{net_circle_dir}/all_netVisual_circle_weights.pdf"), width = 9, height = 9)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength", arrow.width = 2.5, arrow.size = 0.3)
dev.off()

#Due to the complicated cell-cell communication network, we can examine `the signaling sent from each cell group`. 
#Here we also control the parameter edge.weight.max so that we can compare edge weights between differet networks.
mat <- cellchat@net$weight
#par(mar=c(1,1,1.5,1))
for (i in 1:nrow(mat)) {
  pdf(glue::glue("{net_circle_dir}/{rownames(mat)[i]}_weight_netVisual_circle.pdf"), width = 9, height = 9)
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 2.5,
                   arrow.size = 0.3, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  dev.off()
  png(glue::glue("{net_circle_dir}/{rownames(mat)[i]}_weight_netVisual_circl.png"), width = 9, height = 9, units = "in", res = 600)
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 2.5,
                   arrow.size = 0.3, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  dev.off()
}

mat <- cellchat@net$count
#par(mar=c(1,1,1.5,1))
for (i in 1:nrow(mat)) {
  pdf(glue::glue("{net_circle_dir}/{rownames(mat)[i]}_count_netVisual_circle.pdf"), width = 9, height = 9)
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ] 
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 2.5,
                   arrow.size = 0.3, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  dev.off()
  png(glue::glue("{net_circle_dir}/{rownames(mat)[i]}_count_netVisual_circle.png"), width = 9, height = 9, units = "in", res = 600)
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 2.5,
                   arrow.size = 0.3, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  dev.off()
}

#Visualization of cell-cell communication network
#Visualize each signaling pathway using Hierarchy plot, Circle plot or Chord diagram
#All the signaling pathways showing significant communications can be accessed by 
cellchat@netP$pathways
#pathways.show <- c("VEGF") 
if(F){
  # Hierarchy plot
  levels(cellchat@idents) #check the order of the cell cluster
  # Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
  vertex.receiver = seq(4,7) # a numeric vector. 指定靶细胞的索引
  #netVisual_aggregate(cellchat, signaling = pathways.show, vertex.weight = 0.5, vertex.receiver = vertex.receiver,layout="hierarchy")
  for (i in 1:length(cellchat@netP$pathways)) {
    pdf(glue::glue("{cellchat@netP$pathways[i]}_hierarchy.pdf"), width = 14, height = 9)
    netVisual_aggregate(cellchat, signaling = cellchat@netP$pathways[i], vertex.weight = 0.5,  vertex.receiver = vertex.receiver,layout="hierarchy")
    dev.off()
    png(glue::glue("{cellchat@netP$pathways[i]}_hierarchy.png"), width = 14, height = 9, units = "in", res = 600)
    netVisual_aggregate(cellchat, signaling = cellchat@netP$pathways[i], vertex.weight = 0.5,  vertex.receiver = vertex.receiver,layout="hierarchy")
    dev.off()
  }
}

pathways_dir <- jhtools::checkdir(glue::glue("{project_dir}/pathways"))
# Circle plot
#par(mfrow=c(1,1))
#netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
for (i in 1:length(cellchat@netP$pathways)) {
  pdf(glue::glue("{pathways_dir}/{cellchat@netP$pathways[i]}_circle.pdf"), width = 4.5, height = 6)
  netVisual_aggregate(cellchat, signaling = cellchat@netP$pathways[i], vertex.weight = 0.5, arrow.width = 2.5,
                      arrow.size = 0.3, layout="circle")
  dev.off()
  png(glue::glue("{pathways_dir}/{cellchat@netP$pathways[i]}_circle.png"), width = 4.5, height = 6, units = "in", res = 600)
  netVisual_aggregate(cellchat, signaling = cellchat@netP$pathways[i], vertex.weight = 0.5, arrow.width = 2.5,
                      arrow.size = 0.3, layout="circle")
  dev.off()
}

# Chord diagram
#par(mfrow=c(1,1))
#netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
for (i in 1:length(cellchat@netP$pathways)) {
  pdf(glue::glue("{pathways_dir}/{cellchat@netP$pathways[i]}_chord.pdf"), width = 5, height = 5)
  netVisual_aggregate(cellchat, signaling = cellchat@netP$pathways[i], vertex.weight = 0.5, layout="chord")
  dev.off()
  png(glue::glue("{pathways_dir}/{cellchat@netP$pathways[i]}_chord.png"), width = 5, height = 5, units = "in", res = 600)
  netVisual_aggregate(cellchat, signaling = cellchat@netP$pathways[i], vertex.weight = 0.5, layout="chord")
  dev.off()
}

# Heatmap
#par(mfrow=c(1,1))
#netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

for (i in 1:length(cellchat@netP$pathways)) {
  pdf(glue::glue("{pathways_dir}/{cellchat@netP$pathways[i]}_heatmap.pdf"), width = 5, height = 5)
  netVisual_heatmap(cellchat, signaling = cellchat@netP$pathways[i], color.heatmap=c("#2166ac", "#b2182b"))
  dev.off()
  png(glue::glue("{pathways_dir}/{cellchat@netP$pathways[i]}_heatmap.png"), width = 5, height = 5, units = "in", res = 600)
  netVisual_heatmap(cellchat, signaling = cellchat@netP$pathways[i], color.heatmap=c("#2166ac", "#b2182b"))
  dev.off()
}

if(F){
  # Chord diagram
  table(cellchat@idents)
  group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # exp:grouping cell clusters into fibroblast, DC and TC cells (more general cell groups
  names(group.cellType) <- levels(cellchat@idents)
  netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))
  #> Plot the aggregated cell-cell communication network at the signaling pathway level
}


# find LR pairs with high contribution for pathways and visualize
# pathways.show <- c("VEGF") 
netAnalysis_contribution(cellchat, signaling = pathways.show)

#可视化由单个配体受体对调节的细胞-细胞通信.
#函数extractEnrichedLR来提取给定信号通路的所有重要相互作用（L-R对）和相关信号基因。
pairLR <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR[1,] # show one ligand-receptor pair

if(F){
  # Hierarchy plot
  vertex.receiver = seq(1,4) # a numeric vector
  netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
  netVisual(cellchat, signaling = pathways.show, vertex.receiver = vertex.receiver, layout = "hierarchy")
}

# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
# Chord diagram
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.

# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
for (i in 1:length(pathways.show.all)) {
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}







