pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse",
          "CellChat", "ggplot2", "patchwork", "ggalluvial", "svglite", "Seurat")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}

work_dir <- c('/cluster/home/yjliu_jh/projects/mef2d/analysis/zwn/home/cellchat/')
project_dir <- file.path(work_dir,"tumor_cell_types")
setwd(project_dir)

sce_tumor2 <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/data/sce_tumor_final.rds")
sce_tumor <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/data/merged_4_tumors.rds")

#Idents(sce_tumor_fil) <- "cell_types_all"
Idents(sce_tumor2) <- "orig.ident"
sce_TCF3 <- subset(x = sce_tumor2, subset = orig.ident %in% "E1")
sce_MEF2D <- subset(x = sce_tumor2, subset = orig.ident %notin% "E1")









# after naming the clusters
# save plots both from TCF3 and MEF2D to one 
# first, overall cellchat
dat_tcf <- sce_TCF3@assays$RNA@data
meta_tcf <- data.frame(all = sce_TCF3$cell_types_all, broad = sce_TCF3$cell_types_broad)
dat_mef <- sce_MEF2D@assays$RNA@data
meta_mef <- data.frame(all = sce_MEF2D$cell_types_all, broad = sce_MEF2D$cell_types_broad)

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

readr::write_rds(cc_tcf, "/cluster/home/yjliu_jh/projects/mef2d/output/cellchat_TCF3_1.rds")
readr::write_rds(cc_mef, "/cluster/home/yjliu_jh/projects/mef2d/output/cellchat_MEF2D_1.rds")


comm_net_tcf <- subsetCommunication(cc_tcf)
comm_net_mef <- subsetCommunication(cc_mef)
jhtools::df2excel(comm_net_mef, "inferred_cell_communications_mef2d.xlsx")
jhtools::df2excel(comm_net_mef, "inferred_cell_communications_tcf3.xlsx")



cc_all <- mergeCellChat(list(cc_tcf, cc_mef), add.names = c("TCF3", "MEF2D")) 
readr::write_rds(cc_all, "/cluster/home/yjliu_jh/projects/mef2d/output/cellchat_all_1.rds")
cc_all <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/output/cellchat_all_1.rds")


compareInteractions(cc_all, show.legend = F, group = c(1, 2))
compareInteractions(cc_all, show.legend = F, group = c(1, 2), measure = "weight")



liftCellChat_ADAPTED <- function(object, group.new = NULL) {
  if (object@options$mode == "merged") {
    idents <- object@idents[1:(length(object@idents)-1)]
    if (is.null(group.new)) {
      group.max.all <- unique(unlist(sapply(idents, levels)))
      group.num <- sapply(idents, nlevels)
      group.num.max <- max(group.num)
      group.max <- levels(idents[[which(group.num == group.num.max)]])
      if (length(group.max) != length(group.max.all)) {
        stop("CellChat object cannot lift up due to the missing cell groups in any dataset. Please define the parameter `group.new`!")
      }
    } else {
      group.max <- group.new
      group.num.max <- length(group.new)
    }
    message(paste0("The CellChat object will be lifted up using the cell labels ", paste(group.max, collapse=", ")))
    for (i in 1:length(idents)) {
      cat("Update slots object@net, object@netP, object@idents in dataset ", names(object@idents)[i],'\n')
      # cat("Update slot object@net...", '\n')
      net <- object@net[[i]]
      group.i <- levels(idents[[i]])
      group.existing <- group.max[group.max %in% group.i]
      # CHANGED THIS LINE: in order to get the indices of the cell types of group.new that exist in the orignal levels we changed it to this: 
      group.existing.index <- unlist(lapply(group.i, function(x) which(group.max %in% x)))
      for (net.j in names(net)) {
        values <- net[[net.j]]
        if (net.j %in% c("prob","pval")) {
          values.new <- array(data = 0, dim = c(group.num.max, group.num.max, dim(values)[3]),
                              dimnames = list(group.max, group.max, dimnames(values)[[3]]))
          values.new[group.existing.index, group.existing.index, ] <- values
          net[[net.j]] <- values.new
        }
        if (net.j %in% c("count","sum","weight")) {
          values.new <- array(data = 0, dim = c(group.num.max, group.num.max),
                              dimnames = list(group.max, group.max))
          values.new[group.existing.index, group.existing.index] <- values
          net[[net.j]] <- values.new
        }
        if (net.j %in% c("pairwiseRank")) {
          for (k in 1:length(values)) {
            values.new1 <- vector("list", group.num.max)
            values.new1[group.existing.index] <- values[[k]]
            temp <- values[[k]][[1]]
            temp$prob  <- 0; temp$pval <- 1
            for (kk in setdiff(1:group.num.max, group.existing.index)) {
              values.new1[[kk]] <- temp
            }
            names(values.new1) <- group.max
            values[[k]] <- values.new1
          }
          values.new <- vector("list", group.num.max)
          values.new[group.existing.index] <- values
          temp <- lapply(values.new1, function(x) {
            x$prob  <- 0; x$pval <- 1
            return(x)
          })
          for (kk in setdiff(1:group.num.max, group.existing.index)) {
            values.new[[kk]] <- temp
          }
          names(values.new) <- group.max
        }
        net[[net.j]] <- values.new
      }
      object@net[[i]] <- net
      
      # cat("Update slot object@netP...", '\n')
      netP <- object@netP[[i]]
      for (netP.j in names(netP)) {
        values <- netP[[netP.j]]
        if (netP.j %in% c("pathways")) {
          values.new <- values
          netP[[netP.j]] <- values.new
        }
        if (netP.j %in% c("prob")) {
          values.new <- array(data = 0, dim = c(group.num.max, group.num.max, dim(values)[3]),
                              dimnames = list(group.max, group.max, dimnames(values)[[3]]))
          values.new[group.existing.index, group.existing.index, ] <- values
          netP[[netP.j]] <- values.new
        }
        if (netP.j %in% c("centr")) {
          for (k in 1:length(values)) {
            values.new <- lapply(values, function(x) {
              values.new2 <- lapply(x, function(x) {
                values.new1 = as.vector(matrix(0, nrow = 1, ncol = group.num.max))
                values.new1[group.existing.index] <- x
                names(values.new1) <- group.max
                return(values.new1)
              })
              names(values.new2) <- names(x)
              return(values.new2)
            })
            names(values.new) <- names(values)
          }
        }
        netP[[netP.j]] <- values.new
      }
      object@netP[[i]] <- netP
      # cat("Update slot object@idents...", '\n')
      # idents[[i]] <- factor(group.max, levels = group.max)
      idents[[i]] <- factor(idents[[i]], levels = group.max)
    }
    object@idents[1:(length(object@idents)-1)] <- idents
  } else {
    if (is.null(group.new)) {
      stop("Please define the parameter `group.new`!")
    } else {
      group.max <- as.character(group.new)
      group.num.max <- length(group.new)
      message(paste0("The CellChat object will be lifted up using the cell labels ", paste(group.max, collapse=", ")))
    }
    cat("Update slots object@net, object@netP, object@idents in a single dataset...", '\n')
    # cat("Update slot object@net...", '\n')
    net <- object@net
    idents <- object@idents
    group.i <- levels(idents)
    group.existing <- group.max[group.max %in% group.i]
    # CHANGED THIS LINE: in order to get the indices of the cell types of group.new that exist in the orignal levels we changed it to this: 
    group.existing.index <- unlist(lapply(group.i, function(x) which(group.max %in% x)))
    for (net.j in names(net)) {
      values <- net[[net.j]]
      if (net.j %in% c("prob","pval")) {
        values.new <- array(data = 0, dim = c(group.num.max, group.num.max, dim(values)[3]),
                            dimnames = list(group.max, group.max, dimnames(values)[[3]]))
        values.new[group.existing.index, group.existing.index, ] <- values
        net[[net.j]] <- values.new
      }
      if (net.j %in% c("count","sum","weight")) {
        values.new <- array(data = 0, dim = c(group.num.max, group.num.max),
                            dimnames = list(group.max, group.max))
        values.new[group.existing.index, group.existing.index] <- values
        net[[net.j]] <- values.new
      }
      if (net.j %in% c("pairwiseRank")) {
        for (k in 1:length(values)) {
          values.new1 <- vector("list", group.num.max)
          values.new1[group.existing.index] <- values[[k]]
          temp <- values[[k]][[1]]
          temp$prob  <- 0; temp$pval <- 1
          for (kk in setdiff(1:group.num.max, group.existing.index)) {
            values.new1[[kk]] <- temp
          }
          names(values.new1) <- group.max
          values[[k]] <- values.new1
        }
        values.new <- vector("list", group.num.max)
        values.new[group.existing.index] <- values
        temp <- lapply(values.new1, function(x) {
          x$prob  <- 0; x$pval <- 1
          return(x)
        })
        for (kk in setdiff(1:group.num.max, group.existing.index)) {
          values.new[[kk]] <- temp
        }
        names(values.new) <- group.max
      }
      net[[net.j]] <- values.new
    }
    object@net <- net
    
    
    # cat("Update slot object@netP...", '\n')
    netP <- object@netP
    for (netP.j in names(netP)) {
      values <- netP[[netP.j]]
      if (netP.j %in% c("pathways")) {
        values.new <- values
        netP[[netP.j]] <- values.new
      }
      if (netP.j %in% c("prob")) {
        values.new <- array(data = 0, dim = c(group.num.max, group.num.max, dim(values)[3]),
                            dimnames = list(group.max, group.max, dimnames(values)[[3]]))
        values.new[group.existing.index, group.existing.index, ] <- values
        netP[[netP.j]] <- values.new
      }
      if (netP.j %in% c("centr")) {
        for (k in 1:length(values)) {
          values.new <- lapply(values, function(x) {
            values.new2 <- lapply(x, function(x) {
              values.new1 = as.vector(matrix(0, nrow = 1, ncol = group.num.max))
              values.new1[group.existing.index] <- x
              names(values.new1) <- group.max
              return(values.new1)
            })
            names(values.new2) <- names(x)
            return(values.new2)
          })
          names(values.new) <- names(values)
        }
      }
      netP[[netP.j]] <- values.new
    }
    object@netP <- netP
    
    # cat("Update slot object@idents...", '\n')
    idents <- factor(idents, levels = group.max)
    object@idents <- idents
  }
  
  return(object)
}





cc_all <- liftCellChat_ADAPTED(cc_all, c("CLP", "erythroid_cell", "immatureB", "matureB",  
                                 "myeloid", "NK_T", "pre_proB", "preB_I", "preB_II", "proB"))



netVisual_diffInteraction(cc_all, weight.scale = T)
netVisual_diffInteraction(cc_all, weight.scale = T, measure = "weight")
netVisual_heatmap(cc_all)
netVisual_heatmap(cc_all, measure = "weight")

rankNet(cc_all, mode = "comparison", stacked = T, do.stat = TRUE)

# combined visualization  

mat2 <- cc_mef@net$weight
mat1 <- cc_tcf@net$weight
cc_tcf <- liftCellChat_ADAPTED(cc_tcf, unique(c(rownames(mat2), rownames(mat1))))
cc_mef <- liftCellChat_ADAPTED(cc_mef, unique(c(rownames(mat2), rownames(mat1))))


groupSize_tcf <- as.numeric(table(cc_tcf@idents))
groupSize_mef <- as.numeric(table(cc_mef@idents))
net_circle_dir <- jhtools::checkdir(glue::glue("{project_dir}/cellTypes_net_circle2"))


for (i in 1:nrow(mat1)) {
  mat1_u <- matrix(0, nrow = nrow(mat1), ncol = ncol(mat1), dimnames = dimnames(mat1))
  mat1_u[i, ] <- mat1[i, ]
  mat2_u <- matrix(0, nrow = nrow(mat2), ncol = ncol(mat2), dimnames = dimnames(mat2))
  mat2_u[i, ] <- mat2[i, ]
  png(glue::glue("{net_circle_dir}/{rownames(mat1)[i]}_weight_netVisual_circ_tcf.png"), width = 5, height = 5, units = "in", res = 600)
  netVisual_circle(mat1_u, vertex.weight = groupSize_tcf, weight.scale = T, arrow.width = 2.5,
                   arrow.size = 0.3, edge.weight.max = max(mat1), title.name = rownames(mat1)[i])
  dev.off()
  png(glue::glue("{net_circle_dir}/{rownames(mat1)[i]}_weight_netVisual_circ_mef.png"), width = 5, height = 5, units = "in", res = 600)
  netVisual_circle(mat2_u, vertex.weight = groupSize_mef, weight.scale = T, arrow.width = 2.5,
                   arrow.size = 0.3, edge.weight.max = max(mat2), title.name = rownames(mat2)[i])
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










cc_tcf <- createCellChat(object = dat_tcf, meta = meta_tcf, group.by = "all")
cc_tcf@DB <- CellChatDB.human ## get DB data
cc_tcf <- subsetData(cc_tcf) ## filter for signaling genes
cc_mef <- createCellChat(object = dat_mef, meta = meta_mef, group.by = "all")
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

readr::write_rds(cc_tcf, "/cluster/home/yjliu_jh/projects/mef2d/output/cellchat_TCF3.rds")
readr::write_rds(cc_mef, "/cluster/home/yjliu_jh/projects/mef2d/output/cellchat_MEF2D.rds")




