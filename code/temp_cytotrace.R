pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr", "igraph",
          "Seurat", "CytoTRACE", "reticulate", "monocle", "numDeriv"#, "monocle3"
)
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}


ddr <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/output/ddrtree_object_1.rds")

ddr_test <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/output/ddrtree_object_test.rds")

new_anno_col <- na.omit(unique(ddr_test@phenoData@data$anno))
rownames(ddr_test@phenoData@data) <- ddr_test@phenoData@data$cell_id
ddr_test@phenoData@data$cell_types2 <- ddr_test@phenoData@data$anno
ddr_test@phenoData@data[, c("cell_id", "anno")] <- NULL



colors_x <- list()
for(i in 1:length(new_anno_col)){
  colors_x[[new_anno_col[i]]] <- rep("gray", length(new_anno_col))
  colors_x[[new_anno_col[i]]][i] <- "green"
  names(colors_x[[new_anno_col[i]]]) <- new_anno_col
}



for (i in 1:length(colors_x)){
  ddr_test@phenoData@data$anno2 <- ddr_test@phenoData@data$anno
  ddr_test@phenoData@data$anno2[ddr_test@phenoData@data$anno2 %notin% new_anno_col[[i]]] <- NA
  pm1 <- monocle::plot_cell_trajectory(ddr_test, color_by = "Pseudotime", size = 1, show_backbone = TRUE)
  pm2 <- monocle::plot_cell_trajectory(ddr_test, color_by = "anno2", alpha = 0.3, size = 1, show_backbone = TRUE) +
    scale_color_discrete(colors_x[[i]]) + theme(legend.title = element_blank())
  pm_all <- pm1 + pm2 + plot_layout(nrow = 1)
  pdf(glue::glue("/cluster/home/yjliu_jh/projects/temp/ddrtree_{new_anno_col[[i]]}.pdf"), width = 10, height = 5)
  print(pm_all)
  dev.off()
}


i <- 1
ddrx <- ddr_test[, rownames(ddr_test@phenoData@data)[ddr_test@phenoData@data$cell_types2 %notin% new_anno_col[[i]]]]
pm2 <- monocle::plot_cell_trajectory(ddrx, color_by = "cell_types2", alpha = 0.3, size = 1, show_backbone = TRUE) +
  theme(legend.title = element_blank())



i <- 1
ddrx <- ddr_test[, rownames(ddr_test@phenoData@data)[ddr_test@phenoData@data$cell_types2 %notin% new_anno_col[[i]]]]





pm1 <- monocle::plot_cell_trajectory(ddr, color_by = "Pseudotime", size = 1, show_backbone = TRUE)
pm2 <- monocle::plot_cell_trajectory(ddr, color_by = "cell_types2", size = 1, show_backbone = TRUE) +
  scale_color_discrete(new_color) + theme(legend.title = element_blank())
pm_all <- pm1 + pm2 + plot_layout(nrow = 1)

pdf("/cluster/home/yjliu_jh/projects/temp/whateeeee.pdf", width = 10, height = 5)
pm_all
dev.off()





# ====== repeat the test ======


ddr_test <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/output/ddrtree_object_test.rds")
temp <- ddr_test@phenoData@data
tempanno <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/output/mta.rds")
colnames(tempanno)[2] <- "annox"
temp <- left_join(temp, tempanno)
temp$labels_leu1[is.na(temp$labels_leu1)] <- "NA"
temp$labels_leu1[temp$labels_leu1 %in% c("GMP")] <- "GMP-like"
temp$labels_leu1[temp$labels_leu1 %in% c("earlyEry", "lateEry")] <- "Ery"
temp$labels_leu1[temp$labels_leu1 %in% c("cDC-like", "pDC", "ProMono")] <- "DC-like"
temp$labels_leu1[temp$labels_leu1 %in% c("NK")] <- "T"



ddr_test@phenoData@data$cell_types2 <- temp$labels_leu1
ddr_test@phenoData@data$cell_types2[is.na(ddr_test@phenoData@data$cell_types2)] <- "NA"
  
samples <- unique(temp$labels_leu1)
sample_color <- c("#ee5d5a", "#9f2f6a", "#34956c", "#6b6498","#b29bc9", "#71a3a2", 
                           "#c88978",  "#a15217", "#ce662f", "#cc9b32", "#53096a", "#6569c2", "#66676c")[1:length(samples)]
names(sample_color) <- samples
                           
pm2 <- monocle::plot_cell_trajectory(ddr_test, color_by = "cell_types2", alpha = 0.3, size = 1, show_backbone = TRUE) +
  theme(legend.title = element_blank()) + scale_color_discrete(sample_color) + facet_wrap("~cell_types2", nrow = 3)



monocle::plot_cell_trajectory(ddr_test, color_by = "cell_types", alpha = 0.3, size = 1, show_backbone = TRUE) +
  theme(legend.title = element_blank()) + scale_color_discrete(sample_color) + facet_wrap("~cell_types", nrow = 3)



