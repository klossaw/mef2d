pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "cowplot",
          "Seurat", "vroom", "harmony", "SingleR", "scuttle", "Matrix")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}

gse148_files <- list.files("/cluster/home/yjliu_jh/projects/mef2d/data/public/GSE148218_RAW", full.names = T)

bar_files <- gse148_files[grepl("barcode", gse148_files)]
fea_files <- gse148_files[grepl("features", gse148_files)]
mtx_files <- gse148_files[grepl("matrix", gse148_files)]


matrices <- lapply(mtx_files, readMM)
mtx_names <- sub("_matrix.mtx.gz", "", basename(mtx_files))
names(matrices) <- sub(".*?\\_", "", mtx_names)

features <- lapply(fea_files, read.delim, header = F)
barcodes <- lapply(bar_files, read.delim, header = F)

for (i in 1:6) {  ## only E/R samples without treatment are selected
  colnames(matrices[[i]]) <- paste0(names(matrices)[i], ":", barcodes[[i]][[1]])
  rownames(matrices[[i]]) <- features[[i]][[2]]
  matrices[[i]] <- as(matrices[[i]], "matrix")
}

er_matrix <- do.call(cbind, matrices[1:6])
er <- CreateSeuratObject(er_matrix)
er@meta.data$orig.ident <- str_extract(rownames(er@meta.data), "[^:]+")


sce_hd <- readr::read_rds("/cluster/home/jhuang/projects/mef2d/analysis/zwn/human/tenx/seurat/merge/sce_hd.rds")
sce_tumor <- readr::read_rds("/cluster/home/jhuang/projects/mef2d/analysis/zwn/human/tenx/seurat/merge/merge.rds")


obj <- NormalizeData(er)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)

obj <- FindNeighbors(obj, dims = 1:30, reduction = "pca")
obj <- FindClusters(obj, resolution = 2, cluster.name = "unintegrated_clusters")
obj <- RunUMAP(obj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")



pdf("/cluster/home/yjliu_jh/projects/temp/test_dimplot_erdata_01.pdf", height = 10, width = 18)
DimPlot(obj, reduction = "umap.unintegrated", group.by = c("orig.ident", "seurat_clusters"))
dev.off()




# and harmony and dotplot and 


