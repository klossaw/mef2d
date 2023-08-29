pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr", "igraph",
          "Seurat", "anndata", "stats", "pheatmap", "pracma", "SeuratDisk"
)
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}

sce_tumor2 <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/data/sce_tumor_final.rds")
Idents(sce_tumor2) <- "orig.ident"
sce_tumor <- subset(sce_tumor2, `orig.ident` %notin% "B069_Dx_CD19")

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes


# extensive tests
# current conclusions: scaledata will eliminate the similarities between M samples
# less dims lead to better data merge
# runHarmony



sce_test <- subset(sce_tumor2, `orig.ident` %notin% "B069_Dx_CD19")
sce_test <- sce_test %>% RunUMAP(dims = 1:15, reduction.name = "umap15")
sce_test <- sce_test %>% RunHarmony("group") %>%
  RunUMAP(dims = 1:10, reduction.name = "umap10") %>%
  RunUMAP(dims = 1:20, reduction.name = "umap20") %>% #ScaleData() %>%
  RunUMAP(reduction = "harmony", dims = 1:10, reduction.name = "umap10_h") %>%
  RunUMAP(reduction = "harmony", dims = 1:20, reduction.name = "umap20_h") 

sce_test <- sce_test %>% RunUMAP(features = gene01, reduction.name = "umap01") %>%
  RunUMAP(features = gene03, reduction.name = "umap03") %>%
  RunUMAP(features = gene003, reduction.name = "umap003")

sce_test <- sce_test %>% RunUMAP(reduction = "harmony", features = gene01, reduction.name = "umap01_h") %>%
  RunUMAP(reduction = "harmony", features = gene03, reduction.name = "umap03_h") %>%
  RunUMAP(reduction = "harmony", features = gene003, reduction.name = "umap003_h")
  
sce_test <- sce_test %>% RunUMAP(reduction = "harmony", features = gene02, reduction.name = "umap02_h")
sce_test <- sce_test %>% RunUMAP(reduction = "harmony", features = gene03[1:30], reduction.name = "umap04_h")


DimPlot(sce_test, reduction = "umap10", label = T, group.by = "orig.ident")
DimPlot(sce_test, reduction = "umap20_h", label = T, group.by = "orig.ident")
DimPlot(sce_test, reduction = "umap03", label = T, group.by = "cell_types_final")
DimPlot(sce_test, reduction = "umap01_h", group.by = "orig.ident", split.by = "orig.ident")
DimPlot(sce_test, reduction = "umap02_h", group.by = "orig.ident", split.by = "orig.ident")
DimPlot(sce_test, reduction = "umap02_h", group.by = "orig.ident")

DimPlot(sce_test, reduction = "umap03_h", label = T, group.by = "cell_types_final")
DimPlot(sce_test, reduction = "umap02_h", label = T, group.by = "cell_types_final")
DimPlot(sce_test, reduction = "umap04_h", label = T, group.by = "cell_types_final")




# scaled seurat will 
sce_tumor <- sce_tumor %>% 
  CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE) 
sce_tumor$cc_diff <- sce_tumor$S.Score - sce_tumor$G2M.Score
sce_tumor <- SCTransform(sce_tumor, vars.to.regress = "cc_diff", verbose = FALSE)
# sce_tumor <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/data/sce_tumor_4.rds")
sce_tumor <- sce_tumor %>% RunUMAP(dims = 1:15, reduction.name = "umap15") %>%
  RunUMAP(dims = 1:10, reduction.name = "umap10") %>%
  RunUMAP(dims = 1:20, reduction.name = "umap20") %>%
  RunUMAP(dims = 1:30, reduction.name = "umap30")
sce_tumor <- sce_tumor %>% RunHarmony("group") %>%
  RunUMAP(reduction = "harmony", dims = 1:15, reduction.name = "umap15_h") %>%
  RunUMAP(reduction = "harmony", dims = 1:20, reduction.name = "umap20_h") %>%
  RunUMAP(reduction = "harmony", dims = 1:30, reduction.name = "umap30_h")
sce_tumor <- sce_tumor %>% RunUMAP(reduction = "harmony", dims = 1:10, reduction.name = "umap10_h")

sce_tumor <- sce_tumor %>% RunUMAP(reduction = "harmony", dims = 1:10, min.dist = 0.15, reduction.name = "umap10_h_m015")

sce_tumor <- sce_tumor %>% RunUMAP(reduction = "harmony", features = gene01, reduction.name = "umap01_h") %>%
  RunUMAP(reduction = "harmony", features = gene03, reduction.name = "umap03_h") %>%
  RunUMAP(reduction = "harmony", features = gene02, reduction.name = "umap02_h")

sce_tumor <- sce_tumor %>% RunUMAP(reduction = "harmony", features = gene03[1:33], reduction.name = "umap01_h") 


DimPlot(sce_tumor, reduction = "umap10_h_m015", label = T, group.by = "cell_types_final")
DimPlot(sce_tumor, reduction = "umap01_h", label = T, group.by = "orig.ident")




marrow <- RunPCA(marrow, features = VariableFeatures(marrow), nfeatures.print = 10)




  
  RunHarmony("group") %>%
  RunUMAP(reduction = "harmony", dims = 1:15, reduction.name = "umap15_hb") %>%
  RunUMAP(reduction = "harmony", dims = 1:15, reduction.name = "umap15_hb") %>%
sce_tumor <- sce_tumor %>% 
  FindVariableFeatures(verbose = FALSE) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  FindNeighbors() %>% 
  FindClusters() %>% 
  identity()


sce_tumor <- sce_tumor %>% RunUMAP(reduction = "harmony", dims = 1:15)
sce_tumor <- RunUMAP(sce_tumor, dims = 1:15, reduction.name = "umap15")
sce_tumor <- RunUMAP(sce_tumor, dims = 1:15, reduction.name = "umap10")
sce_tumor <- RunUMAP(sce_tumor, reduction = "harmony", dims = 1:15, reduction.name = "umap15_h")
sce_tumor <- RunUMAP(sce_tumor, reduction = "harmony", dims = 1:30, reduction.name = "umap30_h")
sce_tumor <- RunUMAP(sce_tumor, reduction = "harmony", dims = 1:20, reduction.name = "umap20_h")
sce_tumor <- RunUMAP(sce_tumor, reduction = "harmony", dims = 1:10, reduction.name = "umap10_h")

DimPlot(sce_tumor, reduction = "umap20_h", label = T, group.by = "cell_types_final")

DimPlot(sce_tumor, reduction = "umap10", label = T, group.by = "orig.ident")



Idents(sce_tumor) <- "cell_types"
DimPlot(sce_tumor, reduction = "umap", label = T, group.by = "seurat_clusters")
DimPlot(sce_tumor, reduction = "umap", group.by = "cell_types", 
        split.by = "cell_types")



# re-anno the cells according to harmony umap




  


# clusters













# read new 



















