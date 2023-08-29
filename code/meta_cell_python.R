pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr", "igraph",
          "Seurat", "anndata", "stats", "pheatmap", "pracma", "SeuratDisk"
)
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}

# 1 assist the final annotation of cells
# 2 may also aid for later analysis 
# 3 test if this can be incorporated in the reference annotation pipeline (when with spare time)

sce <- readr::read_rds("~/projects/mef2d/analysis/zwn/human/tenx/seurat/merged/mef2d.rds")

# step 1: Slim down a Seurat object. So you get raw counts, lognorm counts
seu_sct = DietSeurat(
  sce,
  counts = TRUE, # so, raw counts save to adata.raw.X 
  data = TRUE, # so, log1p counts save to adata.X
  scale.data = FALSE, # set to false, or else will save to adata.X
  features = rownames(sce), # export all genes, not just top highly variable genes
  assays = "SCT",
  dimreducs = c("pca","umap"),
  graphs = c("RNA_nn", "RNA_snn"), # to RNA_nn -> distances, RNA_snn -> connectivities
  misc = TRUE
)


# step 2: factor to character, or else your factor will be number in adata 
i <- sapply(seu_sct@meta.data, is.factor)
seu_sct@meta.data[i] <- lapply(seu_sct@meta.data[i], as.character)

# step 3: convert 
SaveH5Seurat(seu_sct, filename = "/cluster/home/yjliu_jh/projects/mef2d/data/mef2d_all_sct.h5seurat", overwrite = TRUE)
Convert("/cluster/home/yjliu_jh/projects/mef2d/data/mef2d_all_sct.h5seurat", 
        "/cluster/home/yjliu_jh/projects/mef2d/data/mef2d_all_sct.h5ad", assay = "SCT", overwrite = TRUE)


# === RNA assay 
DefaultAssay(sce) <- "RNA"
seu_rna = DietSeurat(
  sce, counts = TRUE, data = TRUE, scale.data = FALSE, features = rownames(sce),
  assays = "RNA",
  dimreducs = c("pca","umap"), graphs = c("RNA_nn", "RNA_snn"), misc = TRUE
)


# step 2: factor to character, or else your factor will be number in adata 
i <- sapply(seu_rna@meta.data, is.factor)
seu_rna@meta.data[i] <- lapply(seu_rna@meta.data[i], as.character)

# step 3: convert 
SaveH5Seurat(seu_rna, filename = "/cluster/home/yjliu_jh/projects/mef2d/data/mef2d_all_rna.h5seurat", overwrite = TRUE)
Convert("/cluster/home/yjliu_jh/projects/mef2d/data/mef2d_all_rna.h5seurat", 
        "/cluster/home/yjliu_jh/projects/mef2d/data/mef2d_all_rna.h5ad", assay = "RNA", overwrite = TRUE)



# === only tumor, RNA assay

DefaultAssay(sce_tumor) <- "RNA"
stm_rna = DietSeurat(
  sce_tumor, counts = TRUE, data = TRUE, scale.data = FALSE, features = rownames(sce_tumor),
  assays = "RNA",
  dimreducs = c("pca","umap"), graphs = c("RNA_nn", "RNA_snn"), misc = TRUE
)


# step 2: factor to character, or else your factor will be number in adata 
i <- sapply(stm_rna@meta.data, is.factor)
stm_rna@meta.data[i] <- lapply(stm_rna@meta.data[i], as.character)

# step 3: convert 
SaveH5Seurat(stm_rna, filename = "/cluster/home/yjliu_jh/projects/mef2d/data/mef2d_tumor_rna.h5seurat", overwrite = TRUE)
Convert("/cluster/home/yjliu_jh/projects/mef2d/data/mef2d_tumor_rna.h5seurat", 
        "/cluster/home/yjliu_jh/projects/mef2d/data/mef2d_tumor_rna.h5ad", assay = "RNA", overwrite = TRUE)






