



sce_tumor_fil <- FindVariableFeatures(sce_tumor_fil)
sce_tumor_fil <- ScaleData(sce_tumor_fil)
sce_tumor_fil <- RunPCA(sce_tumor_fil)
mef_list <- SplitObject(sce_tumor_fil2, split.by = "group2")
mef_anchors <- FindIntegrationAnchors(object.list = mef_list, reduction = "rpca",
                                         anchor.features = 2000, dims = 1:30)

sce_tumor_fil3 <- IntegrateData(anchorset = mef_anchors, dims = 1:30)
DefaultAssay(sce_tumor_fil3) <- "integrated"
Idents(sce_tumor_fil2)

sce_tumor_fil3 <- quick_anno(sce_tumor_fil3, ancol4$cluster)




sce_tumor_fil3 <- RunUMAP(sce_tumor_fil3, dims = 1:30)
DimPlot(sce_tumor_fil3, group.by = 'orig.ident', reduction = 'umap')
Idents(sce_tumor_fil) <- "anno"
DimPlot(sce_tumor_fil2, group.by = 'anno', reduction = 'umap')
DoHeatmap(sce_tumor_fil, features = c(gene03, "MME", "CD19"))


sce_tumor_fil3 <- ScaleData(sce_tumor_fil3, do.center = T, do.scale = F)
sce_tumor_fil3 <- RunPCA(sce_tumor_fil3, npcs = 30, ndims.print = 1:5, nfeatures.print = 5)
sce_tumor_fil3 <- RunUMAP(sce_tumor_fil3, dims = 1:30, reduction = "pca", n.neighbors = 15,
                  min.dist = 0.5, spread = 1, metric = "euclidean", seed.use = 1)  

sce_tumor_fil3 <- RunUMAP(sce_tumor_fil3, dims = 1:30, reduction = "pca", n.neighbors = 15,
                          min.dist = 0.5, spread = 1, metric = "euclidean", seed.use = 1)  

sce_tumor_fil3 <- RunTSNE(sce_tumor_fil3)

sce_tumor_fil3 <- FindNeighbors(sce_tumor_fil3, reduction = "pca", dims = 1:20, k.param = 20)
sce_tumor_fil3 <- FindClusters(sce_tumor_fil3, resolution = 0.8, algorithm = 1, random.seed = 100)

DefaultAssay(sce_tumor_fil3) <- "RNA"  

# Visualize the Louvain clustering and the batches on the UMAP. 
# Remember, the clustering is stored in @meta.data in column seurat_clusters 
# and the technology is stored in the column tech. Remember you can also use DimPlot.
p1 <- DimPlot(sce_tumor_fil3, reduction = "umap", group.by = "seurat_clusters")

DimPlot(sce_tumor_fil3, reduction = "umap", group.by = "orig.ident")
DimPlot(sce_tumor_fil3, reduction = "umap", group.by = "cell_types_all")

DimPlot(sce_tumor_fil3, reduction = "tsne", group.by = "cell_types_all")

Idents(sce_tumor_fil3) <- "seurat_clusters"
DotPlot(object = sce_tumor_fil3, features = gene01) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))

readr::write_rds(sce_tumor_fil3, "/cluster/home/yjliu_jh/projects/mef2d/data/sce_tumor_fil_rpca.rds")

