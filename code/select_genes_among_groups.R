

# credit to the author of MUDAN package
# no, so slow, don't use that function


# first construct trees for select pairs for comparison
# abm_assay <- t(abm_assay)
# abm_celltypes <- abm_coldata$cell_label
# names(abm_celltypes) <- abm_coldata$cname
# 
# Variations <- matrixStats::rowVars(as.matrix(abm_assay))
# gc()
# 
# # don't cbind! avoid to generate large data 
# Variations <- matrixStats::rowVars(as.matrix(abm_assay))
# var_order <- order(Variations)
# 
# abm_assay_sub <- abm_assay[tail(var_order, 5000), ]
# #abm_assay_sub <- abm_assay[tail(var_order, 1000), ]
# 
# 
# abm_group <- do.call(cbind, lapply(unique(abm_celltypes), function(ct) { rowMeans(abm_assay_sub[, abm_celltypes == ct]) }))
# colnames(abm_group) <- unique(abm_celltypes)
# 
# hc_abm <- hclust(dist(t(abm_group)), method = "complete")
# plot(hc_abm)
# 
# dend_abm <- as.dendrogram(hc_abm)
# 
# pv.sig.all <- c()
# 
# pv.recur <- function(dend) {
#   g1 <- labels(dend[[1]])
#   g2 <- labels(dend[[2]])
#   #print(g1)
#   #print(g2)
#   ingroup <- names(abm_celltypes)[abm_celltypes %in% g1]
#   outgroup <- names(abm_celltypes)[abm_celltypes %in% g2]
#   pv <- sapply(1:5000, function(i) {
#     t.test(abm_assay_sub[i, ingroup], abm_assay_sub[i, outgroup])$p.value
#   })
#   names(pv) <- rownames(abm_assay_sub)
#   pv.sig <- names(pv)[pv < 0.05/5000/length(hc_abm$height)] ## bonferonni
#   #print(pv.sig)
#   pv.sig.all <<- c(pv.sig.all, pv.sig) ## save
#   
#   ## recursion to go down tree if not leaf
#   if(!is.leaf(dend[[1]])) {
#     pv.recur(dend[[1]])
#   }
#   if(!is.leaf(dend[[2]])) {
#     pv.recur(dend[[2]])
#   }
# }
# 
# pv.recur(dend_abm)
# readr::write_rds(pv.sig.all, "/cluster/home/yjliu_jh/projects/temp/sig_genes_among_groups_230202.rds")




# turns out that method is too slow... try to use Seurat's FindAllMarkers()



abm_assay <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/data/public/rds/adultbm_cells_matrix_sub.rds")
abm_coldata <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/data/public/rds/adultbm_cells_metadata.rds")
abm_assay <- t(abm_assay)
abm_assay <- abm_assay[, !duplicated(colnames(abm_assay))]
abm_coldata <- abm_coldata %>% dplyr::select("cname", everything())
rownames(abm_coldata) <- abm_coldata$cname
abm <- CreateSeuratObject(counts = abm_assay, project = "abm_anno", meta.data = abm_coldata,
                          min.cells = 3, min.features = 200)
#abm <- SetAllIdent(object = abm, id = "cell_label")
Idents(abm) <- abm@meta.data$cell_identity
abm <- NormalizeData(abm)
all_markers_abm <- FindAllMarkers(abm)
readr::write_rds(all_markers_abm, "/cluster/home/yjliu_jh/projects/mef2d/data/public/rds/all_markers_abm_1.rds")

# the results could be used also for generating signature matrix for cibersort


fbm_assay <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/data/public/rds/cordblood_cells_matrix_sub.rds")
fbm_coldata <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/data/public/rds/cordblood_cells_metadata.rds")
fbm_assay <- t(fbm_assay)
fbm_assay <- fbm_assay[, colnames(fbm_assay) %notin% colnames(fbm_assay)[duplicated(colnames(fbm_assay))]]
fbm_coldata <- fbm_coldata %>% dplyr::select("cname", everything())
rownames(fbm_coldata) <- fbm_coldata$cname
fbm <- CreateSeuratObject(counts = fbm_assay, project = "fbm_anno", meta.data = fbm_coldata,
                          min.cells = 3, min.features = 200)
#abm <- SetAllIdent(object = abm, id = "cell_label")
Idents(fbm) <- fbm@meta.data$cell_identity
fbm <- NormalizeData(fbm)
all_markers_fbm <- FindAllMarkers(fbm)
readr::write_rds(all_markers_fbm, "/cluster/home/yjliu_jh/projects/mef2d/data/public/rds/all_markers_fbm_1.rds")




