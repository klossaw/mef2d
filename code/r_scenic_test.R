pkgs <- c("SCENIC", "SCopeLoomR", "patchwork", "ggalluvial", "svglite", "Seurat")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}


defaultDbNames[["hgnc"]][1] <- "hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather"
defaultDbNames[["hgnc"]][2] <- "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather"

scenicOptions <- initializeScenic(org = "hgnc", nCores = 20,  dbs = defaultDbNames[["hgnc"]], 
                                  dbDir = "/cluster/home/yjliu_jh/projects/mef2d/data/public/cistarget/v1")


sce_tumor <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/output/sce_tumor.rds")
Idents(sce_tumor) <- "cluster_final"
sce_tumor_fil <- subset(x = sce_tumor, subset = cluster_final %in% 1:18)
exp_mat <- sce_tumor_fil@assays$RNA@data
cell_info <- sce_tumor@meta.data




lincs_dir <- "/cluster/home/yjliu_jh/projects/mef2d/data/public/LINCS"
setwd(lincs_dir)

PrepareReference(cell.info="GSE70138_Broad_LINCS_cell_info_2017-04-28.txt",
                 gene.info="GSE70138_Broad_LINCS_gene_info_2017-03-06.txt",
                 GSE70138.sig.info = "GSE70138_Broad_LINCS_sig_info_2017-03-06.txt",
                 GSE92742.sig.info = "GSE92742_Broad_LINCS_sig_info.txt",
                 GSE70138.gctx = "GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx",
                 GSE92742.gctx = "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx",
                 Output.Dir = "DrugReference/")


