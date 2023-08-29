pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr", "igraph",
          "Seurat", "Asgard", "celldex", "cmapR", "pracma", "SeuratDisk",
          "cowplot", "edgeR"
)
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}

# set paths 
lincs_dir <- "/cluster/home/jhuang/reference/scell/ASGARD/LINCS"
output_dir <- "/cluster/home/yjliu_jh/projects/mef2d/output/drug"
setwd(output_dir)
set.seed(2023)

PrepareReference(cell.info = glue("{lincs_dir}/GSE70138_Broad_LINCS_cell_info.txt"),
                 gene.info = glue("{lincs_dir}/GSE70138_Broad_LINCS_gene_info.txt"),
                 GSE70138.sig.info = glue("{lincs_dir}/GSE70138_Broad_LINCS_sig_info.txt"),
                 GSE92742.sig.info = glue("{lincs_dir}/GSE92742_Broad_LINCS_sig_info.txt"),
                 GSE70138.gctx = glue("{lincs_dir}/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328.gctx"),
                 GSE92742.gctx = glue("{lincs_dir}/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"),
                 Output.Dir = "output_dir"
)


# construct merged sce object
sce_tumor <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/data/merged_4_tumors.rds")
sce_hd <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/data/healthy_donors.rds")
sce_all <- readr::read_rds("/cluster/home/jhuang/projects/mef2d/analysis/zwn/human/tenx/seurat/merged/all_merged.rds")
new_cell_id <- str_sub(colnames(sce_all), 1 + str_locate(colnames(sce_all), ":")[, 1] / 2)
sce_all <- RenameCells(sce_all, new.names = new_cell_id)

# set basic 
case_type = "MEF2D"
case = c("M2", "M3", "M4")
target_cell = setdiff(unique(sce_hd$cell_types), "erythroid_cell")
control = c("H1_CD19", "H1_CD34pCD19n", "H2_CD19", "H2_CD34pCD19n")

mef_id <- colnames(sce_tumor)[sce_tumor$orig.ident %in% case]
tcf_id <- colnames(sce_tumor)[sce_tumor$orig.ident %in% "E1"]
hd_id <- colnames(sce_hd)

sce_all_mef <- sce_all[, c(mef_id, hd_id)]
sce_all_tcf <- sce_all[, c(tcf_id, hd_id)]

sce <- sce_all_mef
sce@meta.data$celltype <- "other"
sce$celltype[sce$orig.ident %in% case] <- sce_tumor$cell_types_broad[sce_tumor$orig.ident %in% case]
sce$orig.ident[sce$orig.ident %in% "Healthy"] <- sce_hd$orig.ident
sce$celltype[sce$orig.ident %in% control] <- sce_hd$cell_types[sce_hd$orig.ident %in% control]

sce@meta.data$type <- ifelse(sce@meta.data$orig.ident %in% c("B069_Dx_CD19", "M2", "M3", "M4"), "MEF2D", "normal")
sce@meta.data$sample <- sce@meta.data$orig.ident



sce <- sce_all_tcf
sce@meta.data$celltype <- "other"
sce$celltype[sce$orig.ident %in% "E1"] <- sce_tumor$cell_types_broad[sce_tumor$orig.ident %in% "E1"]
sce$orig.ident[sce$orig.ident %in% "Healthy"] <- sce_hd$orig.ident
sce$celltype[sce$orig.ident %in% control] <- sce_hd$cell_types[sce_hd$orig.ident %in% control]

sce@meta.data$type <- ifelse(sce@meta.data$orig.ident %in% c("E1"), "MEF2D", "normal")
sce@meta.data$sample <- sce@meta.data$orig.ident




data <- as.data.frame(sce@meta.data)
celltypes <- unique(data$celltype)
final_table <- data.frame()
for(i in c("normal", case_type)){
  sub_data <- subset(data, type == i)
  severity <- unique(sub_data$type)
  celltype_freq <- round(100 * table(sub_data$celltype) / nrow(sub_data), 2)
  celltype_freq <- celltype_freq[celltypes]
  final_table <- rbind(final_table, celltype_freq)
}
colnames(final_table) <- celltypes
rownames(final_table) <- c("normal", case_type)



# calculate differential expressed genes
# limma
DefaultAssay(sce) <- "RNA"
min.cells <- 3
gene_list <- list()
c_names <- NULL
for (i in unique(sce@meta.data$celltype)) {
  Idents(sce) <- "celltype"
  c_cells <- subset(sce, celltype == i)
  Idents(c_cells) <- "type"
  samples <- c_cells@meta.data
  controlsample <- row.names(subset(samples, sample %in% control))
  casesample <- row.names(subset(samples, sample %in% case))
  if(length(controlsample)>min.cells & length(casesample)>min.cells){
    expr <- c_cells@assays$RNA@data
    new_expr <- expr[, c(casesample, controlsample)]
    new_sample <- data.frame(sample_id = c(casesample, controlsample),
                             type = c(rep("case", length(casesample)), rep("control", length(controlsample))))
    row.names(new_sample) <- paste(new_sample$sample_id, row.names(new_sample), sep = "_")
    expr <- new_expr
    bad <- which(rowSums(expr > 0) < 3)
    if(length(bad) > 0){
      expr <- expr[-bad, ]
    }
    mm <- model.matrix(~0 + type, data = new_sample)
    fit <- lmFit(expr, mm)
    contr <- makeContrasts(typecase - typecontrol, levels = colnames(coef(fit)))
    tmp <- contrasts.fit(fit, contrasts = contr)
    tmp <- eBayes(tmp)
    c_data <- topTable(tmp, sort.by = "P",n = nrow(tmp))
    c_data_for_drug <- data.frame(row.names = row.names(c_data), score = c_data$t, 
                                  adj.P.Val = c_data$adj.P.Val, P.Value = c_data$P.Value)
    gene_list[[i]] <- c_data_for_drug
    c_names <- c(c_names, i)
  }
}
names(gene_list) <- c_names
readr::write_rds(gene_list, file = "MEF2D_genelist_limma.rds")


# Seurat
gene_list <- list()
c_names <- NULL
for (i in unique(sce@meta.data$celltype)){
  Idents(sce) <- "celltype"
  c_cells <- subset(sce, celltype == i)
  Idents(c_cells) <- "type"
  if (length(table(c_cells$type)) > 1) {
    c_data <- FindMarkers(c_cells, ident.1 = case_type, ident.2 = "normal")
    c_data_for_drug <- data.frame(row.names = row.names(c_data), score = c_data$avg_log2FC,
                                  adj.P.Val = c_data$p_val_adj, P.Value = c_data$p_val)
    gene_list[[i]] <- c_data_for_drug
    c_names <- c(c_names, i)
  }

}
names(gene_list) <- c_names
readr::write_rds(gene_list, file = "MEF2D_genelist_seurat.rds")


# Seurat DESeq2
gene_list <- list()
c_names <- NULL
for(i in unique(sce@meta.data$celltype)){
  Idents(sce) <- "celltype"
  c_cells <- subset(sce, celltype == i)
  Idents(c_cells) <- "type"
  c_data <- FindMarkers(c_cells, ident.1 = case_type, ident.2 = "normal", test.use = "DESeq2")
  c_data_for_drug <- data.frame(row.names = row.names(c_data), score = c_data$avg_log2FC,
                                adj.P.Val = c_data$p_val_adj, P.Value = c_data$p_val)
  gene_list[[i]] <- c_data_for_drug
  c_names <- c(c_names, i)
}
names(gene_list) <- c_names
readr::write_rds(gene_list, file = "MEF2D_genelist_DESeq2.rds")
# there is a zero problem, see here: https://www.biostars.org/p/440379/
# currently we can ignore this part


# edgeR
gene_list <- list()
c_names <- NULL
for (i in unique(sce@meta.data$celltype)){
  Idents(sce) <- "celltype"
  c_cells <- subset(sce, celltype == i)
  Idents(c_cells) <- "type"
  samples = c_cells@meta.data
  controlsample <- row.names(subset(samples, sample %in% control))
  casesample <- row.names(subset(samples, sample %in% case))
  if(length(controlsample) > min.cells & length(casesample) > min.cells){
    expr <- as.matrix(c_cells@assays$RNA@data)
    new_expr <- as.matrix(expr[, c(casesample, controlsample)])
    new_sample <- data.frame(samples=c(casesample, controlsample), 
                             type = c(rep("case", length(casesample)), rep("control", length(controlsample))))
    row.names(new_sample) <- paste(new_sample$samples, row.names(new_sample), sep = "_")
    expr <- new_expr
    bad <- which(rowSums(expr > 0) < 3)
    expr <- expr[-bad, ]
    group <- new_sample$type
    dge <- DGEList(counts = expr, group = group)
    group_edgeR <- factor(group, levels = c("control", "case"))
    design <- model.matrix(~ group_edgeR)
    dge <- estimateDisp(dge, design = design)
    fit <- glmFit(dge, design)
    res <- glmLRT(fit)
    c_data <- res$table
    c_data_for_drug <- data.frame(row.names = row.names(c_data), score = c_data$logFC,
                                  adj.P.Val = p.adjust(c_data$PValue,method = "BH"), P.Value = c_data$PValue)
    gene_list[[i]] <- c_data_for_drug
    c_names <- c(c_names,i)
  }
}
names(gene_list) <- c_names
#Save data
readr::write_rds(gene_list, file="MEF2D_genelist_edgeR.rds")



#Drug repurposing
for (i in c("limma", "seurat")){
  gene_list <- readr::read_rds(paste0("MEF2D_genelist_", i, ".rds"))
  my_gene_info <- read.table(file = "haematopoietic-and-lymphoid-tissue_gene_info.txt", sep="\t",header = T,quote = "")
  my_drug_info <- read.table(file = "haematopoietic-and-lymphoid-tissue_drug_info.txt", sep="\t",header = T,quote = "")
  cmap.ref.profiles <- GetDrugRef(drug.response.path = 'haematopoietic-and-lymphoid-tissue_rankMatrix.txt',
                                  probe.to.genes = my_gene_info, drug.info = my_drug_info)
  Drug.ident.res <- GetDrug(gene.data = gene_list, drug.ref.profiles = cmap.ref.profiles, 
                            repurposing.unit = "drug", connectivity = "negative", drug.type = "FDA")
  readr::write_rds(Drug.ident.res, paste0("MEF2D_drugs_FDA_", i, ".rds"))
}


# drug score
GSE92742.gctx.path <- glue::glue("{lincs_dir}/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx")
GSE70138.gctx.path <- glue::glue("{lincs_dir}/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328.gctx")
Tissue <- "haematopoietic and lymphoid tissue"
for(i in c("limma", "seurat")){
  Gene.list <- read_rds(paste0("MEF2D_genelist_", i, ".rds"))
  Drug.ident.res <- read_rds(paste0("MEF2D_drugs_FDA_", i, ".rds"))  
  Drug.score <- DrugScore(SC.integrated = sce,
                        Gene.data = gene_list,
                        Cell.type = NULL,
                        Drug.data = Drug.ident.res,
                        FDA.drug.only = T,
                        Case = case,
                        Tissue = Tissue,
                        GSE92742.gctx = GSE92742.gctx.path,
                        GSE70138.gctx = GSE70138.gctx.path)
  readr::write_rds(Drug.score, file = paste0("MEF2D_drugscore_FDA_", i, ".rds"))
}







#===================================

case <- "E1"

DefaultAssay(sce) <- "RNA"
min.cells <- 3
gene_list <- list()
c_names <- NULL
for (i in unique(sce@meta.data$celltype)) {
  Idents(sce) <- "celltype"
  c_cells <- subset(sce, celltype == i)
  Idents(c_cells) <- "type"
  samples <- c_cells@meta.data
  controlsample <- row.names(subset(samples, sample %in% control))
  casesample <- row.names(subset(samples, sample %in% case))
  if(length(controlsample)>min.cells & length(casesample)>min.cells){
    expr <- c_cells@assays$RNA@data
    new_expr <- expr[, c(casesample, controlsample)]
    new_sample <- data.frame(sample_id = c(casesample, controlsample),
                             type = c(rep("case", length(casesample)), rep("control", length(controlsample))))
    row.names(new_sample) <- paste(new_sample$sample_id, row.names(new_sample), sep = "_")
    expr <- new_expr
    bad <- which(rowSums(expr > 0) < 3)
    if(length(bad) > 0){
      expr <- expr[-bad, ]
    }
    mm <- model.matrix(~0 + type, data = new_sample)
    fit <- lmFit(expr, mm)
    contr <- makeContrasts(typecase - typecontrol, levels = colnames(coef(fit)))
    tmp <- contrasts.fit(fit, contrasts = contr)
    tmp <- eBayes(tmp)
    c_data <- topTable(tmp, sort.by = "P",n = nrow(tmp))
    c_data_for_drug <- data.frame(row.names = row.names(c_data), score = c_data$t, 
                                  adj.P.Val = c_data$adj.P.Val, P.Value = c_data$P.Value)
    gene_list[[i]] <- c_data_for_drug
    c_names <- c(c_names, i)
  }
}
names(gene_list) <- c_names
readr::write_rds(gene_list, file = "TCF3_genelist_limma.rds")


# Seurat


#Drug repurposing
for (i in c("limma")){
  gene_list <- readr::read_rds(paste0("TCF3_genelist_", i, ".rds"))
  my_gene_info <- read.table(file = "haematopoietic-and-lymphoid-tissue_gene_info.txt", sep="\t",header = T,quote = "")
  my_drug_info <- read.table(file = "haematopoietic-and-lymphoid-tissue_drug_info.txt", sep="\t",header = T,quote = "")
  cmap.ref.profiles <- GetDrugRef(drug.response.path = 'haematopoietic-and-lymphoid-tissue_rankMatrix.txt',
                                  probe.to.genes = my_gene_info, drug.info = my_drug_info)
  Drug.ident.res <- GetDrug(gene.data = gene_list, drug.ref.profiles = cmap.ref.profiles, 
                            repurposing.unit = "drug", connectivity = "negative", drug.type = "FDA")
  readr::write_rds(Drug.ident.res, paste0("TCF3_drugs_FDA_", i, ".rds"))
}


# drug score
GSE92742.gctx.path <- glue::glue("{lincs_dir}/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx")
GSE70138.gctx.path <- glue::glue("{lincs_dir}/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328.gctx")
Tissue <- "haematopoietic and lymphoid tissue"
for(i in c("limma")){
  Gene.list <- read_rds(paste0("TCF3_genelist_", i, ".rds"))
  Drug.ident.res <- read_rds(paste0("TCF3_drugs_FDA_", i, ".rds"))  
  Drug.score <- DrugScore(SC.integrated = sce,
                          Gene.data = gene_list,
                          Cell.type = NULL,
                          Drug.data = Drug.ident.res,
                          FDA.drug.only = T,
                          Case = case,
                          Tissue = Tissue,
                          GSE92742.gctx = GSE92742.gctx.path,
                          GSE70138.gctx = GSE70138.gctx.path)
  readr::write_rds(Drug.score, file = paste0("TCF3_drugscore_FDA_", i, ".rds"))
}



Drug.score <- readr::read_rds(paste0("TCF3_drugscore_FDA_", "limma", ".rds"))

#Drug score plot
Score.list<-data.frame(Patient = "E1", Drug = row.names(Drug.score),
                       DrugScore = Drug.score$Drug.therapeutic.score,
                       Pvalue = Drug.score$P.value, FDR=Drug.score$FDR)


Score.list <- Score.list[Score.list$Pvalue < 0.03, ]
Score.list$Drug<-capitalize(Score.list$Drug)
Drug.order<-Score.list$Drug[order(Score.list$DrugScore, decreasing = F)]
Score.list$Drug<-factor(Score.list$Drug, levels = Drug.order)





pdf(file = "Drug_individual_TCF3.pdf",width = 3,height = 3)
ggplot(Score.list,aes(x=Patient, y=Drug, size=DrugScore, color=Patient)) +
  geom_point(alpha=1) +
  scale_size(name="DrugScore")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+ guides(color=FALSE)
dev.off()








# 

Drug.score <- readr::read_rds(paste0("MEF2D_drugscore_FDA_", "limma", ".rds"))
Combined.Drug.score <- readr::read_rds(paste0("MEF2D_drugscore_FDA_", "limma", ".rds"))
Gene.list<-readRDS("MEF2D_genelist_limma.rds")
Drug.ident.res<-readRDS("MEF2D_drugs_FDA_limma.rds")

#Drug score plot
Score.list<-data.frame(Patient = "Overall", Drug = row.names(Drug.score),
                       DrugScore = Drug.score$Drug.therapeutic.score,
                       Pvalue = Drug.score$P.value, FDR=Drug.score$FDR)

j <- 0
for(i in case){
  j=j+1
  Drug.score<-DrugScore(SC.integrated=sce,
                        Gene.data=Gene.list,
                        Cell.type=NULL,
                        Drug.data=Drug.ident.res,
                        FDA.drug.only=T,
                        Case=i,
                        Tissue=Tissue,
                        GSE92742.gctx=GSE92742.gctx.path,
                        GSE70138.gctx=GSE70138.gctx.path)
  Temp=data.frame(Patient=paste0("Patient", j+1),Drug=row.names(Drug.score),DrugScore=Drug.score$Drug.therapeutic.score,Pvalue=Drug.score$P.value,FDR=Drug.score$FDR)
  Score.list=rbind(Score.list,Temp)
}

library(Hmisc)
Score.list$Drug<-capitalize(Score.list$Drug)
sig.list<-subset(Score.list,DrugScore>quantile(Score.list$DrugScore, 0.9,na.rm=T))
Drug.order<-row.names(Combined.Drug.score)[order(Combined.Drug.score$Drug.therapeutic.score,decreasing = F)]
Drug.order<-capitalize(Drug.order)
sig.list$Drug<-factor(sig.list$Drug,levels = Drug.order)
library('ggplot2')
pdf(file = "Drug_individual_MEF2D.pdf",width = 3,height = 3)
ggplot(sig.list,aes(x=Patient, y=Drug, size=DrugScore, color=Patient)) +
  geom_point(alpha=1) +
  scale_size(name="DrugScore")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+ guides(color=FALSE)
dev.off()











#  myeloid (MPO, CSF1R), T-cell (CD7, CD244) and stem cell (SPINK2, PROM1) 


merged <- read_rds("/cluster/home/jhuang/projects/mef2d/analysis/zwn/human/tenx/seurat/merged/mef2d.rds")

temp_id <- rownames(merged@meta.data)[]

merged@meta.data$celltype[]

