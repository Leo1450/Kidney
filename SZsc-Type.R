setwd("E:/sc-Type/Setareh/")

library(dplyr)
library(Seurat)
library(SeuratData)
library(openxlsx)
library(HGNChelper)
library(pheatmap)

obj<- readRDS("F:/PostDoc/GSE183904_obj_harmony_preannot.rds")

DimPlot(obj, reduction = 'tsne', raster=FALSE, label = T, split.by = "status")

source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"; tissue = "Kidney"

gs_list = gene_sets_prepare(db_, tissue)

#es.max = sctype_score(scRNAseqData = sobj[["RNA"]]@scale.data, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
es.max = sctype_score(scRNAseqData = GetAssayData(obj, layer = "scale.data"), scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

cL_resutls = do.call("rbind", lapply(unique(obj@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(obj@meta.data[obj@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(obj@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores) 

sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

obj@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  obj@meta.data$customclassif[obj@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}
DimPlot(obj, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif', raster = F, split.by = "status")  
DimPlot(obj, reduction = "tsne", label = TRUE, repel = TRUE, group.by = 'customclassif', raster = F)  

tab <- table(Assigned=obj$customclassif, Clusters=obj$seurat_clusters)
write.csv(tab, "Table of Annotation (GSE156632).csv")

pheatmap(log10(tab+10), color = colorRampPalette(c('white','blue'))(10))

#lapply(c("ggraph","igraph","tidyverse", "data.tree"), library, character.only = T)
