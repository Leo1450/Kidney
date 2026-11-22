setwd("E:/scRNA/SZ/Results/harmony/")

library(Seurat)
library(tidyverse)
library(BiocManager)
library(metap)
library(harmony)
library(patchwork)
library(cowplot)
library(chunked)
library(dplyr)

# Read
CTRL1<- read.csv("E:/scRNA/SZ/para-tumor/GSM4735365_RCC1n.csv")
### Remove duplicate
CTRL1<- distinct(CTRL1, Symbol, .keep_all= TRUE)
row.names(CTRL1)<- CTRL1[ ,1]
CTRL1<- CTRL1[,-1]

CTRL2<- read.csv("E:/scRNA/SZ/para-tumor/GSM4735367_RCC2n.csv")
### Remove duplicate
CTRL2<- distinct(CTRL2, Symbol, .keep_all= TRUE)
row.names(CTRL2)<- CTRL2[ ,1]
CTRL2<- CTRL2[,-1]

CTRL3<- read.csv("E:/scRNA/SZ/para-tumor/GSM4735369_RCC3n.csv")
### Remove duplicate
CTRL3<- distinct(CTRL3, Symbol, .keep_all= TRUE)
row.names(CTRL3)<- CTRL3[ ,1]
CTRL3<- CTRL3[,-1]

CTRL4<- read.csv("E:/scRNA/SZ/para-tumor/GSM4735371_RCC4n.csv")
### Remove duplicate
CTRL4<- distinct(CTRL4, Symbol, .keep_all= TRUE)
row.names(CTRL4)<- CTRL4[ ,1]
CTRL4<- CTRL4[,-1]

CTRL5<- read.csv("E:/scRNA/SZ/para-tumor/GSM4735373_RCC5n.csv")
### Remove duplicate
CTRL5<- distinct(CTRL5, Symbol, .keep_all= TRUE)
row.names(CTRL5)<- CTRL5[ ,1]
CTRL5<- CTRL5[,-1]

CTRL6<- read.csv("E:/scRNA/SZ/para-tumor/GSM4735375_RCC7t.csv")
### Remove duplicate
CTRL6<- distinct(CTRL6, Symbol, .keep_all= TRUE)
row.names(CTRL6)<- CTRL6[ ,1]
CTRL6<- CTRL6[,-1]


RCC1<- read.csv("E:/scRNA/SZ/tumor/GSM4735364_RCC1t.csv")
### Remove duplicate
RCC1<- distinct(RCC1, Symbol, .keep_all= TRUE)
row.names(RCC1)<- RCC1[ ,1]
RCC1<- RCC1[,-1]

RCC2<- read.csv("E:/scRNA/SZ/tumor/GSM4735366_RCC2t.csv")
### Remove duplicate
RCC2<- distinct(RCC2, Symbol, .keep_all= TRUE)
row.names(RCC2)<- RCC2[ ,1]
RCC2<- RCC2[,-1]

RCC3<- read.csv("E:/scRNA/SZ/tumor/GSM4735368_RCC3t.csv")
### Remove duplicate
RCC3<- distinct(RCC3, Symbol, .keep_all= TRUE)
row.names(RCC3)<- RCC3[ ,1]
RCC3<- RCC3[,-1]

RCC4<- read.csv("E:/scRNA/SZ/tumor/GSM4735370_RCC4t.csv")
### Remove duplicate
RCC4<- distinct(RCC4, Symbol, .keep_all= TRUE)
row.names(RCC4)<- RCC4[ ,1]
RCC4<- RCC4[,-1]

RCC5<- read.csv("E:/scRNA/SZ/tumor/GSM4735372_RCC5t.csv")
### Remove duplicate
RCC5<- distinct(RCC5, Symbol, .keep_all= TRUE)
row.names(RCC5)<- RCC5[ ,1]
RCC5<- RCC5[,-1]

RCC6<- read.csv("E:/scRNA/SZ/tumor/GSM4735374_RCC6t.csv")
### Remove duplicate
RCC6<- distinct(RCC6, Symbol, .keep_all= TRUE)
row.names(RCC6)<- RCC6[ ,1]
RCC6<- RCC6[,-1]


#CTRL1 <- GetAssayData(CTRL1, assay = "RNA")
# Create Seurat Object
CTRL1 <- CreateSeuratObject(counts = CTRL1, project = "CTRL", min.cells = 3, min.features = 200)
CTRL1 <- PercentageFeatureSet(CTRL1, pattern = "^MT-", col.name = "percent.mt")

CTRL2 <- CreateSeuratObject(counts = CTRL2, project = "CTRL", min.cells = 3, min.features = 200)
CTRL2 <- PercentageFeatureSet(CTRL2, pattern = "^MT-", col.name = "percent.mt")

CTRL3 <- CreateSeuratObject(counts = CTRL3, project = "CTRL", min.cells = 3, min.features = 200)
CTRL3 <- PercentageFeatureSet(CTRL3, pattern = "^MT-", col.name = "percent.mt")

CTRL4 <- CreateSeuratObject(counts = CTRL4, project = "CTRL", min.cells = 3, min.features = 200)
CTRL4 <- PercentageFeatureSet(CTRL4, pattern = "^MT-", col.name = "percent.mt")

CTRL5 <- CreateSeuratObject(counts = CTRL5, project = "CTRL", min.cells = 3, min.features = 200)
CTRL5 <- PercentageFeatureSet(CTRL5, pattern = "^MT-", col.name = "percent.mt")

CTRL6 <- CreateSeuratObject(counts = CTRL6, project = "CTRL", min.cells = 3, min.features = 200)
CTRL6 <- PercentageFeatureSet(CTRL6, pattern = "^MT-", col.name = "percent.mt")


RCC1 <- CreateSeuratObject(counts = RCC1, project = "RCC", min.cells = 3, min.features = 200)
RCC1 <- PercentageFeatureSet(RCC1, pattern = "^MT-", col.name = "percent.mt")

RCC2 <- CreateSeuratObject(counts = RCC2, project = "RCC", min.cells = 3, min.features = 200)
RCC2 <- PercentageFeatureSet(RCC2, pattern = "^MT-", col.name = "percent.mt")

RCC3 <- CreateSeuratObject(counts = RCC3, project = "RCC", min.cells = 3, min.features = 200)
RCC3 <- PercentageFeatureSet(RCC3, pattern = "^MT-", col.name = "percent.mt")

RCC4 <- CreateSeuratObject(counts = RCC4, project = "RCC", min.cells = 3, min.features = 200)
RCC4 <- PercentageFeatureSet(RCC4, pattern = "^MT-", col.name = "percent.mt")

RCC5 <- CreateSeuratObject(counts = RCC5, project = "RCC", min.cells = 3, min.features = 200)
RCC5 <- PercentageFeatureSet(RCC5, pattern = "^MT-", col.name = "percent.mt")

RCC6 <- CreateSeuratObject(counts = RCC6, project = "RCC", min.cells = 3, min.features = 200)
RCC6 <- PercentageFeatureSet(RCC6, pattern = "^MT-", col.name = "percent.mt")

# Filter low quality cells
VlnPlot(CTRL1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
CTRL1<- subset(CTRL1, subset = nFeature_RNA < 3000 & nCount_RNA < 10000 & percent.mt < 15)

VlnPlot(CTRL2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
CTRL2<- subset(CTRL2, subset = nFeature_RNA < 5000 & nCount_RNA < 20000 & percent.mt < 15)

VlnPlot(CTRL3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
CTRL3<- subset(CTRL3, subset = nFeature_RNA < 2000 & nCount_RNA < 10000 & percent.mt < 10)

VlnPlot(CTRL4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
CTRL4<- subset(CTRL4, subset = nFeature_RNA < 2000 & nCount_RNA < 20000 & percent.mt < 10)

VlnPlot(CTRL5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
CTRL5<- subset(CTRL5, subset = nFeature_RNA < 2000 & nCount_RNA < 5000 & percent.mt < 10)

VlnPlot(CTRL6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
CTRL6<- subset(CTRL6, subset = nFeature_RNA < 5000 & nCount_RNA < 20000 & percent.mt < 15)


VlnPlot(RCC1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
RCC1 <- subset(RCC1, subset = nFeature_RNA < 4000 & nCount_RNA < 15000 & percent.mt < 15)

VlnPlot(RCC2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
RCC2 <- subset(RCC2, subset = nFeature_RNA < 5000 & nCount_RNA < 20000 & percent.mt < 15)

VlnPlot(RCC3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
RCC3 <- subset(RCC3, subset = nFeature_RNA < 4000 & nCount_RNA < 15000 & percent.mt < 7.5)

VlnPlot(RCC4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
RCC4 <- subset(RCC4, subset = nFeature_RNA < 7500 & nCount_RNA < 50000 & percent.mt < 10)

VlnPlot(RCC5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
RCC5 <- subset(RCC5, subset = nFeature_RNA < 4000 & nCount_RNA < 15000 & percent.mt < 10)

VlnPlot(RCC6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
RCC6 <- subset(RCC6, subset = nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 10)


CTRL1 <- NormalizeData(CTRL1, verbose = FALSE)
CTRL1 <- FindVariableFeatures(CTRL1, selection.method = "vst", nfeatures = 2000)

CTRL2 <- NormalizeData(CTRL2, verbose = FALSE)
CTRL2 <- FindVariableFeatures(CTRL2, selection.method = "vst", nfeatures = 2000)

CTRL3 <- NormalizeData(CTRL3, verbose = FALSE)
CTRL3 <- FindVariableFeatures(CTRL3, selection.method = "vst", nfeatures = 2000)

CTRL4 <- NormalizeData(CTRL4, verbose = FALSE)
CTRL4 <- FindVariableFeatures(CTRL4, selection.method = "vst", nfeatures = 2000)

CTRL5 <- NormalizeData(CTRL5, verbose = FALSE)
CTRL5 <- FindVariableFeatures(CTRL5, selection.method = "vst", nfeatures = 2000)

CTRL6 <- NormalizeData(CTRL6, verbose = FALSE)
CTRL6 <- FindVariableFeatures(CTRL6, selection.method = "vst", nfeatures = 2000)


RCC1 <- NormalizeData(RCC1, verbose = FALSE)
RCC1 <- FindVariableFeatures(RCC1, selection.method = "vst", nfeatures = 2000)

RCC2 <- NormalizeData(RCC2, verbose = FALSE)
RCC2 <- FindVariableFeatures(RCC2, selection.method = "vst", nfeatures = 2000)

RCC3 <- NormalizeData(RCC3, verbose = FALSE)
RCC3 <- FindVariableFeatures(RCC3, selection.method = "vst", nfeatures = 2000)

RCC4 <- NormalizeData(RCC4, verbose = FALSE)
RCC4 <- FindVariableFeatures(RCC4, selection.method = "vst", nfeatures = 2000)

RCC5 <- NormalizeData(RCC5, verbose = FALSE)
RCC5 <- FindVariableFeatures(RCC5, selection.method = "vst", nfeatures = 2000)

RCC6 <- NormalizeData(RCC6, verbose = FALSE)
RCC6 <- FindVariableFeatures(RCC6, selection.method = "vst", nfeatures = 2000)


CTRL<- merge(CTRL1 ,c(CTRL2, CTRL3, CTRL4, CTRL5, CTRL6), 
             add.cell.ids = c("CTRL1", "CTRL2", "CTRL3", "CTRL4",
                              "CTRL5", "CTRL6"))


RCC<- merge(RCC1 ,c(RCC2, RCC3, RCC4, RCC5, RCC6),
            add.cell.ids = c("RCC1", "RCC2", "RCC3", "RCC4",
                             "RCC5", "RCC6"))

rm(CTRL1,CTRL2, CTRL3, CTRL4, CTRL5, CTRL6, RCC1 ,RCC2, RCC3, RCC4, RCC5, RCC6)

# Add cell type labels
CTRL$status <- "CTRL"
RCC$status <- "RCC"

# merge into one obj
obj <- merge(x= CTRL, y=list(RCC))


###Now let’s do the same as we did before:
sobj <- NormalizeData(obj, verbose = F)
sobj <- FindVariableFeatures(sobj, selection.method = "vst", nfeatures = 2000, verbose = F)
sobj <- ScaleData(sobj, verbose = F)
sobj <- RunPCA(sobj, npcs = 30, verbose = F)
sobj <- RunUMAP(sobj, reduction = "pca", dims = 1:30, verbose = F)

###Let’s plot the “before” plot again:

DimPlot(sobj, raster = FALSE, reduction = "umap") + plot_annotation(title = "CTRL, RCC and Rel, before integration")

sobj <- sobj %>% RunHarmony("orig.ident", plot_convergence = T)

###Check the generated embeddings:
harmony_embeddings <- Embeddings(sobj, 'harmony')
harmony_embeddings[1:5, 1:5]

###Check the PCA plot after:
p1 <- DimPlot(object = sobj, raster = FALSE, reduction = "harmony", pt.size = .1, group.by = "orig.ident") + NoLegend()
p2 <- VlnPlot(object = sobj, raster = FALSE, features = "harmony_1", group.by = "orig.ident", pt.size = .1) + NoLegend()
plot_grid(p1,p2)

###Do UMAP and clustering:
sobj <- sobj %>% 
  RunUMAP(reduction = "harmony", dims = 1:30, verbose = F) %>% 
  FindNeighbors(reduction = "harmony", k.param = 10, dims = 1:30) %>% 
  FindClusters(resolution = 0.4) %>% 
  identity()

###Finally, so same UMAP plots of integrated datasets as above:
sobj <- SetIdent(sobj,value = "orig.ident")
DimPlot(sobj, raster=FALSE, reduction = "umap") + plot_annotation(title = "CTRL, RCC and Rel, after integration")

DimPlot(sobj, raster=FALSE,reduction = "umap", group.by = "orig.ident", pt.size = .1, split.by = 'orig.ident') + NoLegend()

###These look a bit worse than ones from Seurat:
sobj <- SetIdent(sobj,value = "seurat_clusters")
DimPlot(sobj,raster=FALSE,label = T) + NoLegend()


DimPlot(sobj, raster=FALSE, reduction = "umap", split.by = "status", label = T)

saveRDS(sobj, "sobj(Harmony).RDS")
sobj <- readRDS("sobj(Harmony).RDS")

length(levels(sobj))
length(levels(sobj$seurat_clusters))

count_table <- table(sobj@meta.data$seurat_clusters, sobj@meta.data$orig.ident)
count_table
write.csv(count_table, "count_table.csv")

head(sobj)
tail(sobj)

sobj <- RenameIdents(sobj, `0` = "Endothelial", `1` = "Hematopoietic", `2` = "Proximal tubule", `3` = "Loop of Henle",
                     `4` = "Proximal tubule", `5` = "Renal interstitium",
                     `6` = "Immune cells", `7` = "Endothelial", `8` = "Hematopoietic", `9` = "Proximal tubule",
                     `10` = "Hematopoietic", `11` = "Proximal tubule", `12` = "Cap mesenchyme",
                     `13` = "α-intercalated cells",  `14` = "Endothelial", `15` = "Proximal tubule", 
                     `16` = "Hematopoietic", `17` = "Stromal cells", `18` = "Hematopoietic",  `19` = "Endothelial")


length(levels(sobj))
length(levels(sobj$seurat_clusters))

DimPlot(sobj, raster=FALSE, label = TRUE)
DimPlot(sobj, raster=FALSE, reduction = "umap", split.by = "status", label = T)

###Identify differential expressed genes across 2 samples
head(sobj)
tail(sobj)

sobj<- JoinLayers(sobj)
markers <- FindMarkers(sobj, ident.1 = "RCC", ident.2 = "CTRL", group.by = "status", min.pct = 0.1, logfc.threshold = 0.25)
head(markers, n = 10)
write.csv(markers, "DEG-total(RCC).csv")


#Identify differential expressed genes across conditions
sobj$celltype.status <- paste(Idents(sobj), sobj$status, sep = "_")
sobj$celltype <- Idents(sobj)
Idents(sobj) <- "celltype.status"

length(levels(sobj$celltype))

head(sobj)
tail(sobj)

saveRDS(sobj, "RCC_for_DEGs.RDS")
sobj <- readRDS("TumorCNTR_forSingleR.RDS")

RCCCTRL1 <- FindMarkers(sobj, ident.1 = "Endothelial_RCC", ident.2 = "Endothelial_CTRL", verbose = FALSE)
head(RCCCTRL1, n = 10)
write.csv(RCCCTRL1, "DEG-Endothelial.csv")

RCCCTRL2 <- FindMarkers(sobj, ident.1 = "Hematopoietic_RCC", ident.2 = "Hematopoietic_CTRL", verbose = FALSE)
head(RCCCTRL2, n = 10)
write.csv(RCCCTRL2, "DEG-Hematopoietic.csv")

RCCCTRL3 <- FindMarkers(sobj, ident.1 = "Proximal tubule_RCC", ident.2 = "Proximal tubule_CTRL", verbose = FALSE)
head(RCCCTRL3, n = 10)
write.csv(RCCCTRL3, "DEG-Proximal tubule.csv")

RCCCTRL4 <- FindMarkers(sobj, ident.1 = "Loop of Henle_RCC", ident.2 = "Loop of Henle_CTRL", verbose = FALSE)
head(RCCCTRL4, n = 10)
write.csv(RCCCTRL4, "DEG-Loop of Henle.csv")

RCCCTRL5 <- FindMarkers(sobj, ident.1 = "Renal interstitium_RCC", ident.2 = "Renal interstitium_CTRL", verbose = FALSE)
head(RCCCTRL5, n = 10)
write.csv(RCCCTRL5, "DEG-Renal interstitium.csv")

RCCCTRL6 <- FindMarkers(sobj, ident.1 = "Immune cells_RCC", ident.2 = "Immune cells_CTRL", verbose = FALSE)
head(RCCCTRL6, n = 10)
write.csv(RCCCTRL6, "DEG-Immune cells.csv")

RCCCTRL7 <- FindMarkers(sobj, ident.1 = "Cap mesenchyme_RCC", ident.2 = "Cap mesenchyme_CTRL", verbose = FALSE)
head(RCCCTRL7, n = 10)
write.csv(RCCCTRL7, "DEG-Cap mesenchyme.csv")

RCCCTRL8 <- FindMarkers(sobj, ident.1 = "α-intercalated cells_RCC", ident.2 = "α-intercalated cells_CTRL", verbose = FALSE)
head(RCCCTRL8, n = 10)
write.csv(RCCCTRL8, "DEG-α-intercalated.csv")

RCCCTRL9 <- FindMarkers(sobj, ident.1 = "Stromal cells_RCC", ident.2 = "Stromal cells_CTRL", verbose = FALSE)
head(RCCCTRL9, n = 10)
write.csv(RCCCTRL9, "DEG-Stromal cells.csv")
