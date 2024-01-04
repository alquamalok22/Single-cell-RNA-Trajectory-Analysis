---
title: "Single rna seq"
author: "Alquama"
date: "01/04/2024"
output: html_document
---

```{r setup, include=FALSE}
```

```{r}
set.seed(1234)

library(monocle3)
library(SeuratWrappers)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(Matrix)
```

## Read the data
```{r}

NML_1 = Read10X(data.dir = "../Desktop/Bioinfo Projects/Monocle 3 data/NML1/")
NML_2 = Read10X(data.dir = "../Desktop/Bioinfo Projects/Monocle 3 data/NML2/")
NML_3 = Read10X(data.dir = "../Desktop/Bioinfo Projects/Monocle 3 data/NML3/")

## Seurat object ##
NML_1 = CreateSeuratObject(counts = NML_1, project = "NML_1",min.cells = 3,min.features = 200)
NML_1 = PercentageFeatureSet(NML_1, pattern = "^MT-",col.name = "percent.mt")

NML_2 = CreateSeuratObject(counts = NML_2, project = "NML_2",min.cells = 3,min.features = 200)
NML_2 = PercentageFeatureSet(NML_2, pattern = "^MT-",col.name = "percent.mt")

NML_3 = CreateSeuratObject(counts = NML_3, project = "NML_3",min.cells = 3,min.features = 200)
NML_3 = PercentageFeatureSet(NML_3, pattern = "^MT-",col.name = "percent.mt")

NML <- merge(NML_1,y=c(NML_2,NML_3))
VlnPlot(NML, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

NML <- subset(NML,nFeature_RNA < 5000 & nFeature_RNA < 2000 & 
                   percent.mt < 10)
View(NML@meta.data)

NML.list <- SplitObject(NML, split.by = "orig.ident")
## Normalisation part ##
NML.list <- lapply(X = NML.list, FUN = function(X) {
 X <- NormalizeData(X)
 X <- FindVariableFeatures(X, selection.method = "vst" , nfeatures = 2000)
})

## Features ##
features <- SelectIntegrationFeatures(object.list = NML.list)
features <- SelectIntegrationFeatures(object.list = NML.list)
NML.anchors <- FindIntegrationAnchors(object.list = NML.list, anchor.features = features)

## Integrated data assay ##
NML.combined <- IntegrateData(anchorset = NML.anchors)

## Downstream Analysis ##
DefaultAssay(NML.combined) <- "integrated"

## Cell clustering ##
NML.combined <- ScaleData(NML.combined, verbose = FALSE)
NML.combined <- RunPCA(NML.combined, npcs =30 , verbose = FALSE)
NML.combined <- RunUMAP(NML.combined, reduction = "pca" , dims = 1:30)
NML.combined <- FindNeighbors(NML.combined, reduction = "pca" , dims = 1:30)
NML.combined <- FindClusters(NML.combined, resolution = 0.5)

## Visualisation ##
DimPlot(NML.combined, reduction = "umap" , label = "TRUE")

# Four Cluster Markers ###
FeaturePlot(NML.combined, features = c("COL1A2" , "LUM" , "PDGFRA" , "PDGFRB"),
            label = T , cols = c('lightgrey', 'blue'))
FeaturePlot(NML.combined, features = c("COL1A2" , "LUM" , "PDGFRA" , "PDGFRB"),
             label = T , cols = c('lightgrey', 'blue'))
FeaturePlot(NML.combined, features = c("EPCAM" , "SFTPC" , "AGER" , "SCGB1A1"),
             label = T , cols = c('lightgrey', 'blue'))
FeaturePlot(NML.combined, features = c("CLDN5" , "CCL21" , "PECAM1" , "EMCN"),
             label = T , cols = c('lightgrey', 'blue'))
FeaturePlot(NML.combined, features = c("PTPRC" , "CD52" , "AIF1" , "TRBC2"),
             label = T , cols = c('lightgrey', 'blue'))
FeaturePlot(NML.combined, features = c("MSLN" , "CALB2" , "HP" , "PRG4"),
             label = T , cols = c('lightgrey', 'blue'))

NML_Mesenchymal <- subset(x = NML.combined, idents = c("0","1","2","3" ,"4" ,"5" ,"6","7" ,"8" ,"9"))
DimPlot(NML_Mesenchymal, reduction = "umap" , label = "TRUE")


## Standard flow for the clustering ##
NML_Mesenchymal <- ScaleData(NML_Mesenchymal, verbose = FALSE)
NML_Mesenchymal <- RunPCA(NML_Mesenchymal, npcs =30 , verbose = FALSE)
NML_Mesenchymal <- RunUMAP(NML_Mesenchymal, reduction = "pca" , dims = 1:30)
NML_Mesenchymal <- FindNeighbors(NML_Mesenchymal, reduction = "pca" , dims = 1:30)
NML_Mesenchymal <- FindClusters(NML_Mesenchymal, resolution = 0.5)
DimPlot(NML_Mesenchymal, reduction = "umap" , label = "TRUE")

NML_Mesenchymal <- subset(x = NML_Mesenchymal, idents = c("0","1","2","3" ,"4" ,"5" ,"6"))
FeaturePlot(NML_Mesenchymal, features = c("ACTA2" ,"RGS5" , "LUM" , "ACTG2"),
            label = T , cols = c('lightgrey', 'blue'))
FeaturePlot(NML.combined, features = c("COL1A1" , "CTHRC1" , "COL3A1" , "POSTN"),
             label = T , cols = c('lightgrey', 'blue'))

NML_Mesenchymal <- RenameIdents(NML_Mesenchymal , '0' = "Pericytes" , '4'= "Pericytes", '1'= "SMCs" , '3'= "MyoF" , '2'= "AlvF", '5'= "AlvF" , '6'= "AirwayF")
DimPlot(NML_Mesenchymal, reduction = "umap" , label = "TRUE")
saveRDS(NML_Mesenchymal, file ="../Desktop/Bioinfo Projects/Monocle 3 data/NML_Mesenchymal.RDS")
#rm("features" , "NML.anchors", "NML.combined" , "NML" ,"NML.list","NML_1" ,"NML_2" ,"NML_3")
```


2. Convert Seruat object into cell data set ....
```{r}
cds <- as.cell_data_set(NML_Mesenchymal, group.by = 'ident')
DefaultAssay(NML_Mesenchymal)
DefaultAssay(NML_Mesenchymal) <- "RNA"
colData(cds)
#plot#
plot_cells(cds, color_cells_by = 'seurat_clusters', label_groups_by_cluster = FALSE,group_label_size = 5) + theme(legend.position = "right")

plot_cells(cds, color_cells_by = 'ident', label_groups_by_cluster = FALSE,group_label_size = 5) + theme(legend.position = "right")

### trajectory graph ###
cds <- learn_graph(cds, use_partition = FALSE)
cds <- cluster_cells(cds, resolution = 1e-5)
plot_cells(cds, color_cells_by = 'ident', label_groups_by_cluster = FALSE,group_label_size = 5) + theme(legend.position = "right")
```

3. Trajectory Analysis...

```{r}
cds <- learn_graph(cds, use_partition = FALSE)
plot_cells(cds, color_cells_by = 'ident', label_groups_by_cluster = FALSE,group_label_size = 5) + theme(legend.position = "right")
plot_cells(cds, genes = c("CD34" , "PI16" , "SCARA5" , "MFAP5"))

## Gene Annotation file ##
fData(cds)
rownames(fData(cds))[1:10]
fData(cds)$gene_short_name<-rownames(fData(cds))
fData(cds)

## Order the cells by pseudotime
cds <- order_cells(cds)
genes_to_check <- c("CD34", "PI16", "SCARA5", "MFAP5")
genes_in_cds <- colData(cds)$gene_name  # Adjust the column name if needed
missing_genes <- genes_to_check[!(genes_to_check %in% genes_in_cds)]

if (length(missing_genes) > 0) {
    cat("The following genes are missing from the cds object:", missing_genes, "\n")
 }

## The following genes are missing from the cds object: CD34 PI16 SCARA5 MFAP#
```

