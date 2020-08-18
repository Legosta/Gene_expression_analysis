if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("MAST", quietly = TRUE)) BiocManager::install("MAST")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("Matrix", quietly = TRUE)) install.packages("Matrix")
library(Seurat)
library(Matrix)
library(MAST)
library(ggplot2)
library(dplyr)


#upload the data

df <- Read10X("~/GSE114986")
dim(df)

#UMI_distribution

plotData <- data.frame(
  umis <- colSums(df)
)
ggplot(data=plotData, aes(x=umis)) +
  geom_histogram() + theme_bw()

#filtering genes and barcodes

seurat <- CreateSeuratObject(df, min.cells = 10, min.features = 10)
dim(seurat)

seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^mt:")
FeatureScatter(seurat, "nCount_RNA", "nFeature_RNA") + scale_x_log10() + scale_y_log10()

FeatureScatter(seurat, "nCount_RNA", "percent.mt") + scale_x_log10() 
FeatureScatter(seurat, "nFeature_RNA", "percent.mt") + scale_x_log10()

seurat <- subset(seurat, subset = nFeature_RNA > 800)
dim(seurat)    

#normalizaton

seurat <- SCTransform(seurat, vars.to.regress = "percent.mt", verbose = FALSE)

#PCA

seurat <- RunPCA(seurat, verbose = FALSE)
ElbowPlot(seurat, ndims = 50)

#UMAP

seurat <- RunUMAP(seurat, dims=1:30)
DimPlot(seurat, eduction = "umap") + NoLegend()

#clusters

seurat <- FindNeighbors(seurat, dims = 1:30, verbose = FALSE)
seurat <- FindClusters(seurat, resolution=0.6, verbose = FALSE)
DimPlot(seurat, reduction = "umap", label = TRUE) + NoLegend()

FeaturePlot(seurat, c("E93", "stg"), cols=c("grey", "red"), reduction="umap", ncol=3)
