#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(tidyverse)
    library(future)
    library(Seurat)
    library(qs)
})
# Project specific Functions
options(future.globals.maxSize = 50 * 1024^3)
plan("multicore")
source("$workflow.projectDir/functions.R")

seurat_object <- qread("raw_seurat.qs")

# Plotting of certain QC metrics:
# General QC: nCount,nFeature,percent.mt
dir.create(file.path("QC_Plots"), showWarnings = F, recursive = T)
ungrouped_violins <- c(
    VlnPlot(object = seurat_object, features = c("nCount_RNA"), pt.size = 0, combine = F, group.by = "patient", log = T),
    VlnPlot(object = seurat_object, features = c("nFeature_RNA", "percent.mt"), pt.size = 0, combine = F, group.by = "patient", log = F)
)
ungrouped_violins <- purrr::map(ungrouped_violins, ~ . + stat_summary(fun.y = mean, geom = "point", size = 1.5, colour = "black") + NoLegend() + ggplot2::theme_minimal())

pdf(file.path("QC_Plots", "QC_Violin.pdf"))
CombinePlots(ungrouped_violins, ncol = 3, legend = "right")
dev.off()

# General QC but grouped by orig.ident
grouped_violins <- c(
    VlnPlot(object = seurat_object, features = c("nCount_RNA"), pt.size = 0, combine = F, group.by = "orig.ident", log = T),
    VlnPlot(object = seurat_object, features = c("nFeature_RNA", "percent.mt"), pt.size = 0, combine = F, group.by = "orig.ident", log = F)
)
grouped_violins <- purrr::map(grouped_violins, ~ . + stat_summary(fun.y = mean, geom = "point", size = 1.5, colour = "black") + NoLegend() + ggplot2::theme_minimal())

pdf(file.path("QC_Plots", "QC_Violin_grouped.pdf"))
CombinePlots(grouped_violins, ncol = 3, legend = "right")
dev.off()

# Scatter of QC-metrics
plot1 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf(file.path("QC_Plots", "QC_FeatureScatter.pdf"))
CombinePlots(plots = list(plot1, plot2), legend = "bottom")
dev.off()


# remove low quality cells
seurat_object <- subset(x = seurat_object, subset = nFeature_RNA > 400 & percent.mt < 10)


## Identify variable genes:
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_object), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat_object)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

pdf(file.path("QC_Plots", "QC_VariableFeatures.pdf"))
CombinePlots(plots = list(plot1, plot2), legend = "bottom")
dev.off()

## BLAS Control threads:
library(RhpcBLASctl)
blas_set_num_threads(1)
omp_set_num_threads(1)

## Run SCTransfrom with standard parameters
seurat_object <- suppressWarnings(SCTransform(seurat_object, verbose = TRUE, return.only.var.genes = FALSE))
 
# Run PCA and UMAP
seurat_object <- RunPCA(seurat_object, verbose = FALSE)
seurat_object <- RunUMAP(seurat_object, dims = 1:30, verbose = FALSE)

# Do a kNN and find clusters
seurat_object <- FindNeighbors(seurat_object, dims = 1:30, verbose = FALSE)
seurat_object <- FindClusters(seurat_object, verbose = FALSE)


# UMAP Plots of clusters and QC-metrics
umap1 <- DimPlot(seurat_object, label = TRUE) + NoLegend()
umap2 <- DimPlot(seurat_object, label = TRUE, group.by = "orig.ident") + NoLegend()
umap3 <- FeaturePlot(seurat_object, "nCount_RNA") + NoLegend()
umap4 <- FeaturePlot(seurat_object, "nFeature_RNA") + NoLegend()
pdf(file.path("QC_Plots", "QC_UMAP.pdf"), height = 8, width = 8)
CombinePlots(list(umap1, umap2, umap3, umap4))
dev.off()


# Plot marker genes - genes_of_interest in "functions.R"
ngenesplot(seurat_object, genes_of_interest, max = 9, file.path("QC_Plots"))


qsave(seurat_object, file.path(paste0("seurat_patient.qs")), nthread = 4)
