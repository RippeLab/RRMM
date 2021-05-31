#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(Seurat)
    library(qs)
    library(future)
    library(tidyverse)
})
plan("multicore")
options(future.globals.maxSize = 100 * 1024^3)

seurat_object <- qread("all_patients_merged.qs", nthread = 4)


library(RhpcBLASctl)
blas_set_num_threads(1)
omp_set_num_threads(1)

## Run SCTransfrom with standard parameters
seurat_object <- suppressWarnings(SCTransform(seurat_object, verbose = TRUE, return.only.var.genes = FALSE,	conserve.memory = TRUE))
 
 
# Run PCA and UMAP
seurat_object <- RunPCA(seurat_object, npcs = 70, verbose = TRUE)
seurat_object <- RunUMAP(seurat_object, dims = 1:70, verbose = TRUE)

# Do a kNN and find clusters
seurat_object <- FindNeighbors(seurat_object, dims = 1:70, verbose = TRUE)
seurat_object <- FindClusters(seurat_object, verbose = FALSE)

merged_plots_path <- "SCT_MergedPlots"

dir.create(merged_plots_path, showWarnings = F, recursive = T)

# UMAP Plots of clusters and QC-metrics
umap1 <- DimPlot(seurat_object, label = TRUE) + NoLegend()
umap2 <- DimPlot(seurat_object, label = TRUE, group.by = "orig.ident") + NoLegend()
umap3 <- FeaturePlot(seurat_object, "nCount_RNA") + NoLegend()
umap4 <- FeaturePlot(seurat_object, "nFeature_RNA") + NoLegend()
pdf(file.path(merged_plots_path, "QC_UMAP.pdf"), height = 8, width = 8)
CombinePlots(purrr::map(list(umap1, umap2, umap3, umap4),~AugmentPlot(.)))
dev.off()

pdf(file.path(merged_plots_path, "Patient_umap.pdf"), height = 8, width = 8)
AugmentPlot(DimPlot(seurat_object, label = TRUE, group.by = "patient") + NoLegend())
dev.off()

pdf(file.path(merged_plots_path, "Sorting_umap.pdf"), height = 8, width = 8)
AugmentPlot(DimPlot(seurat_object, label = TRUE, group.by = "sorting") + NoLegend())
dev.off()

pdf(file.path(merged_plots_path, "Sorting_hyperdiploidy.pdf"), height = 8, width = 8)
AugmentPlot(DimPlot(seurat_object, label = TRUE, group.by = "hyperdiploidy") + NoLegend())
dev.off()


pdf(file.path(merged_plots_path, "Celltype_1_Predicted_UMAP.pdf"), height = 7, width = 7)
AugmentPlot(DimPlot(seurat_object, group.by = "celltype_1_predicted", label = T) +
    NoLegend() +
    ggtitle(paste0(seurat_object[["orig.ident"]][[1]], " Assigned celltypes")))

dev.off()


pdf(file.path(merged_plots_path, "DoubletFinder_Doublet_UMAP.pdf"), height = 7, width = 7)
AugmentPlot(DimPlot(seurat_object, group.by = "celltype_1_predicted", cells.highlight = names(colnames(seurat_object)[which(is.na(seurat_object[["celltype_1_predicted"]][[1]]))]), label = T) +
    NoLegend() +
    ggtitle(paste0(seurat_object[["orig.ident"]][[1]], " Not assigned cells")))
dev.off()



umap1 <- DimPlot(seurat_object, label = TRUE) + NoLegend() +
    ggtitle(paste0(seurat_object[["orig.ident"]][[1]], " Cluster"))

umap_df <- DimPlot(seurat_object, group.by = "DoubletStatus", label = T) + NoLegend() +
    ggtitle(paste0(seurat_object[["orig.ident"]][[1]], " Detected Doublets"))

umap_df_pann <- FeaturePlot(seurat_object, "pANN") + NoLegend() +
    ggtitle(paste0(seurat_object[["orig.ident"]][[1]], " pANN Scores"))


pdf(file.path(merged_plots_path, "DoubletFinder_Doublet_UMAP.pdf"), height = 7, width = 14)
CombinePlots(list(umap1, umap_df))
dev.off()

pdf(file.path(merged_plots_path, "DoubletFinder_Doublet_UMAP_pANN.pdf"), height = 7, width = 7)
umap_df_pann
dev.off()



qsave(seurat_object, "all_patients_merged_sct_umap.qs", nthread = 4)
