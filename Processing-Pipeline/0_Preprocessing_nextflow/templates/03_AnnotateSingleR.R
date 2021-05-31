#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(SingleR)
    library(Seurat)
    library(qs)
    library(future)
    library(tidyverse)
})
# Project specific Functions
options(future.globals.maxSize = 50 * 1024^3)
#plan("multicore")

seurat_object <- qread("patient_seurat.qs")
hca_small <- qread("HCA_subsample.qs")
hca_counts <- hca_small@assays[["SCT"]]@data

patient_sct_counts <- seurat_object@assays[["SCT"]]@data



all_2_pred <- SingleR(  test = patient_sct_counts,
                        ref = hca_counts,
                        labels = hca_small[["celltype_1"]][[1]],
                        de.method="wilcox",
                        BPPARAM=MulticoreParam($task.cpus))


seurat_object <- AddMetaData(seurat_object, metadata = all_2_pred[["pruned.labels"]], col.name = "celltype_1_predicted")
scores <- data.frame(all_2_pred[["scores"]])
colnames(scores) <- paste0(colnames(scores),"_pred_score")
rownames(scores) <- colnames(seurat_object)

seurat_object <- AddMetaData(seurat_object, scores)

dir.create("SingleR_Plots", showWarnings = F, recursive = T)
pdf(file.path("SingleR_Plots", "Celltype_1_Predicted_UMAP.pdf"), height = 7, width = 7)
DimPlot(seurat_object, group.by = "celltype_1_predicted", label = T) + 
        NoLegend() + 
        ggtitle(paste0(seurat_object[["orig.ident"]][[1]], " Assigned celltypes"))


dev.off()


pdf(file.path("SingleR_Plots", "DoubletFinder_Doublet_UMAP.pdf"), height = 7, width = 7)
DimPlot(seurat_object, group.by = "celltype_1_predicted",cells.highlight = which(is.na(seurat_object[["celltype_1_predicted"]][[1]])), label = T) +
    NoLegend() +
    ggtitle(paste0(seurat_object[["orig.ident"]][[1]], " Not assigned cells"))
dev.off()

seurat_object <- subset(x = seurat_object, subset = Scrublet_Score < 0.4)


qsave(seurat_object, file.path(paste0("singler_seurat.qs")), nthread = 4)