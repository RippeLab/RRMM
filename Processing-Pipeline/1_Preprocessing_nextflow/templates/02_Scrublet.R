#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(tidyverse)
    library(future)
    library(Seurat)
    library(qs)
    library(DoubletFinder)
})
# Project specific Functions

seurat_object <- qread("patient_seurat.qs")

# Run Scrublet on the seurat object: Creates temporary file output.

mat <- GetAssayData(object = seurat_object, assay = "RNA", slot = "counts")
# write out a mtx file
tf <- "count_matrix"
dtf <- paste(tf, "doubletScores", sep = ".")
Matrix::writeMM(mat, file=tf)
cmd <- paste0("python -c 'import sys; import pandas; import scrublet;import scipy.io; df = scipy.io.mmread(\"", tf, "\").T.tocsc(); scrub = scrublet.Scrublet(df,expected_doublet_rate=0.10); doublet_scores, predicted_doublets = scrub.scrub_doublets(); pandas.DataFrame(doublet_scores).to_csv(\"", dtf, "\");'")
tmp <- system(cmd, intern = T)

x <- as.numeric(as.data.frame(data.table::fread(dtf, sep = ",", header = F, skip = 1))[, 2])
names(x) <- colnames(mat)
file.remove(tf)
file.remove(dtf)

seurat_object <- AddMetaData(seurat_object, x, "Scrublet_Score")

umap_scr_score <- FeaturePlot(seurat_object, "Scrublet_Score") + NoLegend() +
    ggtitle(paste0(seurat_object[["orig.ident"]][[1]], " Scrublet Scores"))

dir.create("Scrublet_Plots",showWarnings=F,recursive=T)
pdf(file.path("Scrublet_Plots", "Scrublet_Doublet_UMAP.pdf"), height = 7, width = 7)
umap_scr_score
dev.off()

#readr::write_csv(seurat_object@meta.data %>% dplyr::select(c("DoubletStatus","pANN")), "DF_DoubletStatus.txt")


qsave(seurat_object,file.path(paste0("scrublet_seurat.qs")),nthread=4)