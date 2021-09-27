#!/usr/bin/env Rscript
# Setup R Session
library(Seurat)
library(qs)
library(tidyverse)
library(future)
library(RhpcBLASctl)
library(scrattch.io)

# Options for SCT
options(future.globals.maxSize = 50 * 1024^3)
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)
plan("multiprocess")


# Read seurat.qs file 
seuobj <- qread("patient.qs", nthreads = $task.cpus)

# Subset with Celltype_Cutoff
ct_cutoff <- $ct_cutoff
table(seuobj[["celltype_1"]])
seuobj<- subset(seuobj,subset = celltype_1 %in% names(which(table(seuobj[["celltype_1"]])>ct_cutoff)))

#SCTransform with Standard settings
seuobj <- SCTransform(seuobj, return.only.var.genes = F)

countmatrix <- GetAssayData(object = seuobj, assay = "SCT", slot = "data")



# Celltype Metadata
seu_meta <- FetchData(object = seuobj, vars = c("celltype_1"))

seu_meta <- tibble::rownames_to_column(seu_meta)
colnames(seu_meta) <- c("cell", "celltype")

# Writing out
write_dgCMatrix_csv(countmatrix, filename = "countmatrix_sct.csv", col1_name = "Gene_Name")
write_csv(seu_meta,"celltype_1.csv")
