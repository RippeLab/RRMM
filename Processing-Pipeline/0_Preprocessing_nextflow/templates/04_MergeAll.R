#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(Seurat)
    library(qs)
    library(future)
    library(tidyverse)
})
plan("multicore")


seurat_file_list <- list.files()
seu_patient_list <- sapply(seurat_file_list,function(seurat_file){
    seurat_object <- qread(seurat_file)
    print(paste("Read in patient",seurat_object[["patient"]][[1]][[1]]))
    seurat_object
})

seurat_object <- merge(seu_patient_list[[1]], y = seu_patient_list[c(2:length(seu_patient_list))], merge.data = TRUE)


qsave(seurat_object, "all_patients_merged.qs", nthread = 4)
