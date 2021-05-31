#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(tidyverse)
    library(Seurat)
    library(qs)
})
# Project specific Functions
future::plan("multicore")
patient_annotation <- readxl::read_excel("/media/ag-rippe/NGS_Stephan/HIPO2_K43R/Analysis_v3/_HIPO_K43R_Patient_Metadata.xlsx")


genome <- "$genome"
patient_10x <- "$patientkey"

#Filter patient_annotation for this patient_id
metadata <- patient_annotation %>% filter(orig.ident == patient_10x)


print(patient_10x)
pat_10x_data <- Read10X(file.path(patient_10x, genome))
# Rename cell names
colnames(pat_10x_data) <- paste0(patient_10x, "_", colnames(pat_10x_data))

# Get and store counts of Ig genes in independent file
Ig_genes <- read.table("Ig_genes.csv", col.names = FALSE)[, 1]
pat_10x_data_Ig <- pat_10x_data[Ig_genes, ]
write.csv(pat_10x_data_Ig, file = file.path(paste0(patient_10x, "_ig_genes.csv")))


# Remove Ig genes from count matrix
pat_10x_data_noIg <- pat_10x_data[!row.names(pat_10x_data) %in% Ig_genes, ]


namesfield <- seq_len(length(unlist(str_split(patient_10x, "_"))))
# Create Seurat Object

seurat_object <- CreateSeuratObject(
    counts = pat_10x_data_noIg,
    min.cells = 0,
    min.features = 200,
    project = patient_10x,
    names.field = namesfield
)
 
# Add Metadata per Patient
seurat_object <- AddMetaData(seurat_object, as.list(metadata))

# Add IGG-Counts:
seurat_object <- AddMetaData(object = seurat_object, metadata = as.data.frame(Matrix::t(pat_10x_data_Ig)))


# Calculate percentage of mitochondrial reads based on genes starting with MT-
seurat_object[["percent.mt"]] <- PercentageFeatureSet(object = seurat_object, pattern = "^MT-")

write(ncol(seurat_object), "cellno.txt")

qsave(seurat_object, file.path(paste0("raw_seurat.qs")), nthread = 4)