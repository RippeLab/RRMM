# Subclone specific microenvironmental impact and drug responses in refractory multiple myeloma revealed by single cell transcriptomics

Nextflow pipelines and R code for scRNA-seq analysis of relapsed and refractory multiple myeloma samples
Folders contain scripts to create figures and a pipeline to process the raw data. 

## Figures:
Contains numbered R markdown notebooks for all figures in the manuscript.

## Processing-Pipeline:
Contains templates for the nextflow pipeline to:
1. Do quality Control and Cell-type annotation
2. Infer Clonalilty tumor cells with scRNAseq (InferCNV)
3. Predict cellular interactions with CellPhoneDB