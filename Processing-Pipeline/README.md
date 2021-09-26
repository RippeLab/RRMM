# Nextflow pipelines developed for the RRMM project
This folder contains Nextflow pipelines and analysis scripts to conduct the following main steps:

1. Preprocessing, qualtiy control (qc) and cell type annotation with scrublet and SingleR
   - Input: 10x mtx files and raw HCA data
   - Requirements: annotated cell-type reference: HCA BM (also healthy reference for tumor cell-type)
   - Established thresholds for qc-metrics and estimates for doublet fraction
2. Create a summary table and a merged object of the final processed file
3. Identification of tumor subclones from the scRNA-seq data with InferCNV
4. Prediction of interactions between different cell types with CellPhoneDB
