# Specialized nextflow pipelines for the RRMM project
Contains analysis scripts and pipelines that should enable someone to reproduce the key steps nextflow pipeline.

# nextflow pipeline with 3 key steps:
1. Preprocessing, qc and cell type annotation (scrublet, SingleR)
   - Input: 10x mtx files and raw HCA data
   - Requirements: annotated cell-type reference: HCA BM (also healthy reference for tumor cell-type)
   - Established thresholds for qc-metrics and estimates for doublet fraction
2. Create a summary table and a merged object of the final processed file
3. Tumor clonality analysis with InferCNV
4. Cellular Interaction analysis with CellPhoneDB
