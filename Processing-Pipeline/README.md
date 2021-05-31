# Specialized nextflow Pipelines for the RRMM project
Containing analysis scripts and pipelines for the analysis.  
Should enable someone to reproduce the key steps nextflow pipeline.

# nextflow pipeline with 3 key steps:
1. preprocessing and qc & Cell type annotation  (scrublet,SingleR)
	- input: 10x mtx files + raw HCA data
	- requirements: annotated cell-type reference: HCA BM (also healthy reference for tumor cell-type)
	- Established thresholds for QC-metrics and estimates for doublet fraction
	- 2. create a summary table and a merged object of the final processed file
	
3. Tumor clonality analysis with InferCNV
4. Cellular Interaction analysis with CellPhoneDB
