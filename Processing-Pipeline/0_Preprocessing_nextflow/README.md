# Preprocessing workflow

## Input: 
- 10x_Cellranger matrix files for all samples

## Steps:
 - Create Seurat Objects
 - Remove and Filter genes/cells
 - Preproc and Cluster by patient
 - Doublet removal (scrublet)
 - Annotation with SingleR
 - Merge and Process Merged object

## Output:  
Folders 01-04 by steps as depicted in the template folders (04-06 in folder 04_Merged)
Resulting in two merged datasets 
 - all_patients_merged.qs (no SCT, No Embedding)
 - all_patients_merged_sct_umap.qs (SCT + UMAP with all patients)

