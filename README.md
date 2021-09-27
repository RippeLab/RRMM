# Nextflow data processing pipelines and R code for single cell RNA sequencing analysis of relapsed and refractory multiple myeloma

The RRMM (relapsed and refractory multiple myeloma) repository contains Nextflow data processing pipelines and R code for single cell RNA sequencing (scRNA-seq) analysis. The code was used in the study *"Subclone-specific microenvironmental impact and drug responses in refractory multiple myeloma revealed by single cell transcriptomics"* by Stephan M. Tirier, Jan-Philipp Mallm, Simon Steiger, Alexandra M. Poos, Mohamed H. S. Awwad, Nicola Giesen, Nicola Casiraghi, Hana Susak, Katharina Bauer, Anja Baumann, Lukas John, Anja Seckinger, Dirk Hose, Carsten Müller-Tidow, Hartmut Goldschmidt, Oliver Stegle, Michael Hundemer, Niels Weinhold, Marc S. Raab and Karsten Rippe.   
Plots for main and supplementary figures of the study can be created with the R-scripts provided in the RRMM project and the data available at accession number [GSE161801](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161801) of Gene Expression Omnibus . 

## Figures
Contains R markdown notebooks to recreate figures from the above study.
### RRMM samples and scRNA-seq dataset
- [HTML Preview](http://htmlpreview.github.io/?https://raw.githubusercontent.com/RippeLab/RRMM/main/Figures/K43R_code_1.nb.html)

### scRNA-seq analysis of RRMM tumor cells
- [HTML Preview](http://htmlpreview.github.io/?https://raw.githubusercontent.com/RippeLab/RRMM/main/Figures/K43R_code_2.nb.html)

### Treatment response of +1q subclones
- [HTML Preview](http://htmlpreview.github.io/?https://raw.githubusercontent.com/RippeLab/RRMM/main/Figures/K43R_code_3.nb.html)

### Analysis of cellular interaction of myeloma and BME cells
- [HTML Preview](http://htmlpreview.github.io/?https://raw.githubusercontent.com/RippeLab/RRMM/main/Figures/K43R_code_4.nb.html)
- [HTML Preview Interactions](http://htmlpreview.github.io/?https://raw.githubusercontent.com/RippeLab/RRMM/main/Figures/K43R_code_4_interaction.nb.html)

### T-cell heterogeneity in RRMM patient vs. healthy donor samples
- [HTML Preview](http://htmlpreview.github.io/?https://raw.githubusercontent.com/RippeLab/RRMM/main/Figures/K43R_code_5.nb.html)


### CD16+ monocyte heterogeneity in RRMM
- [HTML Preview](http://htmlpreview.github.io/?https://raw.githubusercontent.com/RippeLab/RRMM/main/Figures/K43R_code_6.nb.html)
- [HTML Preview interactions](http://htmlpreview.github.io/?https://raw.githubusercontent.com/RippeLab/RRMM/main/Figures/K43R_code_6_interaction.nb.html)

### BME changes in +1q RRMM
- [HTML Preview](http://htmlpreview.github.io/?https://raw.githubusercontent.com/RippeLab/RRMM/main/Figures/K43R_code_7.nb.html)

## Processing Pipeline
Contains templates for the Nextflow pipeline to conduct the following tasks:
- [Preprocessing, qc and cell type annotation with scrublet and SingleR](https://github.com/RippeLab/RRMM/tree/main/Processing-Pipeline/1_Preprocessing_nextflow)
- [Create a summary table and a merged object of the final processed file](https://github.com/RippeLab/RRMM/tree/main/Processing-Pipeline/2_SummaryStatistics)
- [Identification of tumor subclones from the scRNA-seq data with InferCNV](https://github.com/RippeLab/RRMM/tree/main/Processing-Pipeline/3_InferCNV_nextflow)
- [Prediction of interactions between different cell types with CellPhoneDB](https://github.com/RippeLab/RRMM/tree/main/Processing-Pipeline/4_CellPhoneDB_nextflow)
