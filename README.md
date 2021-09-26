# Nextflow data processing pipelines and R code for single cell RNA sequencing analysis of relapsed and refractory multiple myeloma

The RRMM (relapsed and refractory multiple myeloma)repository contains Nextflow data processing pipelines and R code for single cell RNA sequencing (scRNA-seq) analysis. The codes was used in the study "Subclone-specific microenvironmental impact and drug responses in refractory multiple myeloma revealed by single cell transcriptomics" by Stephan M. Tirier, Jan-Philipp Mallm, Simon Steiger, Alexandra M. Poos, Mohamed H. S. Awwad, Nicola Giesen, Nicola Casiraghi, Hana Susak, Katharina Bauer, Anja Baumann, Lukas John, Anja Seckinger, Dirk Hose, Carsten Müller-Tidow, Hartmut Goldschmidt, Oliver Stegle, Michael Hundemer, Niels Weinhold, Marc S. Raab and Karsten Rippe. Plots for figures 1-7 of the study can be created with the R-scripts provided in the RRMM project and the data available at accession number GSE161801 of Gene Expression Omnibus (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161801). 

## Figures
Contains numbered R markdown notebooks for figures 1-7 from the above study.

### Figure 1
- [R Notebook Figure 1](http://htmlpreview.github.io/?https://raw.githubusercontent.com/RippeLab/RRMM/main/Figures/K43R_code_Fig_1.nb.html)

### Figure 2
- [R Notebook Figure 2](http://htmlpreview.github.io/?https://raw.githubusercontent.com/RippeLab/RRMM/main/Figures/K43R_code_Fig_2.nb.html)

### Figure 3
- [R Notebook Figure 3](http://htmlpreview.github.io/?https://raw.githubusercontent.com/RippeLab/RRMM/main/Figures/K43R_code_Fig_3.nb.html)

### Figure 4
- [R Notebook Figure 4](http://htmlpreview.github.io/?https://raw.githubusercontent.com/RippeLab/RRMM/main/Figures/K43R_code_Fig_4.nb.html)
- [R Notebook Figure 4 Interactions](http://htmlpreview.github.io/?https://raw.githubusercontent.com/RippeLab/RRMM/main/Figures/K43R_code_Fig_4_interaction.nb.html)

### Figure 5
- [R Notebook Figure 5](http://htmlpreview.github.io/?https://raw.githubusercontent.com/RippeLab/RRMM/main/Figures/K43R_code_Fig_5.nb.html)

### Figure 6
- [R Notebook Figure 6](http://htmlpreview.github.io/?https://raw.githubusercontent.com/RippeLab/RRMM/main/Figures/K43R_code_Fig_6.nb.html)
- [R Notebook Figure 6 interactions](http://htmlpreview.github.io/?https://raw.githubusercontent.com/RippeLab/RRMM/main/Figures/K43R_code_Fig_6_interaction.nb.html)

### Figure 7
- [R Notebook Figure 7](http://htmlpreview.github.io/?https://raw.githubusercontent.com/RippeLab/RRMM/main/Figures/K43R_code_Fig_7.nb.html)

## Processing Pipeline
Contains templates for the Nextflow pipeline to conduct the following tasks:
- [Preprocessing, qc and cell type annotation with scrublet and SingleR](https://github.com/RippeLab/RRMM/tree/main/Processing-Pipeline/0_Preprocessing_nextflow)
- [Create a summary table and a merged object of the final processed file](https://github.com/RippeLab/RRMM/tree/main/Processing-Pipeline/1_SummaryStatistics)
- [Identification of tumor subclones from the scRNA-seq data with InferCNV](https://github.com/RippeLab/RRMM/tree/main/Processing-Pipeline/2_InferCNV_nextflow)
- [Prediction of interactions between different cell types with CellPhoneDB](https://github.com/RippeLab/RRMM/tree/main/Processing-Pipeline/3_CellPhoneDB_nextflow)
