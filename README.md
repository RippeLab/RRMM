# Subclone specific microenvironmental impact and drug responses in refractory multiple myeloma revealed by single cell transcriptomics

Nextflow pipelines and R code for scRNA-seq analysis of relapsed and refractory multiple myeloma samples
Folders contain scripts to create figures and a pipeline to process the raw data. 

## Figures:
Contains numbered R markdown notebooks for all figures in the manuscript.
### Figure 1
- [HTML Preview Figure 1](http://htmlpreview.github.io/?https://raw.githubusercontent.com/RippeLab/RRMM/main/Figures/K43R_code_Fig_1.nb.html)

### Figure 2
- [HTML Preview Figure 2](http://htmlpreview.github.io/?https://raw.githubusercontent.com/RippeLab/RRMM/main/Figures/K43R_code_Fig_2.nb.html)

### Figure 3
- [HTML Preview Figure 3](http://htmlpreview.github.io/?https://raw.githubusercontent.com/RippeLab/RRMM/main/Figures/K43R_code_Fig_3.nb.html)

### Figure 4
- [HTML Preview Figure 4](http://htmlpreview.github.io/?https://raw.githubusercontent.com/RippeLab/RRMM/main/Figures/K43R_code_Fig_4.nb.html)
- [HTML Preview Figure 4 Interactions](http://htmlpreview.github.io/?https://raw.githubusercontent.com/RippeLab/RRMM/main/Figures/K43R_code_Fig_4_interaction.nb.html)

### Figure 5
- [HTML Preview Figure 5](http://htmlpreview.github.io/?https://raw.githubusercontent.com/RippeLab/RRMM/main/Figures/K43R_code_Fig_5.nb.html)


### Figure 6
- [HTML Preview Figure 6](http://htmlpreview.github.io/?https://raw.githubusercontent.com/RippeLab/RRMM/main/Figures/K43R_code_Fig_6.nb.html)
- [HTML Preview Figure 6 interactions](http://htmlpreview.github.io/?https://raw.githubusercontent.com/RippeLab/RRMM/main/Figures/K43R_code_Fig_6_interaction.nb.html)

### Figure 7
- [HTML Preview Figure 7](http://htmlpreview.github.io/?https://raw.githubusercontent.com/RippeLab/RRMM/main/Figures/K43R_code_Fig_7.nb.html)


## Processing-Pipeline:
Contains templates for the nextflow pipeline to:
1. Do quality Control and Cell-type annotation
2. Infer Clonalilty tumor cells with scRNAseq (InferCNV)
3. Predict cellular interactions with CellPhoneDB
