#!/bin/bash
# Statistical Analysis
cellphonedb method statistical_analysis celltype_1.csv countmatrix_sct.csv \
--output-path="out" \
--iterations=1000 \
--threads $task.cpus \
--counts-data=gene_name

# Heatmap Plot
cellphonedb plot heatmap_plot  \
  celltype_1.csv 
