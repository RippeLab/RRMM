#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(infercnv)
    library(Seurat)
    library(qs)
})

tumor_sample_test <- qs::qread("tumor.qs")
ref_input <- qread("HCA_PCs.qs")


## prepare input by merging tumor and referenc dataset
CNV_input <- merge(x = tumor_sample_test, y = ref_input)

# annotations (sample id of ref_input = "normal")
cellAnnotations <- data.frame(CNV_input[["sample_id"]]) # can be any other metadata slot


infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = GetAssayData(object = CNV_input, slot = "counts"),
    annotations_file = cellAnnotations,
    gene_order_file = "/Volumes/ag-rippe/NGS_Stephan/General/gene_lists/gencode_v28_CR_ref_gene_pos.txt",
    ref_group_names = c("normal")
)


## run InferCNV (in this case for 10x data, without prior clustering for cell annotations, no HMM and looks for subclusters/clones)
infercnv_obj <- infercnv::run(infercnv_obj,
    cutoff = 0.1,
    out_dir = "./InferCNV_out/",
    cluster_by_groups = F,
    denoise = T,
    HMM = F,
    no_prelim_plot = T,
    analysis_mode = "subclusters",
    num_threads = $task.cpus
)

qs::qsave(infercnv_obj,"infercnv.qs")