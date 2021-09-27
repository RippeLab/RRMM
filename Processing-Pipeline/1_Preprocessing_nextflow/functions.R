createPatientConfig <- function(patient, file.path) {
  configuration_template <- list(
    "General" = list(
      "Patient" = patient,
      "Working_Directory" = file.path,
      "genome" = "hg38"
    ),
    "QC" = list("Under Development"),
    "Seurat" = list(
      "nPCs" = 50,
      "PCs_to_use" = 1:25,
      "Resolution" = 0.8,
      "vars_to_regress" = "percent.mt"
    ),
    "CellAnnotation" = list("not annotated")
  )
  yaml::write_yaml(configuration_template, file.path)
}

Ig_genes <- c("IGLC2", "IGLC3","IGLC4","IGLC5","IGLC6","IGLC7",
                "IGLV7-46", "IGLV6-57", "IGLV4-3","IGLV3-25","IGLV3-21", "IGLV3-19", "IGLV3-10", "IGLV3-1","IGLV2-8", 
                "IGLV1-36", "IGLV1-40", "IGLV1-44", "IGLV1-47", "IGLV1-51", "IGLV2-11","IGLV2-14", "IGLV2-23","IGLV2-28","IGLV3-7",
                "IGKC",
                "IGKV1-13","IGKV1-27", "IGKV1-33", "IGKV1-39", "IGKV1-5", "IGKV1-6","IGKV1D-33", "IGKV1D-39", "IGKV1D-8",
                "IGKV1OR2-108", "IGKV2-28", "IGKV2D-28", "IGKV3-11", "IGKV3-15","IGKV3-20", "IGKV3D-20", "IGKV4-1", 
                "IGHA1", "IGHA2", "IGHG1","IGHG2","IGHG3","IGHG4","IGHM", "JCHAIN")


mergeSeuratList <- function(pathlist = list(), output.location) {
  if (length(pathlist) == 1) {
    warning("List of Length 1 supplied, will output input file at new location")
    file.copy(pathlist[[1]], output.location)
    
  } else {
    
  
  rds_list <- lapply(pathlist, readRDS)
  if(is.null(names(pathlist))){
    names(rds_list) <- 1:length(rds_list)
  } else {
    names(rds_list) <- names(pathlist)
  }
  if(any(duplicated(unlist(lapply(rds_list,colnames))))){
    rds_merge <- Seurat:::merge.Seurat(rds_list[[1]],y=c(rds_list[[2:length(rds_list)]]),add.cell.ids = names(rds_list))
  } else {
    rds_merge <- merge(rds_list[[1]],y=c(rds_list[[2:length(rds_list)]]))  
  }
  
  
  invisible(return(rds_merge))
  }
}


genes_of_interest <- c("SDC1", # MM Plasma/TumorCells == "CD138"
                       "CD38", # MM Plasma/Tumor
                       "TNFRSF17", # "BCMA" - B cell maturation factor
                       # Myeloid Cells: Monocytes/Macrophages/DC
                       "FCGR3A", # monocytes "CD16" - mutually exclusive with CD14 / low
                       "CD14", # monocytes
                       "FCER1A", # Dendritic cells
                       "HBB", # Erythroblasts - Own lineage ! Erythrozytes filtered out in Preproc Lab
                       "CD3D", # T-Cells generally
                       "CD4", # major split t-cells (modulating t-cells)
                       "CD8A", # major split in t-Cells (effector-t-cells) 
                       "PRF1", # Perforin1
                       "GZMB", #
                       "GZMH", #
                       "GZMK", # Memory CD8
                       "NKG7", # NK
                       "NCAM1", # lowly expressed
                       "KLRB1", # CD161 - putative -- Th17 t-cell population
                       "CCR7", # naive CD4 T-cells
                       "IL7R", # mem T-Cells (CD4 or CD8)
                       "MS4A1", # = CD20 B-cell marker
                       "SPINK2", "MYC", "B2M"# broad progenitor marker
)


seurat_aware_preproc <-
  function(Seurat_Paths,
           configuration,
           output_location,
           aware = TRUE,
           workers = 8,
           rerunSCT = FALSE) {
    if (file.exists(output_location) && aware == TRUE) {
      patient_seu_pre <- readRDS(output_location)
      if (all.equal(patient_seu_pre@misc$configuration, configuration) ==
          TRUE) {
        return(patient_seu_pre)
        stop("Already processed")
        # } else {
        #   Seurat_Object <- patient_seu_pre
        #   if (Seurat_Object@misc$configuration$vars_to_regress != configuration$Seurat$vars_to_regress ||
        #       is.null(Seurat_Object@misc$progress$sct)) {
        #     plan(strategy = "multicore", workers = workers)
        #     patient_seu <- suppressWarnings(SCTransform(
        #       Seurat_Object,
        #       vars.to.regress = configuration$Seurat$vars_to_regress,
        #       verbose = F
        #     ))
        #
        #   }
        
      } else {
        Seurat_Object <- patient_seu_pre
      }
      
    } else {
      Seurat_Object <- mergeSeuratList(Seurat_Paths)
    }
    
    
    
    
    Seurat_Object@misc$configuration <- configuration
    if (is.null(Seurat_Object@misc$sct_run) || rerunSCT) {
      print("Running SCTransform")
      plan(strategy = "multicore", workers = workers)
      Seurat_Object <- suppressWarnings(
        SCTransform(
          Seurat_Object,
          vars.to.regress = Seurat_Object@misc$configuration$Seurat$vars_to_regress,
          verbose = F
        )
      )
      Seurat_Object@misc$sct_run <- TRUE
    }
    
    s.genes <- cc.genes$s.genes
    g2m.genes <- cc.genes$g2m.genes
    Seurat_Object <- CellCycleScoring(
      object = Seurat_Object,
      s.features = s.genes,
      g2m.features = g2m.genes,
      set.ident = F
    )
    print("Running PCA")
    Seurat_Object <- RunPCA(
      object = Seurat_Object,
      npcs =  Seurat_Object@misc$configuration$Seurat$nPCs,
      verbose = FALSE
    )
    
    print("Running FindNeighbors")
    Seurat_Object <- FindNeighbors(
      object = Seurat_Object,
      reduction = "pca",
      dims =  Seurat_Object@misc$configuration$Seurat$PCs_to_use
    )
    Seurat_Object <-
      FindClusters(
        Seurat_Object,
        resolution =  Seurat_Object@misc$configuration$Seurat$Resolution,
        verbose = F
      )
    
    print("Running RunUMAP")
    Seurat_Object <- RunUMAP(
      object = Seurat_Object,
      reduction = "pca",
      dims = Seurat_Object@misc$configuration$Seurat$PCs_to_use,
      verbose = F
    )
    
    saveRDS(Seurat_Object, file = output_location)
    invisible(return(Seurat_Object))
  }

#saveRDS(rds_merge, file = output.location)
ngenesplot <- function(patient_seu,genes_of_interest, max = 16, plot.path) {
  d <- seq_along(genes_of_interest)
  genelist <- split(genes_of_interest, ceiling(d / max))
  
  
  lapply(genelist, function(genes) {
    FeaturePlot(patient_seu,
                features = genes,
                cols = c("grey", "firebrick"))
    ggsave(
      paste0(plot.path, "/Markers_", paste(genes, collapse = "_"), ".png"),
      width = 15,
      height = 15
    )
  })
}



clustermarkers <- function(patient_seu,presto_topmarker_auc){
  purrr::map(levels(as.factor(presto_topmarker_auc$group)),function(cluster){
    
    
    genes <- presto_topmarker_auc %>% filter(group==cluster) %>% ungroup() %>% select(feature)  %>% unlist()
    cells <- names(patient_seu$seurat_clusters[patient_seu$seurat_clusters ==cluster])
    testdim <- DimPlot(patient_seu,cells.highlight = cells,combine =F,label=T) %>% 
      purrr::map(~ . + ggtitle(paste("Cluster",cluster)))
    test <- FeaturePlot(patient_seu,genes,combine = F,min.cutoff = 0,max.cutoff = 6,cols = viridis::magma(100))
    
    #test[[1]] + theme_bw()
    CombinePlots(purrr::map(
      c(testdim, test),
      ~ . + theme_minimal() + theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      ) + NoLegend() 
    ))
    dir.create(file.path(params$patient,"Plots/DE"),showWarnings = F)
    ggsave(file.path(params$patient,paste0("Plots/","DE/",
                                           "Markers_Cluster_",sprintf("%02d",as.numeric(cluster)),".png")),
           width = 12,
           height = 12)
    
  })
  
}

scrollingmarkers <- function(presto_topmarker_fc){
  knitr::kable(presto_topmarker_fc)  %>% 
    kable_styling() %>%
    scroll_box(height = "400px")}



clusterde <- function(patient_hca,celltype_sel,deAll){
  metadata <- patient_hca@meta.data
  metadata$cellname <- names(patient_hca$orig.ident)
  metadata <- metadata %>% dplyr::filter(topcelltype == celltype_sel)
  patient_hca_subset <- subset(patient_hca,cells = metadata$cellname)
  
  deSUBSET <-
    wilcoxauc(patient_hca_subset,
              group_by = "patient",
              seurat_assay = "SCT",assay = "data")
  
  auc_diff <-
    deSUBSET %>% left_join(deAll, by = c("feature", "group")) %>% mutate(auc_diff =
                                                                           auc.x - auc.y) %>% arrange(desc(auc_diff))
  
  
  
  
  invisible(return(auc_diff))
}


aucdif_plot <- function(.data,name) {
  ggplot(.data ,
         aes(
           x = auc.x,
           y = auc.y,
           color = auc_diff,
           label = ifelse(auc_diff > max(auc_diff) * 0.7 |
                            auc.x > 0.9, feature, NA)
         )) + geom_point() + ggrepel::geom_label_repel(nudge_x = 0.07) + coord_flip() + scale_color_viridis_c() + ggtitle(name)
  
}
