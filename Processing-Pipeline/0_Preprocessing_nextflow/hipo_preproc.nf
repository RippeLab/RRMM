//patients = ["VP6CSY_N1"]
//patients = ["VP6CSY_1","3S3UDC_2","37A66E_1","6KLSUS_2"]

def memory_fun(String cellno,attempt) {
        if(cellno.toInteger()>15000)  
            300.GB
        else if(cellno.toInteger()>10000)  
            200.GB
        else if(cellno.toInteger()>6000) 
            150.GB * attempt
        else if(cellno.toInteger()>3500) 
            120.GB * attempt
        else if (cellno.toInteger()>3000) 
            90.GB * attempt
        else if (cellno.toInteger()>1000) 
            70.GB * attempt
        else if (cellno.toInteger()>700) 
            40.GB * attempt
        else
            30.GB * attempt
}

// Create input channel from folder 00_input with all patients
input = Channel
    .fromPath("00_input/*/*", type:"dir")
    .map { file ->
        def key = file.name.toString()
        return tuple(key,"$file")//.collect { "$it/hg38" }.join(' '))
    }
    .groupTuple()
    .filter { if(binding.variables["patients"]){ patients.contains(it[0])} else { true } }


process create_seurat {
    cache 'deep'
    tag "cs_$patientkey"
    publishDir "01_Preparation/$patientkey"
    label 'rsing'

    cpus 1
    memory 12.GB

    input:
    val genome from "hg38"
    set patientkey, path(input_path) from input
    path ig_genes from "/media/ag-rippe/NGS_Stephan/General/gene_lists/Ig_genes.csv"

    output:
    stdout creation 
    tuple val(patientkey), path("raw_seurat.qs"), file("cellno.txt") into raw_seurat
    file "*_ig_genes.csv"

    script:
    template '00_CreateSeurat.R'
}

// --- # Preprocessing and first QC # ---
process preproc_seurat {
    cache 'deep'
    tag "sp_$patientkey"
    publishDir "01_Preparation/$patientkey"
    label 'rsing'

    cpus 6
    memory { memory_fun(cellno.text,task.attempt) }

    errorStrategy { task.exitStatus in 137..143 ? 'retry' : 'terminate' }
    maxRetries 3


    input:
    tuple val(patientkey), file("raw_seurat.qs"), cellno from raw_seurat

    output:
    stdout preproc 
    tuple val(patientkey), file("seurat_patient.qs"), cellno into seurat
    path "QC_Plots/*"


    script:
    template '01_ProcessSeurat.R'
}


// --- # DoubletFinder # ---
process scrublet_seurat {
    cache 'lenient'
    validExitStatus 0,1,2 // Might fail?
    tag "scrub_$patientkey"
    conda '/home/bq_ssteiger/miniconda3/envs/scrublet/'
    publishDir "02_Doublets/$patientkey"
    //label 'rsing'
    
    cpus 2
    memory { memory_fun(cellno.text,task.attempt) }
    errorStrategy { task.exitStatus in 137..143 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
    tuple val(patientkey), file("patient_seurat.qs"), cellno from seurat

    output:
    stdout scrublet 
    tuple val(patientkey), file("scrublet_seurat.qs"),cellno into seurat_dfinder
    path "Scrublet_Plots/*"

    script:
    template '02_Scrublet.R'
}



// --- # SingleR Celltype Assignment # ---
process singler_seurat {
    cache 'lenient'
    validExitStatus 0,1,2 // Might fail?
    tag "singler_$patientkey"
    publishDir "03_Annotation/$patientkey", mode: 'copy'
    label 'rsing'
    
    cpus 2
    memory { memory_fun(cellno.text,task.attempt) }
    errorStrategy { task.exitStatus in 137..143 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
    tuple val(patientkey), file("patient_seurat.qs"), cellno from seurat_dfinder
    each file("HCA_subsample.qs") from file("HCA/HCA_subsample.qs")
    output:
    stdout singler
    
    tuple val(patientkey), file("singler_seurat.qs"), cellno into seurat_singler
    path "SingleR_Plots/*"

    script:
    template '03_AnnotateSingleR.R'
}

seurat_singler.multiMap{
    it -> 
    file: it[1]
    sample: it[0]
}.set { seurat_merge }

//seurat_merge.file.collect().subscribe{println it }

// --- # Merge All Seurat Objects # ---
process merge_all {
    cache 'deep'
    validExitStatus 0,1,2 // Might fail?
    tag "MergeAll"
    publishDir "04_Merged/", mode: 'copy'
    label 'rsing'
    
    cpus 4
    memory 370.GB
    time '24h'
    //errorStrategy { task.exitStatus in 137..143 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
    
    file 'patient_seurat_??.qs' from seurat_merge.file.collect()

    output:
    stdout merge
    file("all_patients_merged.qs") into merged_seurat

    script:
    template '04_MergeAll.R'
}

merged_seurat
    .into{ merge_independent_proc; merge_sct_proc }


// --- # Process Merged Object # ---
process merge_proc {
    cache 'lenient'
    validExitStatus 0,1,2 // Might fail?
    tag "MergeAll"
    publishDir "04_Merged/", mode: 'copy'
    label 'rsing'
    
    cpus 4
    memory 370.GB
    time '24h'
    //errorStrategy { task.exitStatus in 137..143 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
    
    file 'all_patients_merged.qs' from merge_independent_proc

    output:
    stdout mergeproc
    file("all_patients_merged_umap.qs") into processed_merge
    path "MergedPlots/*"

    script:
    template '05_ProcMerge.R'
}




// --- # Process Merged Object # ---
// process merge_sct_proc {
//     cache 'lenient'
//     validExitStatus 0,1,2 // Might fail?
//     tag "MergeSCT"
//     publishDir "04_Merged/", mode: 'copy'
//     label 'rsing'
    
//     cpus 1
//     memory 370.GB
//     time '24h'
//     //errorStrategy { task.exitStatus in 137..143 ? 'retry' : 'terminate' }
//     maxRetries 3

//     input:
    
//     file 'all_patients_merged.qs' from merge_sct_proc

//     output:
//     stdout mergesctproc
//     file("all_patients_merged_sct_umap.qs") into processed_sct_merge
//     path "SCT_MergedPlots/*"

//     script:
//     template '06_SCTMerge.R'
// }

merge.subscribe{println it }
