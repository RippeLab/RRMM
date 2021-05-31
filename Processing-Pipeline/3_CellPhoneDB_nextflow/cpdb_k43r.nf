def memory_fun(String cellno,attempt) {
        if(cellno.toInteger()>15000)  
            350.GB
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
            50.GB * attempt
        else
            30.GB * attempt
}

def core_fun(String cellno,attempt) {
        if(cellno.toInteger()>15000)  
            10
        else if(cellno.toInteger()>10000)  
            Math.floor( 15 / attempt ) 
        else if(cellno.toInteger()>6000) 
            Math.floor( 20 / attempt )
        else if(cellno.toInteger()>3500) 
            Math.floor( 20 / attempt )
        else if (cellno.toInteger()>3000) 
            Math.floor( 30 / attempt )
        else
            Math.floor( 50 / attempt )
}



params.folder = 'by_patient'
ct_cutoff = 20

// Create input channel from params.folder with all patients
input = Channel
    .fromPath(params.folder+'/*.qs', type:"file")
    .map { file ->
        def key = file.name.toString().replaceAll(/.qs/, "").replaceAll(/seurat_/, "")
        return tuple(key,"$file")//.collect { "$it/hg38" }.join(' '))
    }
    .groupTuple()


cellno = Channel.fromPath(params.folder + '/cellno.txt')
    .splitCsv(sep:'\t')

input = input.join(cellno,remainder: true)
 

// --- # Prepare CPDB by defined split_Variable # ---
process cpdb_prep_k43r {
    tag "prepCPDB_$patientkey"
    tag "prepCPDB_$patientkey"
    publishDir "$params.folder/out/$patientkey", mode: 'copy'
    module 'R/4.0.2'
    
    cpus { Math.floor( 20 / task.attempt ) }
    memory { memory_fun(cellno,task.attempt) }
    time '12h'
    
    errorStrategy { task.exitStatus in 137..143 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
    tuple val(patientkey), path("patient.qs"), cellno from input

    output:
    stdout preproc 
    tuple val(patientkey), file("countmatrix_sct.csv"),file("celltype_1.csv"), cellno into cpdb

    script:
    template 'cpdb_prep.R'
}


// --- # CPDB and Heatmap Plot for single object # ---
process cpdb_k43r {
    tag "CPDB_$patientkey"
    tag "CPDB_$patientkey"
    publishDir "$params.folder/out/$patientkey", mode: 'copy'
    module 'conda'
    conda '/home/bq_ssteiger/miniconda3/envs/cellphonedb'

    cpus { core_fun(cellno,task.attempt) }
    memory { memory_fun(cellno,task.attempt) }
    time '12h'
    
    errorStrategy { task.exitStatus in 1..143 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
    tuple val(patientkey), file("countmatrix_sct.csv"),file("celltype_1.csv"), cellno from cpdb
    
    output:
    stdout cpdb_stdout
    path "out/*"

    script:
    template 'cpdb.sh'
}