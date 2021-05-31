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


// Create input channel from folder input with all tumor samples merged by Patient 
input = Channel
    .fromPath("input/*", type:"file")
    .map { file ->
        def key = file.name.toString().replaceAll(/.qs/, "")
        return tuple(key,"$file")
    }
    .groupTuple()


cellno = Channel.fromPath('tumorcellno.txt')
    .splitCsv(sep:'\t')

input = input.join(cellno,remainder: true)
 

// --- # InferCNV on tumor samples # ---
process infercnv_k43r {
    tag "infercnv_$patientkey"
    tag "infercnv_$patientkey"
    publishDir "infercnv_out/$patientkey", mode: 'copy'
    label 'rsing'
    
    cpus { Math.floor( 15 / task.attempt ) }
    memory { memory_fun(cellno,task.attempt) }
    time '72h'
    
    errorStrategy { task.exitStatus in 137..143 ? 'retry' : 'terminate' }
    maxRetries 3


    input:
    tuple val(patientkey), path("tumor.qs"), cellno from input
    path reference from "/media/ag-rippe2/NGS_Simon/HIPO_K43R_nextflow_2/reference/HCA_PCs.qs"

    output:
    stdout preproc 
    tuple val(patientkey), file("infercnv.qs"), cellno into infercnv
    path "InferCNV_out"


    script:
    template '01_InferCNV.R'
}
