#!/usr/bin/env nextflow

nextflow.enable.dsl=2

workflow ANALYZE_DEDUPLICATION {

    take:
        assembly
        illumina_reads
        pubDir
        
   
    main:
        quast_result = QUAST(assembly, pubDir)
        busco_result = BUSCO(assembly, pubDir)
        kat_result = KAT(assembly, illumina_reads, pubDir)

    emit:
        quast_result
        busco_result
        kat_result
}

process QUAST {
    tag "quast"
    label 'small'
    publishDir { params.results + pubDir + "/" + "quast" } , mode: "copy"


    input:
        path(assembly)
        val(pubDir)

    output:
        path 'quast_output'

    script:
        """
        quast.py ${assembly} -o quast_output -t $task.cpus
        """
}

process BUSCO {
    tag "busco"
    label 'large'
    publishDir { params.results + pubDir + "/" + "busco" } , mode: "copy"


    input:
        path(assembly)
        val(pubDir)

    output:
        path 'busco_output'

    script:
    """
    busco -i ${assembly} -o busco_output --mode genome --cpu $task.cpus --auto-lineage
    """
}

process KAT {
    tag "kat"
    label 'large'
    publishDir { params.results + pubDir + "/" + "kat" } , mode: "copy"
    conda '/home/groups/ellenyeh/jdoenier/dedup/benchmark/env_kat'

    input:
        path assembly
        tuple path(r1), path(r2) optional true
        path pacbio_reads optional true
        val(pubDir)

    output:
        path '*'

    script:
    """
    if [ -n "${r1}" ] && [ -n "${r2}" ]; then
        kat comp -t $task.cpus -o kat_output '${r1} ${r2}' ${assembly}
    elif [ -n "${pacbio_reads}" ]; then
        kat comp -t $task.cpus -o kat_output ${pacbio_reads} ${assembly}
    else
        echo "Error: No reads provided for KAT" >&2
        exit 1
    fi
    """
}


