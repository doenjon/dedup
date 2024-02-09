
nextflow.enable.dsl=2

process FASTP {
    tag "fastp"
    label 'small'
    publishDir { params.results  + "/fastp" } , mode: "copy"

    input:
        tuple path(r1), path(r2)

    output:
        tuple path("${params.prefix}_R1.fastq.gz"), path("${params.prefix}_R2.fastq.gz"), emit: reads

    script:
        """
        fastp --in1 ${r1} --in2 ${r2} --out1 ${params.prefix}_R1.fastq.gz --out2 ${params.prefix}_R2.fastq.gz 
        """
}

process NANOFILT {
    tag "nanofilt"
    label 'small'
    publishDir { params.results + "/nanofilt" } , mode: "copy"

    input:
        path ont_reads

    output:
        path("${params.prefix}.fastq"), emit: reads

    script:
        """
        cat ${ont_reads} | NanoFilt -l 3000 > ${params.prefix}.fastq
        """
}
