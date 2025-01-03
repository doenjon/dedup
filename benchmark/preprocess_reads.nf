nextflow.enable.dsl=2

process FASTP {
    tag "fastp"
    label 'small'
    publishDir { params.results  + "/fastp" } , mode: "copy"

    input:
        tuple path(r1), path(r2)

    output:
        tuple path("cleaned_R1.fastq"), path("cleaned_R2.fastq"), emit: reads

    script:
        """
        fastp --in1 ${r1} --in2 ${r2} --out1 cleaned_R1.fastq --out2 cleaned_R2.fastq
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

process CONCAT_PACBIO {
    tag "concat_pacbio"
    label 'small'
    publishDir { params.results + "/concat_pacbio" }, mode: "copy"

    input:
        path reads, stageAs: 'reads_*'

    output:
        path "concatenated_pacbio.fastq", emit: reads

    script:
        """
        # Handle both compressed and uncompressed files
        for file in ${reads}; do
            if [[ \$file == *.gz ]]; then
                gunzip -c \$file
            else
                cat \$file
            fi
        done > concatenated_pacbio.fastq
        """
}

process CONCAT_NANOPORE {
    tag "concat_nanopore"
    label 'small'
    publishDir { params.results + "/concat_nanopore" }, mode: "copy"

    input:
        path reads, stageAs: 'reads_*'

    output:
        path "concatenated_nanopore.fastq", emit: reads

    script:
        """
        # Handle both compressed and uncompressed files
        for file in ${reads}; do
            if [[ \$file == *.gz ]]; then
                gunzip -c \$file
            else
                cat \$file
            fi
        done > concatenated_nanopore.fastq
        """
}

process CONCAT_ILLUMINA {
    tag "concat_illumina"
    label 'small'
    publishDir { params.results + "/concat_illumina" }, mode: "copy"

    input:
        tuple path(r1s), path(r2s)  // Receives lists of R1s and R2s separately

    output:
        tuple path("concatenated_R1.fastq"), path("concatenated_R2.fastq"), emit: reads

    script:
        """
        # Concatenate in specified order
        cat ${r1s.join(' ')} > concatenated_R1.fastq
        cat ${r2s.join(' ')} > concatenated_R2.fastq
        """
}
