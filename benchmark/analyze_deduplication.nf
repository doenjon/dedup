#!/usr/bin/env nextflow

nextflow.enable.dsl=2

workflow ANALYZE_DEDUPLICATION {
    take:
        assembly
        illumina_reads
        method_name
        pubDir
        
    main:
        quast_result = QUAST(assembly, pubDir)
        busco_result = BUSCO(assembly, pubDir, method_name)
        kat_result = KAT(assembly, illumina_reads, pubDir)

    emit:
        summary = busco_result.summary
}

process QUAST {
    tag "quast"
    label 'small'
    publishDir { params.results + "/" + pubDir + "/quast" }, mode: "copy"

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
    publishDir { params.results + "/" + pubDir + "/busco" }, mode: "copy"
    errorStrategy = { task.exitStatus == 143 ? 'retry' : 'terminate' }
    maxRetries = 3

    input:
        path(assembly)
        val(pubDir)
        val(method_name)

    output:
        tuple val(method_name), path('busco_output/short_summary.specific.*'), emit: summary

    script:
    """
    set -e  # Exit on error
    busco -i ${assembly} -o busco_output --mode genome --cpu $task.cpus --auto-lineage
    
    # Verify BUSCO completed successfully
    if [ ! -f busco_output/short_summary.specific.* ]; then
        echo "BUSCO failed to produce summary file" >&2
        exit 1
    fi
    """
}

process KAT {
    tag "kat"
    label 'large'
    publishDir { params.results + "/" + pubDir + "/kat" }, mode: "copy"
    conda '/home/groups/ellenyeh/jdoenier/dedup/benchmark/env_kat'

    input:
        path assembly
        tuple path(r1), path(r2)
        val(pubDir)

    output:
        path '*'

    script:
        """
        kat comp -t $task.cpus -o kat_output '${r1} ${r2}' ${assembly}
        """
}

process CONSOLIDATE_SUMMARY {
    tag "consolidate_summary"
    label 'small'
    publishDir { params.results + "/summary" }, mode: "copy"

    input:
        val summaries

    output:
        path "summary.txt"

    script:
    """
    echo -e "Method\tComplete\tSingle-copy\tDuplicated\tFragmented\tMissing" > summary.txt
    
    # Process each summary using simple iteration
    for summary in ${summaries.join(' ')}; do
        method_name=\$(echo "\$summary" | cut -d',' -f1)
        summary_file=\$(echo "\$summary" | cut -d',' -f2)
        
        echo "Processing method: \$method_name with summary file: \$summary_file"
        
        # Extract the summary line that starts with 'C:'
        summary_line=\$(grep "^C:" "\$summary_file")
        
        # Parse the percentages using awk
        stats=\$(echo "\$summary_line" | awk -F'[C:,%\\[\\]]' '{
            complete=\$2
            split(\$3,sd,",")  # Split S:91.0,D:8.0
            single=substr(sd[1],3)  # Remove "S:"
            dupl=substr(sd[2],3)    # Remove "D:"
            fragmented=substr(\$4,3) # Remove "F:"
            missing=substr(\$5,3)    # Remove "M:"
            printf("%s\\t%s\\t%s\\t%s\\t%s", complete, single, dupl, fragmented, missing)
        }')
        
        if [[ -z "\$stats" ]]; then
            echo "Error: Could not parse BUSCO summary line from \$summary_file" >&2
            exit 1
        fi
        
        echo -e "\$method_name\t\$stats" >> summary.txt
    done
    """
}


