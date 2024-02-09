#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process PURGEDUPS {
    tag "purgedups"
    label 'large'
    publishDir { params.results + "/" + "purgedups" } , mode: "copy"

    input:
        path assembly
        path ont_reads

    output:
        path "purged.fa", emit: assembly
        path "*"

    script:
    """
    minimap2 -t $task.cpus -x map-ont ${assembly} ${ont_reads} | gzip -c - > mapped_reads.paf.gz
    pbcstat mapped_reads.paf.gz
    calcuts PB.stat > cutoffs 2>calcults.log
    
    split_fa ${assembly} > ${assembly}.split
    minimap2 -x map-ont -DP ${assembly}.split ${assembly}.split | gzip -c - > ${assembly}.split.self.paf.gz
    purge_dups -2 -T cutoffs -c PB.base.cov ${assembly}.split.self.paf.gz > dups.bed 2> purge_dups.log
    get_seqs -e dups.bed ${assembly}
    """
}

process PURGEHAPLOTIGS {
    tag "purgehaplotigs"
    label 'large'
    publishDir { params.results + "/" + "purgehaplotigs" } , mode: "copy"

    input:
        path assembly
        path ont_reads

    output:
        path "*.curated.fasta", emit: assembly
        path "*"

    script:
    """
    minimap2 -t $task.cpus -ax map-ont ${assembly} ${ont_reads} | samtools view -b | samtools sort -m 1G -@ $task.cpus -o mapped_reads.sorted.bam
    purge_haplotigs hist -b mapped_reads.sorted.bam -g ${assembly} -t $task.cpus
    purge_haplotigs cov -i mapped_reads.sorted.bam.gencov -l $params.l -m $params.m -h $params.h -o coverage_stats.csv
    purge_haplotigs purge -g ${assembly} -c coverage_stats.csv -o ${assembly}.curated -t $task.cpus
    """
}

process FASTPURGE {
    tag "fastpurge"
    label 'large'
    publishDir { params.results + "/" + "fastpurge" } , mode: "copy"

    input:
        path assembly
        tuple path(r1), path(r2)


    output:
        path "deduplicated_contigs.fasta", emit: assembly
        path "*"

    script:
    """
    python3 /home/groups/ellenyeh/jdoenier/dedup/dedup/dedup.py --reads ${r1} --assembly ${assembly}
    """
}