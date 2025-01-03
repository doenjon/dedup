#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { FASTP; NANOFILT } from './preprocess_reads.nf'
include { CONCAT_ILLUMINA; CONCAT_NANOPORE; CONCAT_PACBIO } from './preprocess_reads.nf'
include { ASSEMBLE } from './assembly.nf'
include { PURGEDUPS; PURGEHAPLOTIGS; DEDUP } from './deduplicate.nf'
include { ANALYZE_DEDUPLICATION as ANALYZE_DEDUPLICATION_PURGEDUPS } from './analyze_deduplication.nf'
include { ANALYZE_DEDUPLICATION as ANALYZE_DEDUPLICATION_PURGEHAPLOTIGS } from './analyze_deduplication.nf'
include { ANALYZE_DEDUPLICATION as ANALYZE_DEDUPLICATION_DEDUP } from './analyze_deduplication.nf'
include { ANALYZE_DEDUPLICATION as ANALYZE_DEDUPLICATION_ORIGINAL } from './analyze_deduplication.nf'

def printParams(params) {
    params.each { key, value ->
        println "${key}: ${value}"
    }
}

workflow {
    printParams(params)
    println()

    // Prepare read channels
    illumina_reads = Channel.empty()
    ont_reads = Channel.empty()
    pacbio_reads = Channel.empty()

    if (params.illumina_reads) {
        illumina_reads = Channel
            .fromList(params.illumina_reads)
            .map { pair -> tuple(file(pair[0]), file(pair[1])) }
            .set { illumina_pairs }
        
        FASTP(illumina_pairs)
            .collect()
            .map { reads -> tuple(reads.collect{it[0]}, reads.collect{it[1]}) }
            .set { illumina_collected }
            
        illumina_reads = CONCAT_ILLUMINA(illumina_collected)
    }

    if (params.ont_reads) {
        ont_reads = Channel
            .fromPath(params.ont_reads)
            .collect()
            .set { ont_collected }
            
        ont_reads = CONCAT_NANOPORE(ont_collected) | NANOFILT
    }

    if (params.pacbio_reads) {
        pacbio_reads = Channel
            .fromPath(params.pacbio_reads)
            .collect()
            .set { pacbio_collected }
            
        pacbio_reads = CONCAT_PACBIO(pacbio_collected)
    }

    // Check if a genome is provided
    assembly = params.genome 
        ? Channel.fromPath(params.genome) 
        : ASSEMBLE(ont_reads, illumina_reads).polished_assembly

    // Run deduplication algorithms
    purgedups_result = PURGEDUPS(assembly, ont_reads, pacbio_reads)
    purgehaplotigs_result = PURGEHAPLOTIGS(assembly, ont_reads, pacbio_reads)
    dedup_result = DEDUP(assembly, illumina_reads.mix(pacbio_reads).first())

    // Assay performance with BUSCO and KAT
    ANALYZE_DEDUPLICATION_PURGEDUPS(
        purgedups_result.assembly,
        illumina_reads.mix(pacbio_reads).first(),
        "purgedups"
    )
    
    ANALYZE_DEDUPLICATION_PURGEHAPLOTIGS(
        purgehaplotigs_result.assembly,
        illumina_reads.mix(pacbio_reads).first(),
        "purgehaplotigs"
    )
    
    ANALYZE_DEDUPLICATION_DEDUP(
        dedup_result.assembly,
        illumina_reads.mix(pacbio_reads).first(),
        "dedup"
    )
    
    ANALYZE_DEDUPLICATION_ORIGINAL(
        assembly,
        illumina_reads.mix(pacbio_reads).first(),
        "original"
    )
}

