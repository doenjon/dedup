#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { FASTP; NANOFILT } from './preprocess_reads.nf'
include { PURGEDUPS; PURGEHAPLOTIGS; FASTPURGE } from './deduplicate.nf'
include { ANALYZE_DEDUPLICATION as ANALYZE_DEDUPLICATION1 } from './analyze_deduplication.nf'
include { ANALYZE_DEDUPLICATION as ANALYZE_DEDUPLICATION2 } from './analyze_deduplication.nf'
include { ANALYZE_DEDUPLICATION as ANALYZE_DEDUPLICATION3 } from './analyze_deduplication.nf'
include { ANALYZE_DEDUPLICATION as ANALYZE_DEDUPLICATION4 } from './analyze_deduplication.nf'

def printParams(params) {
    params.each { key, value ->
        println "${key}: ${value}"
    }
}

workflow {
    printParams(params)

    // Preprocess reads
    Channel.fromList(params.illumina_reads).view().println()

    illumina_reads = FASTP(Channel.fromList(params.illumina_reads))
    ont_reads = NANOFILT(params.ont_reads)

    // Run deduplication algorithms
    purgedups_result = PURGEDUPS(params.assembly, ont_reads.reads)
    purgehaplotigs_result = PURGEHAPLOTIGS(params.assembly, ont_reads.reads)
    fastpurge_result = FASTPURGE(params.assembly, illumina_reads.reads)

    // Assay performance with BUSCO and KAT
    ANALYZE_DEDUPLICATION1(purgedups_result.assembly, illumina_reads, "purgedups")
    ANALYZE_DEDUPLICATION2(purgehaplotigs_result.assembly, illumina_reads, "purgehaplotigs")
    ANALYZE_DEDUPLICATION3(fastpurge_result.assembly, illumina_reads, "fastpurge")
    ANALYZE_DEDUPLICATION4(params.assembly, illumina_reads, "original")

    // Print out params
    printParams(params)
}

