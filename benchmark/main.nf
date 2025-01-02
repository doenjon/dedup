#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { FASTP; NANOFILT } from './preprocess_reads.nf'
include { ASSEMBLE } from './assembly.nf'
include { PURGEDUPS; PURGEHAPLOTIGS; FASTPURGE } from './deduplicate.nf'
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

    // Preprocess reads
    illumina_reads = FASTP(Channel.fromList(params.illumina_reads))
    long_reads = NANOFILT(params.long_reads)

    // Perform Genome Assembly
    assembly = ASSEMBLE(long_reads.reads, illumina_reads.reads).polished_assembly

    // Run deduplication algorithms
    purgedups_result = PURGEDUPS(assembly, long_reads.reads)
    purgehaplotigs_result = PURGEHAPLOTIGS(assembly, long_reads.reads)
    dedup_result = DEDUP(assembly, illumina_reads.reads)

    // Assay performance with BUSCO and KAT
    ANALYZE_DEDUPLICATION_PURGEDUPS(purgedups_result.assembly, illumina_reads, "purgedups")
    ANALYZE_DEDUPLICATION_PURGEHAPLOTIGS(purgehaplotigs_result.assembly, illumina_reads, "purgehaplotigs")
    ANALYZE_DEDUPLICATION_DEDUP(dedup_result.assembly, illumina_reads, "fastpurge")
    ANALYZE_DEDUPLICATION_ORIGINAL(assembly, illumina_reads, "original")

}

