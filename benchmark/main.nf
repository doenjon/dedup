#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { FASTP; NANOFILT } from './preprocess_reads.nf'
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
    illumina_reads = params.illumina_reads ? Channel
        .fromFilePairs(params.illumina_reads)
        | FASTP
        | collect()
        | map { reads -> [reads.collect{it[0]}, reads.collect{it[1]}] }  // Group R1s and R2s
        | CONCAT_ILLUMINA : null

    ont_reads = params.ont_reads ? Channel
        .fromList(params.ont_reads)
        .collect()
        | CONCAT_NANOPORE
        | NANOFILT : null

    pacbio_reads = params.pacbio_reads ? Channel
        .fromList(params.pacbio_reads)
        .collect()
        | CONCAT_PACBIO : null

    // Check if a genome is provided
    if (params.genome) {
        assembly = Channel.value(params.genome)
    } else {
        // Perform Genome Assembly
        assembly = ASSEMBLE(ont_reads?.reads, illumina_reads?.reads).polished_assembly
    }

    // Run deduplication algorithms
    purgedups_result = PURGEDUPS(assembly, ont_reads?.reads, pacbio_reads)
    purgehaplotigs_result = PURGEHAPLOTIGS(assembly, ont_reads?.reads, pacbio_reads)
    dedup_result = DEDUP(assembly, illumina_reads?.reads ?: pacbio_reads)

    // Assay performance with BUSCO and KAT
    ANALYZE_DEDUPLICATION_PURGEDUPS(purgedups_result.assembly, illumina_reads ?: pacbio_reads, "purgedups")
    ANALYZE_DEDUPLICATION_PURGEHAPLOTIGS(purgehaplotigs_result.assembly, illumina_reads ?: pacbio_reads, "purgehaplotigs")
    ANALYZE_DEDUPLICATION_DEDUP(dedup_result.assembly, illumina_reads ?: pacbio_reads, "dedup")
    ANALYZE_DEDUPLICATION_ORIGINAL(assembly, illumina_reads ?: pacbio_reads, "original")
}

