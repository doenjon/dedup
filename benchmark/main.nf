#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { FASTP; NANOFILT } from './preprocess_reads.nf'
include { ASSEMBLE } from './assembly.nf'
include { PURGEDUPS; PURGEHAPLOTIGS; DEDUP } from './deduplicate.nf'
include { ANALYZE_DEDUPLICATION; CONSOLIDATE_SUMMARY } from './analyze_deduplication.nf'

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

    // Create a channel for the original assembly
    original_assembly = assembly.map{ it -> [it, "original"] }

    // Run deduplication algorithms and combine with original assembly
    purgedups_result = PURGEDUPS(assembly, long_reads.reads).assembly.map{ it -> [it, "purgedups"] }
    purgehap_result = PURGEHAPLOTIGS(assembly, long_reads.reads).assembly.map{ it -> [it, "purgehaplotigs"] }
    dedup_result = DEDUP(assembly, illumina_reads.reads).assembly.map{ it -> [it, "dedup"] }

    // Combine all results into one channel
    dedup_results = original_assembly
        .mix(purgedups_result)
        .mix(purgehap_result)
        .mix(dedup_result)

    // Analyze each assembly using the workflow
    busco_summaries = dedup_results.map { asm, method_name ->
        ANALYZE_DEDUPLICATION(asm, illumina_reads, method_name, method_name).summary
    }.collect()

    // Consolidate summary
    CONSOLIDATE_SUMMARY(busco_summaries)
}

