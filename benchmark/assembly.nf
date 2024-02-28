
nextflow.enable.dsl=2


workflow ASSEMBLE {
  take:
    long_reads
    illumina_reads  

  main:
    assembly = FLYE(long_reads).asm
    polished_assembly = NEXTPOLISH(assembly, illumina_reads, long_reads).asm

  emit:
    polished_assembly
}

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


process FLYE {
  /* 
  A process to perform genome assembly using flye
  */ 

  publishDir { params.results  + "/Flye" } , mode: "copy"
  label "very_large"

  input:
	path(reads)   // nanopore reads in fastq format

  output:
	path("assembly.fasta"), emit: asm
	path("*")

  script:
	"""
	flye --scaffold --nano-raw ${reads} --threads ${task.cpus} --out-dir ./
	"""
}

process NEXTPOLISH {
  /* 
  A process to perform nextPolish polishing using nanopore reads
  */ 

  publishDir { params.results +  "/" + "nextPolish" } , mode: "copy"
  label "very_large"
  conda '/home/groups/ellenyeh/jdoenier/dedup/benchmark/env_nextpolish'

  input:
    path(assembly)
    tuple path(sr_R1), path(sr_R2)
    path(ont_reads)

  output:
    path("run_dir/*.fasta"), emit: asm
    path("*")

  script:
    """

    cat > run.cfg <<- EOM
    [General]
    job_type = local
    job_prefix = nextPolish
    task = best
    rewrite = yes
    rerun = 3
    parallel_jobs = 4
    multithread_jobs = 7
    genome = ${assembly}
    genome_size = auto
    workdir = ./run_dir
    polish_options = -p 7

    [sgs_option]
    sgs_fofn = ./sgs.fofn
    sgs_options = -max_depth 100 -bwa

    [lgs_option]
    lgs_fofn = ./lgs.fofn
    lgs_options = -min_read_len 1k -max_depth 100
    lgs_minimap2_options = -x map-ont

    EOM
    
    ls ${sr_R1} ${sr_R2} > sgs.fofn
    ls ${ont_reads} > lgs.fofn

    # set +e # don't quit on error...
    /home/groups/ellenyeh/jdoenier/bin/NextPolish/nextPolish run.cfg
 
    """
}