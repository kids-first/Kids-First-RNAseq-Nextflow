#!/usr/bin/env nextflow

include { SAMTOOLS_SORT } from './modules/local/samtools/sort/main'
include { samtools_split } from './modules/local/samtools/split/main'


workflow rnaseq_wf {
    take:
    input_alignment_reads
    output_basename
    threads 
    reference

    main:
    if (input_alignment_reads != "")
        samtools_split(input_alignment_reads, reference, threads)
    SAMTOOLS_SORT(unsorted_bam, output_basename, threads)
    
    emit:
    samtools_split.out ?: SAMTOOLS_SORT.out.sorted_bam
}

workflow {
    input_alignment_reads = Channel.fromPath(params.input_alignment_reads)
    output_basename = Channel.from(params.output_basename)
    threads = Channel.from(params.threads)
    reference = Channel.fromPath(params.reference)
    rnaseq_wf(input_alignment_reads, output_basename, threads, reference)

}