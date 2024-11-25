#!/usr/bin/env nextflow

include { SAMTOOLS_SORT } from './modules/local/samtools/sort/main'
include { SAMTOOLS_SPLIT } from './modules/local/samtools/split/main'
workflow  {
    unsorted_bam = Channel.fromPath(params.unsorted_bam)
    output_basename = Channel.from(params.output_basename)
    threads = Channel.from(params.threads)
    reference = Channel.fromPath(params.reference)

    SAMTOOLS_SORT(unsorted_bam, output_basename, threads)
    SAMTOOLS_SPLIT(SAMTOOLS_SORT.out.sorted_bam, reference, threads)
}