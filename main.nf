#!/usr/bin/env nextflow

include { SAMTOOLS_SORT } from './modules/local/samtools/sort/main'

workflow  {
    unsorted_bam = Channel.fromPath(params.unsorted_bam)
    output_basename = Channel.from(params.output_basename)

    SAMTOOLS_SORT(unsorted_bam, output_basename)
}