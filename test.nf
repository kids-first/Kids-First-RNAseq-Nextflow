#!/usr/bin/env nextflow
include { CUTADAPT } from './modules/local/cutadapt/main'

workflow  {
    reads = Channel.fromPath(params.reads)
    mates = params.mates ? Channel.fromPath(params.mates) : Channel.value([])

    CUTADAPT(reads, mates)
    CUTADAPT.out.fastq_out.view()
}