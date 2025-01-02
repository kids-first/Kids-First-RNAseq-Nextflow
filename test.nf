#!/usr/bin/env nextflow
include { CUTADAPT } from './modules/local/cutadapt/main'

workflow  {
    reads = Channel.fromFilePairs(params.reads)
    reads.view()
    CUTADAPT(reads)
}