#!/usr/bin/env nextflow

include { RMATS } from '../../../modules/local/rmats/tool/main'
include { AWK_JC_FILTER } from '../../../modules/local/rmats/filter/main'

workflow rmats_subworkflow {
    take:
    gtf_annotation // channel: path(GTF)
    sample_1_bams // channel: path(BAM)
    read_length // channel: val(int)
    read_type // channel: val(str)
    strandedness // channel: val(str)

    main:
    RMATS(gtf_annotation, sample_1_bams, read_length, read_type, strandedness)
    jc_files = RMATS.out.alternative_3_prime_splice_sites_jc.concat(
        RMATS.out.alternative_5_prime_splice_sites_jc,
        RMATS.out.mutually_exclusive_exons_jc,
        RMATS.out.retained_introns_jc,
        RMATS.out.skipped_exons_jc
        )
    AWK_JC_FILTER(jc_files)
    
    emit:
    rmats_filtered_jc = AWK_JC_FILTER.out.output // channel: [path(TXT)]

}