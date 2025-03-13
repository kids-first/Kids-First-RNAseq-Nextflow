#!/usr/bin/env nextflow

include { FORMAT_ARRIBA } from '../../../modules/local/annofuse/format_arriba/main'
include { ANNOTATE_ARRIBA } from '../../../modules/local/annofuse/annotate_arriba/main'
include { ANNOFUSE } from '../../../modules/local/annofuse/tool/main'

workflow annofuse_subworkflow {
    take:
    arriba_output_file // channel: path(TSC)
    sample_id // channel: val(str)
    fusion_annotator_tar // channel: path(TAR.GZ)
    rsem_expr_file // channel: path(GZ)
    star_fusion_output_file // channel: path(TSV)
    output_basename // channel: val(str)

    main:
    FORMAT_ARRIBA(arriba_output_file, output_basename)
    ANNOTATE_ARRIBA(FORMAT_ARRIBA.out.formatted_fusion_tsv, fusion_annotator_tar, output_basename)
    ANNOFUSE(ANNOTATE_ARRIBA.out.annotated_tsv, star_fusion_output_file, rsem_expr_file, sample_id, output_basename)

    emit:
    annofuse_filtered_fusions_tsv = ANNOFUSE.out.filtered_fusions_tsv // channel: path(TSV)
}