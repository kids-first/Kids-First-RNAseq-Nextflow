#!/usr/bin/env nextflow

include { SAMTOOLS_SPLIT } from './modules/local/samtools/split/main'
include { SAMTOOLS_HEAD } from './modules/local/samtools/head/main'
include { ALIGNMENT_PAIREDNESS } from './modules/local/python/pairedness/main'
include { SAMTOOLS_FASTQ } from './modules/local/samtools/fastq/main'
include { CUTADAPT } from './modules/local/cutadapt/main'


def build_rgs(rg_list, sample){
    rg_list = rg_list.map { r -> r.replaceFirst(/^@RG\t/,"") }
    rg_list = rg_list.map { r -> r.replaceFirst(/\tSM:.+?\t/, "\tSM:${sample.value}\t") }
    return rg_list
}

def eval_align_pairedness(result_file){
    def result = file(result_file).readLines()
    println result[0]
    if (result[0] =~ /PAIRED/)
        return true
    if (result[0] =~ /SINGLE/)
        return false
    if (result[0] =~ /MIXED/)
        error "Could not accurately determine pairedness of input align file, please manually specify"
}


workflow preprocess_reads {
    take:
    input_alignment_reads
    input_fastq_reads
    input_fastq_mates
    line_filter
    is_paired_end
    max_reads
    output_filename
    sample_id
    threads 
    reference

    main:

    // If reads are from BAM/CRAM, convert to fastq
    input_alignment_reads.view()
    if (params.input_alignment_reads){
        SAMTOOLS_SPLIT(input_alignment_reads, reference, threads)
        SAMTOOLS_HEAD(SAMTOOLS_SPLIT.out.flatten(), line_filter)
        def raw_rg_flist = SAMTOOLS_HEAD.out.collect { v -> v.readLines() }
        star_rg_list = build_rgs(raw_rg_flist.flatten(), sample_id)
        if ( is_paired_end == ""){
            ALIGNMENT_PAIREDNESS(input_alignment_reads, reference, max_reads, output_filename, threads)
            is_paired_end = eval_align_pairedness(ALIGNMENT_PAIREDNESS.out)
        }
        SAMTOOLS_FASTQ(SAMTOOLS_SPLIT.out.flatten(), reference, threads, is_paired_end)
    }
    if (params.cutadapt_r1_adapter || params.cutadapt_r2_adapter || params.cutadapt_min_len || params.cutadapt_quality_base || params.cutadapt_quality_cutoff) {
        CUTADAPT(SAMTOOLS_FASTQ.out.fq1 ?: input_fastq_reads, SAMTOOLS_FASTQ.out.fq2 ?: input_fastq_mates)
    }
    emit:
    fastq_to_align = (params.cutadapt_r1_adapter || params.cutadapt_r2_adapter || params.cutadapt_min_len || params.cutadapt_quality_base || params.cutadapt_quality_cutoff) ? CUTADAPT.out.fastq_out :
        (SAMTOOLS_FASTQ.out.fq1.merge(SAMTOOLS_FASTQ.out.fq2) ?: [input_fastq_reads, input_fastq_mates])
    rg_file = SAMTOOLS_HEAD.out
    star_rg_list
}

workflow {
    main:
    input_alignment_reads = params.input_alignment_reads ? Channel.fromPath(params.input_alignment_reads) : Channel.value([])
    input_fastq_reads = params.input_fastq_reads ? Channel.fromPath(params.input_fastq_reads) : Channel.value([])
    input_fastq_mates = params.input_fastq_mates ? Channel.fromPath(params.input_fastq_mates) : Channel.value([])
    is_paired_end = Channel.value(params.is_paired_end)
    max_reads = Channel.value(params.max_reads)
    output_filename = Channel.value(params.output_filename)
    line_filter = Channel.value(params.line_filter)
    sample_id = Channel.value(params.sample_id)
    threads = Channel.value(params.threads)
    reference = Channel.fromPath(params.reference).first()
    preprocess_reads(input_alignment_reads, input_fastq_reads, input_fastq_mates, line_filter, is_paired_end, max_reads, output_filename, sample_id, threads, reference)
    preprocess_reads.out.fastq_to_align.view()
    preprocess_reads.out.star_rg_list.view()

    publish:
    preprocess_reads.out.fastq_to_align >> 'reads_to_align'
    preprocess_reads.out.star_rg_list >> 'read_groups'
}