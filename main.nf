#!/usr/bin/env nextflow

include { SAMTOOLS_HEAD_RG_CT } from './modules/local/samtools/head_rg_ct/main'
include { SAMTOOLS_SPLIT } from './modules/local/samtools/split/main'
include { SAMTOOLS_HEAD } from './modules/local/samtools/head/main'
include { ALIGNMENT_PAIREDNESS } from './modules/local/python/pairedness/main'
include { SAMTOOLS_FASTQ } from './modules/local/samtools/fastq/main'
include { CUTADAPT } from './modules/local/cutadapt/main'



def build_rgs(rg_list, sample){
    // STAR RG becomes the meta
    rg_list = rg_list.map { _meta, rg, file -> [rg.replaceFirst(/^@RG\t/,""), file] }
    rg_list = rg_list.map { rg, file -> [rg.replaceFirst(/\tSM:.+?\t/, "\tSM:${sample.value}\t").trim(), file] }
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
    input_rg_strs
    line_filter
    is_paired_end
    max_reads
    output_filename
    sample_id
    threads 
    reference

    main:

    // If reads are from BAM/CRAM, convert to fastq
    if (params.input_alignment_reads){

        align_w_meta = input_alignment_reads.map{ file -> tuple("id": file.baseName, file) }
        SAMTOOLS_HEAD_RG_CT(align_w_meta)
        // Must cast as int as stdout or env will still be a string
        SAMTOOLS_HEAD_RG_CT.out.reads.branch { _meta, rg_ct, _file ->
            single: rg_ct.toInteger() == 1
            multi: rg_ct.toInteger() > 1
        }.set { rg_cts }
        single_rg_bams = rg_cts.single.map{ meta, _rg_ct, file -> [meta, file] }
        multi_rg_bams = rg_cts.multi.map{ meta, _rg_ct, file -> [meta, file] }
        SAMTOOLS_SPLIT(multi_rg_bams, reference, threads)
        // Combine any multi RG outputs that have been split with single RG inputs
        // single_rg_bams looks like [[meta, file], [meta, file]]
        // Each file will have one meta
        // SAMTOOLS_SPLIT.out.bams looks like [[meta, [file1, file2, file3]],[meta, [file1, file2, file3]]]
        // Need to transpose the samtools split outputs so each file gets meta
        align_split_w_meta = single_rg_bams.concat(SAMTOOLS_SPLIT.out.bams.transpose())

        // All below here is left a a practice to the reader
        SAMTOOLS_HEAD(align_split_w_meta, line_filter)
        star_rg_list = build_rgs(SAMTOOLS_HEAD.out, sample_id)
        if ( is_paired_end == ""){
            ALIGNMENT_PAIREDNESS(input_alignment_reads, reference, max_reads, output_filename, threads)
            is_paired_end = eval_align_pairedness(ALIGNMENT_PAIREDNESS.out)
        }
        SAMTOOLS_FASTQ(star_rg_list, reference, threads, is_paired_end)
    }

    if (params.cutadapt_r1_adapter || params.cutadapt_r2_adapter || params.cutadapt_min_len || params.cutadapt_quality_base || params.cutadapt_quality_cutoff) {
        in_fq_formatted = params.input_fastq_reads ? Channel.fromList(input_fastq_reads).map{ meta, f -> [meta, f.collect{ file(it,checkIfExists: true) }] }.view() : Channel.empty()
        reads = params.input_alignment_reads ? SAMTOOLS_FASTQ.out.concat(in_fq_formatted): in_fq_formatted
        reads.view()
        CUTADAPT(reads)
    }
    // emit:
    // fastq_to_align = (params.cutadapt_r1_adapter || params.cutadapt_r2_adapter || params.cutadapt_min_len || params.cutadapt_quality_base || params.cutadapt_quality_cutoff) ? CUTADAPT.out.fastq_out :
    //     (SAMTOOLS_FASTQ.out ?: input_fastq_reads)
    // star_rg_list = params.input_alignment_reads ? star_rg_list : input_rg_strs
}

workflow {
    main:
    input_alignment_reads = params.input_alignment_reads ? Channel.fromPath(params.input_alignment_reads) : Channel.value([])
    input_fastq_reads = params.input_fastq_reads ? params.input_fastq_reads : Channel.value([])
    input_rg_strs = params.input_rg_strs ? Channel.value(params.input_rg_strs) : Channel.value([])
    is_paired_end = Channel.value(params.is_paired_end)
    max_reads = Channel.value(params.max_reads)
    output_filename = Channel.value(params.output_filename)
    line_filter = Channel.value(params.line_filter)
    sample_id = Channel.value(params.sample_id)
    threads = Channel.value(params.threads)
    reference = Channel.fromPath(params.reference).first()
    preprocess_reads(input_alignment_reads, input_fastq_reads, input_rg_strs, line_filter, is_paired_end, max_reads, output_filename, sample_id, threads, reference)
    // preprocess_reads.out.fastq_to_align.view()

    // publish:
    // preprocess_reads.out.fastq_to_align >> 'reads_to_align'
    // preprocess_reads.out.star_rg_list >> 'read_groups'
}
