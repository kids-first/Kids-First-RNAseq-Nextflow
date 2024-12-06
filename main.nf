#!/usr/bin/env nextflow

include { samtools_split } from './modules/local/samtools/split/main'
include { samtools_head } from './modules/local/samtools/head/main'
include { alignment_pairedness } from './modules/local/python/pairedness/main'
include { align_to_fastq } from './modules/local/samtools/fastq/main'


// def build_rgs(head_file, sample){
//     def rg_list = head_file.readLines()
//     rg_list = rg_list.replaceAll(/(^@RG\t)(.*)(\tSM:.+?\t)(.*)/) { rgtag, pre_sm, sm, post_sm -> "$pre_sm\t$sample\t$post_sm" }.trim()
//     return rg_list
// }

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


workflow rnaseq_wf {
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
    if (input_alignment_reads != "")
        samtools_split(input_alignment_reads, reference, threads)
        samtools_head(samtools_split.out.flatten(), line_filter)
        // def star_rg = build_rgs(samtools_head.out, sample_id)
        if ( is_paired_end == ""){
            alignment_pairedness(input_alignment_reads, reference, max_reads, output_filename, threads)
            is_paired_end = eval_align_pairedness(alignment_pairedness.out)
        }
        align_to_fastq(samtools_split.out.flatten(), reference, sample_id, threads, is_paired_end)

    
    emit:
    fq1 = align_to_fastq.out.fq1 ?: input_fastq_reads
    fq2 = align_to_fastq.out.fq2 ?: input_fastq_mates
    rg_file = samtools_head.out
}

workflow {
    input_alignment_reads = Channel.fromPath(params.input_alignment_reads)
    input_fastq_reads = Channel.fromPath(params.input_fastq_reads)
    input_fastq_mates = Channel.fromPath(params.input_fastq_mates)
    is_paired_end = Channel.value(params.is_paired_end)
    max_reads = Channel.value(params.max_reads)
    output_filename = Channel.from(params.output_filename)
    line_filter = Channel.from(params.line_filter)
    sample_id = Channel.value(params.sample_id)
    threads = Channel.from(params.threads)
    reference = Channel.fromPath(params.reference)
    rnaseq_wf(input_alignment_reads, input_fastq_reads, input_fastq_mates, line_filter, is_paired_end, max_reads, output_filename, sample_id, threads, reference)

}