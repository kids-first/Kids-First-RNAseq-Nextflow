#!/usr/bin/env nextflow

include { SAMTOOLS_HEAD_RG_CT } from '../../../modules/local/samtools/head_rg_ct/main'
include { SAMTOOLS_SPLIT } from '../../../modules/local/samtools/split/main'
include { SAMTOOLS_HEAD } from '../../../modules/local/samtools/head/main'
include { ALIGNMENT_PAIREDNESS } from '../../../modules/local/python/pairedness/main'
include { SAMTOOLS_FASTQ } from '../../../modules/local/samtools/fastq/main'
include { CUTADAPT } from '../../../modules/local/cutadapt/main'
include { FASTQ_STRANDEDNESS } from '../../../modules/local/fastq_strandedness/main'


def build_rgs(rg_list, sample){
    // STAR RG becomes the meta
    rg_list = rg_list.map { _meta, rg, file -> [rg.replaceFirst(/^@RG\t/,""), file] }
    rg_list = rg_list.map { rg, file -> [rg.replaceFirst(/\tSM:.+?\t/, "\tSM:${sample.value}\t").trim(), file] }
    return rg_list
}


def parse_strandness_len(strand_file, len_out){
    // Read in strand file data and filter on desired fields
    def strand_vals = strand_file.map {
        f -> f.readLines()
    }.flatten().filter( ~/\w+ReadLength.*|^Data is.*/)
    // Parse generic top length values
    def len_vals = len_out.map { info -> ((info =~ /^\s*\d+\s+(\d+)$/)[0][1]).toInteger() }
    return [len_vals, strand_vals]
}


workflow preprocess_reads {
    take:
    input_alignment_reads // channel: path, Optional unless no input_fastq_reads
    input_fastq_reads // channel: [val(rgs), path(fastq | [fastq])], Optional unless no input_alignment_reads
    line_filter // channel: val(string)
    is_paired_end // boolean: optional
    read_length_median // channel: value(int)
    read_length_stddev // channel: value(int)
    strandedness // channel: val(string)
    max_reads // channel: val(int)
    sample_id // channel: val(string)
    reference // channel: path(fasta)
    annotation_gtf // channel: path(gtf)
    kallisto_idx // channel: path(idx)


    main:
    if (is_paired_end instanceof Boolean){
        is_paired_end = Channel.value(is_paired_end)
    }
    // If reads are from BAM/CRAM, convert to fastq
    if (params.input_alignment_reads){
        // initialize with metadata to track input bams
        align_w_meta = input_alignment_reads.map{ file -> tuple("id": file.baseName, file) }
        SAMTOOLS_HEAD_RG_CT(align_w_meta)
        // Must cast as int as stdout or env will still be a string
        SAMTOOLS_HEAD_RG_CT.out.reads.branch { _meta, rg_ct, _file ->
            single: rg_ct.toInteger() == 1
            multi: rg_ct.toInteger() > 1
        }.set { rg_cts }
        single_rg_bams = rg_cts.single.map{ meta, _rg_ct, file -> [meta, file] }
        multi_rg_bams = rg_cts.multi.map{ meta, _rg_ct, file -> [meta, file] }
        SAMTOOLS_SPLIT(multi_rg_bams, reference)
        // Combine any multi RG outputs that have been split with single RG inputs
        // single_rg_bams looks like [[meta, file], [meta, file]]
        // Each file will have one meta
        // SAMTOOLS_SPLIT.out.bams looks like [[meta, [file1, file2, file3]],[meta, [file1, file2, file3]]]
        // Need to transpose the samtools split outputs so each file gets meta
        align_split_w_meta = single_rg_bams.concat(SAMTOOLS_SPLIT.out.bams.transpose())
        // All below here is left a a practice to the reader
        SAMTOOLS_HEAD(align_split_w_meta, line_filter)
        star_rg_list = build_rgs(SAMTOOLS_HEAD.out, sample_id)
        if (is_paired_end == ""){
            ALIGNMENT_PAIREDNESS(input_alignment_reads, reference, max_reads)
            is_paired_end = ALIGNMENT_PAIREDNESS.out.map{ it ==~ /ReadType:PAIRED/ }
        }
        // use collect to is_paired_end to ensure scatter, so [boolean] is created
        SAMTOOLS_FASTQ(star_rg_list, reference, is_paired_end.collect())
    }
    // Combine fastq inputs with output from alignment conversion block if applicable
    reads = params.input_alignment_reads ? SAMTOOLS_FASTQ.out.concat(input_fastq_reads): input_fastq_reads
    // Set PE channel value flag if not set already
    if (is_paired_end == ""){
        is_paired_end = reads.map{ _reads, f ->
            if (f instanceof List){
                return true
            }
            else{
                return false
            }
        }.unique()
    }
    // only run if a param is missing, otherwise skip
    if (!params.read_length_median || !params.read_length_stddev || !params.strandedness){
        FASTQ_STRANDEDNESS(reads, annotation_gtf, kallisto_idx, max_reads)
        // validate read len and strandedness are consistent
        (top_read_len, check_strandedness) = parse_strandness_len(FASTQ_STRANDEDNESS.out.result, FASTQ_STRANDEDNESS.out.top_read_len)
        // Median and stdev only output for single end data!
        check_strandedness.branch { v->
            stdev: v.startsWith("Stddev")
                return v.split(":")[1].toFloat().toInteger()
            strand: v.startsWith("Data is")
                return v.split(" ")[-1]
            other: true
        }.set { strand_info }
        if (!params.read_length_median){ 
            top_read_len.unique().count().map { n ->
                if (n > 1){
                    error("Inconsistent read lengths")
                }
            }
            read_length_median = top_read_len.unique()
        }
        read_length_stddev = !params.read_length_stddev && is_paired_end.map {it} ? strand_info.stdev.unique() : read_length_stddev
    }
    strand_info.strand.unique().count().map { n ->
        if (n > 1){
            error("Inconsistent strandedness")
        }
    }
    // standardize output of strand_info to match unstranded, rf-stranded, or fr-stranded
    strandedness = params.strandedness ? strandedness : strand_info.strand.unique().map { value -> 
        if (value.toString() == "unstranded"){
            return "default"
        }
        if (value.toString().substring(0, 2) == "RF"){
            return "rf-stranded"
        }
        if (value.toString().substring(0, 2) == "FR"){
            return "fr-stranded"
        }
    }
    // collate user input and calculated values into one list
    added_metadata = is_paired_end.first().concat(read_length_median.first(), strandedness.first(), read_length_stddev.first())
    

    if (params.cutadapt_r1_adapter || params.cutadapt_r2_adapter || params.cutadapt_min_len || params.cutadapt_quality_base || params.cutadapt_quality_cutoff) {
        CUTADAPT(reads)
        reads = CUTADAPT.out.fastq_out
    }
    emit:
    fastq_to_align = reads
    cutadapt_stats = params.cutadapt_r1_adapter || params.cutadapt_r2_adapter || params.cutadapt_min_len || params.cutadapt_quality_base || params.cutadapt_quality_cutoff ? CUTADAPT.out.cutadapt_metrics : Channel.value()
    added_metadata = added_metadata
}