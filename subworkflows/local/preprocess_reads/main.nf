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


def qc_pe_values(flag){
    def e_msg = ""
    flag.unique().count().map { n ->
        if (n > 1){
            e_msg += "Inconsistent pairedness\n"
        }
    }
    flag.map{
        if(!it instanceof Boolean){
            e_msg += "Align pairedness could not determine value.\n"
        }
    }
    if (e_msg == ""){
        return flag.unique()
    }
    else{
        if ( params.is_paired_end instanceof Boolean ){
            return params.is_paired_end
        }
    }
    error(e_msg + "Conflicting pairedness detected during pipeline. Please review inputs and/or set is_paired_end flag manually")

}


def qc_strand_values(strand){
    strand.unique().count().map { n ->
        if (n > 1){
            if (params.strandedness != ""){
                return params.strandedness
            }
            else{
                error("Inconsistent strandedness. Revisit inputs and/or set params.strandedness to override")
            }
        }
    }
    return strand.unique()
}


workflow preprocess_reads {
    take:
    input_alignment_reads // channel: path, Optional unless no input_fastq_reads
    input_fastq_reads // channel: [val(rgs), path(FASTQ | [FASTQ])], Optional unless no input_alignment_reads
    line_filter // channel: val(string)
    read_length_median // channel: value(int)
    read_length_stddev // channel: value(int)
    max_reads // channel: val(int)
    sample_id // channel: val(string)
    reference // channel: path(FASTA)
    annotation_gtf // channel: path(GTF)
    kallisto_idx // channel: path(IDX)


    main:
    // If reads are from BAM/CRAM, convert to fastq
    // initialize with metadata to track input bams
    align_w_meta = params.input_alignment_reads ? input_alignment_reads.map{ file -> tuple("id": file.baseName, file) } : Channel.empty()
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
    ALIGNMENT_PAIREDNESS(align_split_w_meta, reference, max_reads)
    is_paired_end = ALIGNMENT_PAIREDNESS.out.result
    // use collect to is_paired_end to ensure scatter, so [boolean] is created
    SAMTOOLS_FASTQ(star_rg_list, reference, is_paired_end.collect())
    // Combine fastq inputs with output from alignment conversion block if applicable
    reads = SAMTOOLS_FASTQ.out.concat(input_fastq_reads)
    // Collate PE flag with FASTQ inputs as a check
    is_paired_end = is_paired_end.concat(
        reads.map{ _rgs, f ->
        if (f instanceof List){
            return true
        }
        else{
            return false
        }
    }
    )
    // Final check of is_paired_end - all should match, if not, need params override
    is_paired_end = qc_pe_values(is_paired_end)
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
            def test = v.split(" ")[-1].toString().substring(0, 2)
                return test == "unstranded" ? "default" : (test == "RF" ? "rf-stranded" : "fr-stranded")
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

    // standardize output of strand_info to match unstranded, rf-stranded, or fr-stranded
    strandedness = params.strandedness ?: strand_info.strand
    strandedness = qc_pe_values(strandedness)
    // collate user input and calculated values into one list

    if (params.cutadapt_r1_adapter || params.cutadapt_r2_adapter || params.cutadapt_min_len || params.cutadapt_quality_base || params.cutadapt_quality_cutoff) {
        CUTADAPT(reads)
        reads = CUTADAPT.out.fastq_out
    }
    emit:
    fastq_to_align = reads // channel: [val(rgs), path(FASTQ | [FASTQ])]
    cutadapt_stats = params.cutadapt_r1_adapter || params.cutadapt_r2_adapter || params.cutadapt_min_len || params.cutadapt_quality_base || params.cutadapt_quality_cutoff ? CUTADAPT.out.cutadapt_metrics : Channel.value() // channel: path(TXT) Optional
    is_paired_end = is_paired_end
    read_length_median = read_length_median
    read_length_stddev = read_length_stddev
    strandedness = strandedness

}