#!/usr/bin/env nextflow

include { SAMTOOLS_HEAD_RG_CT } from './modules/local/samtools/head_rg_ct/main'
include { SAMTOOLS_SPLIT } from './modules/local/samtools/split/main'
include { SAMTOOLS_HEAD } from './modules/local/samtools/head/main'
include { ALIGNMENT_PAIREDNESS } from './modules/local/python/pairedness/main'
include { SAMTOOLS_FASTQ } from './modules/local/samtools/fastq/main'
include { CUTADAPT } from './modules/local/cutadapt/main'
include { STAR_ALIGN } from './modules/local/star/align/main'
include { SAMTOOLS_SORT } from './modules/local/samtools/sort/main'
include { FASTQ_STRANDEDNESS } from './modules/local/fastq_strandedness/main'


def build_rgs(rg_list, sample){
    // STAR RG becomes the meta
    rg_list = rg_list.map { _meta, rg, file -> [rg.replaceFirst(/^@RG\t/,""), file] }
    rg_list = rg_list.map { rg, file -> [rg.replaceFirst(/\tSM:.+?\t/, "\tSM:${sample.value}\t").trim(), file] }
    return rg_list
}


def validate_strandness_len(strand_file, len_out){
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
    input_alignment_reads
    input_fastq_reads
    line_filter
    is_paired_end
    max_reads
    sample_id
    threads 
    reference
    annotation_gtf
    kallisto_idx


    main:

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
        if (is_paired_end == ""){
            ALIGNMENT_PAIREDNESS(input_alignment_reads, reference, max_reads, threads)
            is_paired_end = ALIGNMENT_PAIREDNESS.out.map{ it ==~ /ReadType:PAIRED/ }
        }
        else {
            is_paired_end = Channel.value(is_paired_end)
        }
        // use collect to is_paired_end to ensure scatter, so [boolean] is created
        SAMTOOLS_FASTQ(star_rg_list, reference, threads, is_paired_end.collect())
    }
    // reformat fastq inputs to match output from alignment conversion block
    in_fq_formatted = params.input_fastq_reads ? Channel.fromList(input_fastq_reads).map{ meta, f -> [meta, f instanceof List ? f.collect{ file(it,checkIfExists: true) } : file(f, checkIfExists: true)] } : Channel.empty()
    reads = params.input_alignment_reads ? SAMTOOLS_FASTQ.out.concat(in_fq_formatted): in_fq_formatted
    FASTQ_STRANDEDNESS(reads, annotation_gtf, kallisto_idx, max_reads)
    // validate read len and strandedness are consistent
    check = validate_strandness_len(FASTQ_STRANDEDNESS.out.result, FASTQ_STRANDEDNESS.out.top_read_len)
    check[1].branch { v->
        median: v.startsWith("Median")
            return v.split(":")[1].toFloat().toInteger()
        stdev: v.startsWith("Stddev")
            return v.split(":")[1].toFloat().toInteger()
        strand: v.startsWith("Data is")
            return v.split(" ")[-1]
        other: true
    }.set { strand_info }

    check[0].unique().count().map { n ->
        if (n > 1){
            error("Inconsistent read lengths")
        }
    }
    strand_info.median.view()
    strand_info.stdev.view()
    strand_info.strand.unique().count().map { n ->
        if (n > 1){
            error("Inconsistent strandedness")
        }
    }

    if (params.cutadapt_r1_adapter || params.cutadapt_r2_adapter || params.cutadapt_min_len || params.cutadapt_quality_base || params.cutadapt_quality_cutoff) {
        CUTADAPT(reads)
        reads = CUTADAPT.out.fastq_out
    }
    emit:
    fastq_to_align = reads
    is_paired_end
    median_read_len = strand_info.median.first()
    stdev_read_len = strand_info.stdev.first()
    strand = strand_info.strand.first()
}

workflow align_analyze_rnaseq {
    take:
    genomeDir
    readFilesCommand
    outFileNamePrefix
    readFilesManifest
    samtools_threads


    main:
    STAR_ALIGN(genomeDir, readFilesCommand, readFilesManifest, outFileNamePrefix)
    SAMTOOLS_SORT(STAR_ALIGN.out.genomic_bam_out, outFileNamePrefix, samtools_threads)

    emit:
    genomic_bam_out = STAR_ALIGN.out.genomic_bam_out

}

workflow {
    main:
    input_alignment_reads = params.input_alignment_reads ? Channel.fromPath(params.input_alignment_reads) : Channel.value([])
    input_fastq_reads = params.input_fastq_reads ? params.input_fastq_reads : Channel.value([])
    is_paired_end = params.is_paired_end != "" ? params.is_paired_end : ""
    max_reads = Channel.value(params.max_reads)
    line_filter = Channel.value(params.line_filter)
    sample_id = Channel.value(params.sample_id)
    threads = Channel.value(params.threads)
    reference = Channel.fromPath(params.reference).first()
    // output_basename = Channel.value(params.output_basename)
    gtf_anno = Channel.fromPath(params.gtf_anno).first()
    kallisto_idx = Channel.fromPath(params.kallisto_idx).first()

    // genomeDir = Channel.fromPath(params.genomeDir)
    // readFilesCommand = Channel.value(params.readFilesCommand)

    // samtools_threads = Channel.value(params.samtools_threads)

    preprocess_reads(input_alignment_reads, input_fastq_reads, line_filter, is_paired_end, max_reads, sample_id, threads, reference, gtf_anno, kallisto_idx)
    // preprocess_reads.out.fastq_to_align.view()
    // Create STAR reads manifest from fastq object for multi-read group support
    // star_reads_manifest = preprocess_reads.out.fastq_to_align.map{
    //     rg, fastq -> [fastq instanceof List ? fastq.join('\t'): fastq + '\t-', rg].join('\t')
    // }.collectFile( name: 'star_reads_manifest.txt', newLine: true)
    // align_analyze_rnaseq(genomeDir, readFilesCommand, output_basename, star_reads_manifest, samtools_threads)

}
