#!/usr/bin/env nextflow

include { SAMTOOLS_HEAD_RG_CT } from './modules/local/samtools/head_rg_ct/main'
include { SAMTOOLS_SPLIT } from './modules/local/samtools/split/main'
include { SAMTOOLS_HEAD } from './modules/local/samtools/head/main'
include { ALIGNMENT_PAIREDNESS } from './modules/local/python/pairedness/main'
include { SAMTOOLS_FASTQ } from './modules/local/samtools/fastq/main'
include { CUTADAPT } from './modules/local/cutadapt/main'
include { SAMTOOLS_SORT } from './modules/local/samtools/sort/main'
include { FASTQ_STRANDEDNESS } from './modules/local/fastq_strandedness/main'
include { STAR_ALIGN } from './modules/local/star/align/main'
include { STAR_FUSION } from './modules/local/star/fusion/main'
include { RSEM } from './modules/local/rsem'
include { KALLISTO } from './modules/local/kallisto/main'
include { ARRIBA_FUSION } from './modules/local/arriba/fusion/main'
include { ARRIBA_DRAW } from './modules/local/arriba/draw/main'
include { FORMAT_ARRIBA } from './modules/local/annofuse/format_arriba/main'
include { ANNOTATE_ARRIBA } from './modules/local/annofuse/annotate_arriba/main'
include { ANNOFUSE } from './modules/local/annofuse/tool/main'
include { RMATS } from './modules/local/rmats/tool/main'
include { AWK_JC_FILTER } from './modules/local/rmats/filter/main'
include { RNASEQC } from './modules/local/rnaseqc/main'
include { T1K } from './modules/local/t1k/main'
include { SAMTOOLS_VIEW } from './modules/local/samtools/view/main'
include { TAR_GZ } from './modules/local/ubuntu/tar_gz/main'


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
    input_alignment_reads
    input_fastq_reads
    line_filter
    is_paired_end
    read_length_median
    read_length_stddev
    strandedness
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
    in_fq_formatted = params.input_fastq_reads ? input_fastq_reads.map{ meta, f -> [meta, f instanceof List ? f.collect{ file(it,checkIfExists: true) } : file(f, checkIfExists: true)] } : Channel.empty()
    reads = params.input_alignment_reads ? SAMTOOLS_FASTQ.out.concat(in_fq_formatted): in_fq_formatted
    // Set PE channel value flag if not set already
    if (is_paired_end instanceof Boolean){
        is_paired_end = Channel.value(is_paired_end)
    }
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

workflow annofuse_subworkflow {
    take:
    arriba_output_file
    sample_name
    fusion_annotator_tar
    rsem_expr_file
    star_fusion_output_file
    output_basename

    main:
    FORMAT_ARRIBA(arriba_output_file, output_basename)
    ANNOTATE_ARRIBA(FORMAT_ARRIBA.out.formatted_fusion_tsv, fusion_annotator_tar, output_basename)
    ANNOFUSE(ANNOTATE_ARRIBA.out.annotated_tsv, star_fusion_output_file, rsem_expr_file, sample_name, output_basename)

    emit:
    annofuse_filtered_fusions_tsv = ANNOFUSE.out.filtered_fusions_tsv
}

workflow rmats_subworkflow {
    take:
    gtf_annotation
    sample_1_bams
    read_length
    read_type
    strandedness
    output_basename

    main:
    RMATS(gtf_annotation, sample_1_bams, read_length, read_type, strandedness, output_basename)
    jc_files = RMATS.out.alternative_3_prime_splice_sites_jc.concat(
        RMATS.out.alternative_5_prime_splice_sites_jc,
        RMATS.out.mutually_exclusive_exons_jc,
        RMATS.out.retained_introns_jc,
        RMATS.out.skipped_exons_jc
        )
    AWK_JC_FILTER(jc_files)
    
    emit:
    rmats_filtered_jc = AWK_JC_FILTER.out.output

}
workflow align_analyze_rnaseq {
    take:
    //STAR
    genomeDir
    readFilesCommand
    // Many
    outFileNamePrefix
    input_fastq_reads
    genome_tar
    samtools_threads
    reference_fasta
    reference_index
    gtf_anno
    //RSEM
    RSEMgenome
    // kallisto
    kallisto_idx
    sample_id
    // RNAseqc
    RNAseQC_GTF
    // T1K
    hla_rna_ref_seqs
    hla_rna_gene_coords
    // From preprocess
    added_metadata

    main:
    rgs = input_fastq_reads.map { rg, _fq -> rg}.collect()
    fqs = input_fastq_reads.map { _rg, fq -> fq}.collect()

     // Assign added metadata to named variables for clarity
    (is_paired_end, read_length_median, strandedness) = [added_metadata.first(), added_metadata.take(2).last(), added_metadata.take(3).last()]
    read_length_stddev = is_paired_end ? added_metadata.last() : ""

    STAR_ALIGN(genomeDir, readFilesCommand, rgs, fqs, is_paired_end, outFileNamePrefix)
    STAR_FUSION(genome_tar, STAR_ALIGN.out.chimeric_junctions, outFileNamePrefix)
    SAMTOOLS_SORT(STAR_ALIGN.out.genomic_bam_out, outFileNamePrefix, samtools_threads)
    // Create a value conversion dict as many tools use strand as a param but call it different things

    wf_strand_info = [ "ARRIBA_FUSION": strandedness.map{it.startsWith("rf") ? "reverse" : (it.startsWith("fr") ? "yes" : "auto")},
        "RSEM": strandedness.map{it.startsWith("rf") ? "reverse" : (it.startsWith("fr") ? "forward" : "none")},
        "KALLISTO": strandedness.map{it.startsWith("rf") ? "rf-stranded" : (it.startsWith("fr") ? "fr-stranded" : "")},
        "RNASEQC": strandedness.map{it.startsWith("rf") ? "rf" : (it.startsWith("fr") ? "fr" : "")}
        ]
    KALLISTO(kallisto_idx, wf_strand_info["KALLISTO"], fqs, sample_id, read_length_stddev, read_length_median, is_paired_end)

    sorted_bam_bai = SAMTOOLS_SORT.out.sorted_bam.combine(SAMTOOLS_SORT.out.sorted_bai)
    ARRIBA_FUSION(sorted_bam_bai, reference_fasta, gtf_anno, outFileNamePrefix, wf_strand_info.ARRIBA_FUSION)
    ARRIBA_DRAW(sorted_bam_bai, ARRIBA_FUSION.out.arriba_fusions, gtf_anno)
    RSEM(RSEMgenome, STAR_ALIGN.out.transcriptome_bam_out, outFileNamePrefix, is_paired_end, wf_strand_info.RSEM)
    RNASEQC(RNAseQC_GTF, sorted_bam_bai, wf_strand_info["RNASEQC"], is_paired_end)
    TAR_GZ(RNASEQC.out.Gene_TPM, RNASEQC.out.Gene_count, RNASEQC.out.Exon_count, outFileNamePrefix)
    T1K(sorted_bam_bai, hla_rna_ref_seqs, hla_rna_gene_coords, outFileNamePrefix)

    reference_fai = reference_fasta.combine(reference_index)
    SAMTOOLS_VIEW(reference_fai, sorted_bam_bai)


    emit:
    genomic_bam_out = STAR_ALIGN.out.genomic_bam_out
    STAR_sorted_genomic_cram = SAMTOOLS_VIEW.out.cram
    STAR_chimeric_junctions = STAR_FUSION.out.chimeric_junction_compressed
    STAR_gene_count = STAR_ALIGN.out.gene_counts
    STAR_junctions_out = STAR_ALIGN.out.junctions_out
    STAR_final_log = STAR_ALIGN.out.log_final_out
    STARFusion_results = STAR_FUSION.out.abridged_coding
    arriba_fusion_results = ARRIBA_FUSION.out.arriba_fusions
    arriba_fusion_viz = ARRIBA_DRAW.out.arriba_pdf
    RSEM_isoform = RSEM.out.isoform_out
    RSEM_gene = RSEM.out.gene_out
    RNASeQC_Metrics = RNASEQC.out.Metrics
    RNASeQC_counts = TAR_GZ.out.RNASeQC_counts
    kallisto_Abundance = KALLISTO.out.abundance_out
    t1k_genotype_tsv = T1K.out.genotype_tsv
}

workflow {
    main:
    input_alignment_reads = params.input_alignment_reads ? Channel.fromPath(params.input_alignment_reads) : Channel.value([])
    input_fastq_reads = params.input_fastq_reads ? Channel.fromList(params.input_fastq_reads) : Channel.value([])
    is_paired_end = params.is_paired_end != "" ? params.is_paired_end : ""
    read_length_median = params.read_length_median ? Channel.value(params.read_length_median) : Channel.value([])
    read_length_stddev = params.read_length_stddev ? Channel.value(params.read_length_stddev) : Channel.value([])
    strandedness = params.strandedness ? Channel.value(params.strandedness) : Channel.value([])
    max_reads = Channel.value(params.max_reads)
    line_filter = Channel.value(params.line_filter)
    sample_id = Channel.value(params.sample_id)
    threads = Channel.value(params.threads)
    reference = Channel.fromPath(params.reference).first()
    reference_index = Channel.fromPath(params.reference_index).first()
    output_basename = Channel.value(params.output_basename)
    gtf_anno = Channel.fromPath(params.gtf_anno).first()
    kallisto_idx = Channel.fromPath(params.kallisto_idx).first()
    // STAR
    genomeDir = Channel.fromPath(params.genomeDir)
    readFilesCommand = Channel.value(params.readFilesCommand)
    // STAR Fusion
    FusionGenome = Channel.fromPath(params.FusionGenome)
    // arriba
    fusion_annotator_tar = Channel.fromPath(params.fusion_annotator_tar)
    assembly = Channel.value(params.assembly)

    samtools_threads = Channel.value(params.samtools_threads)
    // RSEM
    RSEM_genome = Channel.fromPath(params.RSEM_genome)
    // RNASeQC
    RNAseQC_GTF = Channel.fromPath(params.RNAseQC_GTF)
    // T1K
    hla_rna_ref_seqs = Channel.fromPath(params.hla_rna_ref_seqs)
    hla_rna_gene_coords = Channel.fromPath(params.hla_rna_gene_coords)

    preprocess_reads(input_alignment_reads, input_fastq_reads, line_filter, is_paired_end, read_length_median, read_length_stddev, strandedness, max_reads, sample_id, threads, reference, gtf_anno, kallisto_idx)

    added_metadata = preprocess_reads.out.added_metadata
    (is_paired_end, read_length_median, strandedness) = [added_metadata.first(), added_metadata.take(2).last(), added_metadata.take(3).last()]
    rmats_strand = strandedness.map{it.startsWith("rf") ? "fr-firststrand" : (it.startsWith("fr") ? "fr-secondstrand" : "fr-unstranded")}

    align_analyze_rnaseq(genomeDir, readFilesCommand, output_basename, preprocess_reads.out.fastq_to_align, FusionGenome, samtools_threads, reference, reference_index, gtf_anno, RSEM_genome, kallisto_idx, sample_id, RNAseQC_GTF, hla_rna_ref_seqs, hla_rna_gene_coords, added_metadata)
    annofuse_subworkflow(align_analyze_rnaseq.out.arriba_fusion_results, sample_id, fusion_annotator_tar, align_analyze_rnaseq.out.RSEM_gene, align_analyze_rnaseq.out.STARFusion_results, output_basename)
    rmats_subworkflow(gtf_anno, align_analyze_rnaseq.out.genomic_bam_out, read_length_median, is_paired_end ? "paired" : "single", rmats_strand, output_basename)

}
