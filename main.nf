#!/usr/bin/env nextflow

include { preprocess_reads } from './subworkflows/local/preprocess_reads/main'
include { align_analyze_rnaseq } from './subworkflows/local/align_analyze_rnaseq/main'
include { annofuse_subworkflow } from './subworkflows/local/annofuse_subworkflow/main'
include { rmats_subworkflow } from './subworkflows/local/rmats_subworkflow/main'

workflow {
    main:
    input_alignment_reads = params.input_alignment_reads ? Channel.fromPath(params.input_alignment_reads) : Channel.value([]) // channel: [path(bam/cram)]
    input_fastq_reads = params.input_fastq_reads ? Channel.fromList(params.input_fastq_reads).map{ meta, f -> [meta, f instanceof List ? f.collect{ file(it,checkIfExists: true) } : file(f, checkIfExists: true)] } : Channel.empty() // channel: [val(rgs), path(fastq | [fastq])], Optional unless no input_alignment_reads
    is_paired_end = params.is_paired_end != "" ? params.is_paired_end : "" // channel: val(boolean) optional
    read_length_median = params.read_length_median ? Channel.value(params.read_length_median) : Channel.value([]) // channel: val(int), optional
    read_length_stddev = params.read_length_stddev ? Channel.value(params.read_length_stddev) : Channel.value([]) // channel: val(int), optional
    strandedness = params.strandedness ? Channel.value(params.strandedness) : Channel.value([]) // channel: val(int), optional
    max_reads = Channel.value(params.max_reads) // channel: val(int)
    line_filter = Channel.value(params.line_filter) // channel: val(str)
    sample_id = Channel.value(params.sample_id) // channel: val(str)
    reference = Channel.fromPath(params.reference).first() // channel: path(FASTA)
    reference_index = Channel.fromPath(params.reference_index).first() // channel: path(FAI)
    output_basename = Channel.value(params.output_basename) // channel: val(str)
    gtf_anno = Channel.fromPath(params.gtf_anno).first() // channel: path(GTF)
    kallisto_idx = Channel.fromPath(params.kallisto_idx).first() // channel: path(IDX)
    // STAR
    genomeDir = Channel.fromPath(params.genomeDir) //channel: path(TAR.GZ)
    readFilesCommand = Channel.value(params.readFilesCommand) // channel: va(str)
    // STAR Fusion
    FusionGenome = Channel.fromPath(params.FusionGenome) // channel: path(TAR.GZ)
    // arriba
    fusion_annotator_tar = Channel.fromPath(params.fusion_annotator_tar) // channel: path(TAR/GZ)
    // RSEM
    RSEM_genome = Channel.fromPath(params.RSEM_genome) // channel: path(TAR.GZ)
    // RNASeQC
    RNAseQC_GTF = Channel.fromPath(params.RNAseQC_GTF) // channel: path(GTF)
    // T1K
    hla_rna_ref_seqs = Channel.fromPath(params.hla_rna_ref_seqs) // channel: path(FASTA)
    hla_rna_gene_coords = Channel.fromPath(params.hla_rna_gene_coords) // channel: path(FASTA)

    preprocess_reads(input_alignment_reads, input_fastq_reads, line_filter, is_paired_end, read_length_median, read_length_stddev, strandedness, max_reads, sample_id, reference, gtf_anno, kallisto_idx)

    added_metadata = preprocess_reads.out.added_metadata
    (is_paired_end, read_length_median, strandedness) = [added_metadata.first(), added_metadata.take(2).last(), added_metadata.take(3).last()]
    rmats_strand = strandedness.map{it.startsWith("rf") ? "fr-firststrand" : (it.startsWith("fr") ? "fr-secondstrand" : "fr-unstranded")}

    align_analyze_rnaseq(genomeDir, readFilesCommand, output_basename, preprocess_reads.out.fastq_to_align, FusionGenome, reference, reference_index, gtf_anno, RSEM_genome, kallisto_idx, sample_id, RNAseQC_GTF, hla_rna_ref_seqs, hla_rna_gene_coords, added_metadata)
    annofuse_subworkflow(align_analyze_rnaseq.out.arriba_fusion_results, sample_id, fusion_annotator_tar, align_analyze_rnaseq.out.RSEM_gene, align_analyze_rnaseq.out.STARFusion_results, output_basename)
    rmats_subworkflow(gtf_anno, align_analyze_rnaseq.out.genomic_bam_out, read_length_median, is_paired_end ? "paired" : "single", rmats_strand, output_basename)

}
