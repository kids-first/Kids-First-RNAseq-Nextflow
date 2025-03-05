#!/usr/bin/env nextflow


include { STAR_ALIGN } from '../../../modules/local/star/align/main'
include { STAR_FUSION } from '../../../modules/local/star/fusion/main'
include { RSEM } from '../../../modules/local/rsem'
include { KALLISTO } from '../../../modules/local/kallisto/main'
include { ARRIBA_FUSION } from '../../../modules/local/arriba/fusion/main'
include { ARRIBA_DRAW } from '../../../modules/local/arriba/draw/main'
include { SAMTOOLS_SORT } from '../../../modules/local/samtools/sort/main'
include { RNASEQC } from '../../../modules/local/rnaseqc/main'
include { T1K } from '../../../modules/local/t1k/main'
include { SAMTOOLS_VIEW } from '../../../modules/local/samtools/view/main'
include { TAR_GZ } from '../../../modules/local/ubuntu/tar_gz/main'


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