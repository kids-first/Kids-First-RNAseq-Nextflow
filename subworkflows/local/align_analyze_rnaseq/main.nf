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
    genomeDir // channel: path(TAR.GZ)
    readFilesCommand // channel: val(str)
    // Many
    outFileNamePrefix // channel: val(str)
    input_fastq_reads // channel: [val(rgs), path(FASTQ | [FASTQ])]
    genome_tar // channel: path(TAR.GZ)
    reference_fasta // channel: path(FASTA)
    reference_index // channel: path(FAI)
    gtf_anno // channel: path(GTF)
    //RSEM
    RSEMgenome // channel: path(TAR.GZ)
    // kallisto
    kallisto_idx // channel: path(IDX)
    sample_id // channel: va(str)
    // RNAseqc
    RNAseQC_GTF // channel: path(GTF)
    // T1K
    hla_rna_ref_seqs // channel: path(FASTA)
    hla_rna_gene_coords // channel: path(FASTA)
    // From preprocess
    is_paired_end
    strandedness
    read_length_median
    read_length_stddev

    main:
    rgs = input_fastq_reads.map { rg, _fq -> rg}.collect()
    fqs = input_fastq_reads.map { _rg, fq -> fq}.collect()

    STAR_ALIGN(genomeDir, readFilesCommand, rgs, fqs, is_paired_end, outFileNamePrefix)
    STAR_FUSION(genome_tar, STAR_ALIGN.out.chimeric_junctions, outFileNamePrefix)
    SAMTOOLS_SORT(STAR_ALIGN.out.genomic_bam_out, outFileNamePrefix)
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
    genomic_bam_out = STAR_ALIGN.out.genomic_bam_out // channel: path(BAM)
    STAR_sorted_genomic_cram = SAMTOOLS_VIEW.out.cram // channel: tuple path(CRAM), path(CRAI)
    STAR_chimeric_junctions = STAR_FUSION.out.chimeric_junction_compressed // channel: path(TSV.GZ)
    STAR_gene_count = STAR_ALIGN.out.gene_counts // channel: path(GZ)
    STAR_junctions_out = STAR_ALIGN.out.junctions_out // channel: path(GZ)
    STAR_final_log = STAR_ALIGN.out.log_final_out //channel: path(txt)
    STARFusion_results = STAR_FUSION.out.abridged_coding // channel: path(TSV)
    arriba_fusion_results = ARRIBA_FUSION.out.arriba_fusions // channel: path(TSV)
    arriba_fusion_viz = ARRIBA_DRAW.out.arriba_pdf // channel: path(PDF)
    RSEM_isoform = RSEM.out.isoform_out // channel: path(GZ)
    RSEM_gene = RSEM.out.gene_out // channel: path(GZ)
    RNASeQC_Metrics = RNASEQC.out.Metrics // channel: path(TXT)
    RNASeQC_counts = TAR_GZ.out.RNASeQC_counts // channel: path(TAR.GZ)
    kallisto_Abundance = KALLISTO.out.abundance_out // channel: path(TSV.GZ)
    t1k_genotype_tsv = T1K.out.genotype_tsv // channel: path(TSV)
}