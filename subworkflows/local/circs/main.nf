#!/usr/bin/env nextflow

include { SAMTOOLS_FASTQ } from './modules/samtools/fastq/main'
include { STAR_ALIGN as STAR_ALIGN_DCC_PAIR } from './modules/star/align/main'
include { STAR_ALIGN as STAR_ALIGN_DCC_R1 } from './modules/star/align/main'
include { STAR_ALIGN as STAR_ALIGN_DCC_R2 } from './modules/star/align/main'
include { DCC_MAIN } from './modules/dcc/main/main'
include { DCC_OUTREADER } from './modules/dcc/outreader/main'
include { BEDTOOLS_WINDOW } form './modules/bedtools/window/main'

workflow {
    input_aligned_reads = Channel.fromPath(params.aligned_reads).map{ [["id": it.baseName, "is_paired_end": params.is_paired_end.toBoolean()], it] }
    cram_reference = Channel.fromPath(params.cram_reference).first()
    star_genome = Channel.fromPath(params.star_genome).first()
    reference = Channel.fromPath(params.reference).first()
    refseq_gtf = Channel.fromPath(param.refseq_gtf).first()
    refseq_bed = Channel.fromPath(param.refseq_bed).first()
    
    SAMTOOLS_FASTQ(input_aligned_reads, cram_reference)

    paired_fq_channel = SAMTOOLS_FASTQ.out.fastq.map { meta, files -> [meta + ["id": meta.id + "_bothreads"], files] }
    r1_fq_channel = SAMTOOLS_FASTQ.out.fastq.map { meta, files -> [["id": meta.id + "_read1", "is_paired_end": false], files[0]] }
    r2_fq_channel = SAMTOOLS_FASTQ.out.fastq.map { meta, files -> [["id": meta.id + "_read2", "is_paired_end": false], files[1]] }

    paired_fq_channel.view()
    r1_fq_channel.view()
    r2_fq_channel.view()

    if (params.run_dcc) {
        STAR_ALIGN_DCC_PAIR(paired_fq_channel, star_genome)
        STAR_ALIGN_DCC_R1(r1_fq_channel, star_genome)
        STAR_ALIGN_DCC_R2(r2_fq_channel, star_genome)

        main_ch = STAR_ALIGN_DCC_PAIR.out.chimeric_junctions.merge(STAR_ALIGN_DCC_R1.out.chimeric_junctions, STAR_ALIGN_DCC_R2.out.chimeric_junctions) 
        main_ch.view()

        //DCC_MAIN()
        // tuple val(meta), path(paired_junctions), path(read1_junctions), path(read2_junctions)
        // path(ref_fasta)
        // path(refseq_bed)
        //BEDTOOLS_WINDOW()
        //    tuple val(meta), path(intervals1), path(intervals2)
        //DCC_OUTREADER()
        //    tuple val(meta), path(annotated_counts), path(cleaned_coordinates)
    }
}
