#!/usr/bin/env nextflow

include { SAMTOOLS_FASTQ } from './modules/samtools/fastq/main'
include { STAR_ALIGN as STAR_ALIGN_DCC_PAIR } from './modules/star/align/main'
include { STAR_ALIGN as STAR_ALIGN_DCC_R1 } from './modules/star/align/main'
include { STAR_ALIGN as STAR_ALIGN_DCC_R2 } from './modules/star/align/main'
include { STAR_ALIGN as STAR_ALIGN_CX } from './modules/star/align/main'
include { DCC_MAIN } from './modules/dcc/main/main'
include { DCC_OUTREADER } from './modules/dcc/outreader/main'
include { BEDTOOLS_WINDOW } from './modules/bedtools/window/main'
include { BEDTOOLS_WINDOW as BEDTOOLS_WINDOW_FC } from './modules/bedtools/window/main'
include { CIRCEXPLORER_STARPARSE } from './modules/circexplorer/starparse/main'
include { CIRCEXPLORER_MAIN } from './modules/circexplorer/main/main'
include { CIRCEXPLORER_OUTREADER } from './modules/circexplorer/outreader/main'
include { BOWTIE2_ALIGN as BOWTIE2_ALIGN_PRIMARY } from './modules/bowtie2/align/main'
include { SAMTOOLS_VIEW } from './modules/samtools/view/main'
include { SAMTOOLS_SPLITFASTA } from './modules/samtools/splitfasta/main'
include { FINDCIRC_ANCHORS } from './modules/findcirc/anchors/main'
include { FINDCIRC_MAIN } from './modules/findcirc/main/main'
include { FINDCIRC_FILTER } from './modules/findcirc/filter/main'
include { FINDCIRC_OUTREADER } from './modules/findcirc/outreader/main'

workflow {
    input_aligned_reads = Channel.fromPath(params.aligned_reads).map{ [["id": it.baseName, "is_paired_end": params.is_paired_end.toBoolean()], it] }
    cram_reference = Channel.fromPath(params.cram_reference).first()
    star_genome = Channel.fromPath(params.star_genome).first()
    reference = Channel.fromPath(params.reference).first()
    reference_index = Channel.fromPath(params.reference_index).first()
    refseq_gtf = Channel.fromPath(params.refseq_gtf).first()
    refseq_bed = Channel.fromPath(params.refseq_bed).first()
    refseq_annot = Channel.fromPath(params.refseq_annot).first()
    bowtie_index = Channel.fromPath(params.bowtie_index).first()
    
    SAMTOOLS_FASTQ(input_aligned_reads, cram_reference)

    if (params.run_dcc) {
        paired_fq_channel = SAMTOOLS_FASTQ.out.fastq.map { meta, files -> [meta + ["id": meta.id + "_bothreads"], files] }
        r1_fq_channel = SAMTOOLS_FASTQ.out.fastq.map { meta, files -> [["id": meta.id + "_read1", "is_paired_end": false], files[0]] }
        r2_fq_channel = SAMTOOLS_FASTQ.out.fastq.map { meta, files -> [["id": meta.id + "_read2", "is_paired_end": false], files[1]] }

        STAR_ALIGN_DCC_PAIR(paired_fq_channel, star_genome)
        STAR_ALIGN_DCC_R1(r1_fq_channel, star_genome)
        STAR_ALIGN_DCC_R2(r2_fq_channel, star_genome)

        a = STAR_ALIGN_DCC_PAIR.out.chimeric_junctions.map { _, file -> [["id": "test"], file] }
        b = STAR_ALIGN_DCC_R1.out.chimeric_junctions.map { _, file -> [["id": "test"], file] }
        c = STAR_ALIGN_DCC_R2.out.chimeric_junctions.map { _, file -> [["id": "test"], file] }
        juntions = a.join(b).join(c)
        x = STAR_ALIGN_DCC_PAIR.out.junctions_out.map { _, file -> [["id": "test"], file] }
        y = STAR_ALIGN_DCC_R1.out.junctions_out.map { _, file -> [["id": "test"], file] }
        z = STAR_ALIGN_DCC_R2.out.junctions_out.map { _, file -> [["id": "test"], file] }
        tabs = x.mix(y,z).groupTuple()
        dcc_main_ch = juntions.join(tabs)

        DCC_MAIN(dcc_main_ch, reference, refseq_gtf)
        BEDTOOLS_WINDOW(DCC_MAIN.out.counts, refseq_bed)
        outreader_ch = BEDTOOLS_WINDOW.out.windows.join(DCC_MAIN.out.coordinates)
        DCC_OUTREADER(outreader_ch)
    }

    if (params.run_cx) {
        println "Running CX"
        cx_fq_channel = SAMTOOLS_FASTQ.out.fastq.map { meta, files -> [meta + ["id": meta.id + "_cx"], files] }

        STAR_ALIGN_CX(cx_fq_channel, star_genome)
        chimeric_junc_channel = STAR_ALIGN_CX.out.chimeric_junctions.map { _, file -> [["id": "test"], file] }

        CIRCEXPLORER_STARPARSE(chimeric_junc_channel)

        CIRCEXPLORER_MAIN(CIRCEXPLORER_STARPARSE.out.fusion_junctions, reference, refseq_annot)
        CIRCEXPLORER_OUTREADER(CIRCEXPLORER_MAIN.out.circs)
    }
    if (params.run_fc) {
        println "Running FC"
        fc_fq_channel = SAMTOOLS_FASTQ.out.fastq.map { meta, files -> [meta + ["id": meta.id + "_fc"], files] }
        BOWTIE2_ALIGN_PRIMARY(fc_fq_channel, bowtie_index, cram_reference, Channel.value(false), Channel.value(true)) 
        primary_alignment_channel = BOWTIE2_ALIGN_PRIMARY.out.bam.join(BOWTIE2_ALIGN_PRIMARY.out.csi).map{ meta, file, index -> [meta + ["id": meta.id + "_unmapped"], file, index] }
        primary_alignment_channel.view()
        SAMTOOLS_VIEW(primary_alignment_channel, Channel.value([[],[]]), Channel.value([]), Channel.value([])) 
        FINDCIRC_ANCHORS(SAMTOOLS_VIEW.out.bam)
        SAMTOOLS_SPLITFASTA(reference.merge(reference_index))
        FINDCIRC_MAIN(FINDCIRC_ANCHORS.out.anchors.map{ meta, file -> [meta + ["single_end": true], file]}, bowtie_index, SAMTOOLS_SPLITFASTA.out.split_fasta, Channel.value(false))
        FINDCIRC_MAIN.out.bed_ci.view()
        FINDCIRC_FILTER(FINDCIRC_MAIN.out.bed_ci)
        BEDTOOLS_WINDOW_FC(FINDCIRC_FILTER.out.candidates, refseq_bed)
        FINDCIRC_OUTREADER(BEDTOOLS_WINDOW_FC.out.windows)
    }
}
