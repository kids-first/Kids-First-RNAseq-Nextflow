#!/usr/bin/env nextflow

include { SAMTOOLS_FASTQ } from './modules/samtools/fastq/main'
include { STAR_ALIGN as STAR_ALIGN_DCC_PAIR } from './modules/star/align/main'
include { STAR_ALIGN as STAR_ALIGN_DCC_R1 } from './modules/star/align/main'
include { STAR_ALIGN as STAR_ALIGN_DCC_R2 } from './modules/star/align/main'
include { STAR_ALIGN as STAR_ALIGN_CX } from './modules/star/align/main'
include { DCC_MAIN } from './modules/dcc/main/main'
include { DCC_OUTREADER } from './modules/dcc/outreader/main'
include { BEDTOOLS_WINDOW as BEDTOOLS_WINDOW_DCC } from './modules/bedtools/window/main'
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
include { CIRCSNAKE_MATRIXMAKER as CIRCSNAKE_MATRIXMAKER_DCC } from './modules/circsnake/matrixmaker/main'
include { CIRCSNAKE_MATRIXMAKER as CIRCSNAKE_MATRIXMAKER_FC } from './modules/circsnake/matrixmaker/main'
include { CIRCSNAKE_MATRIXMAKER as CIRCSNAKE_MATRIXMAKER_CX } from './modules/circsnake/matrixmaker/main'
include { CIRCSNAKE_MATRIXTWO as CIRCSNAKE_MATRIXTWO_DCC } from './modules/circsnake/matrixtwo/main'
include { CIRCSNAKE_MATRIXTWO as CIRCSNAKE_MATRIXTWO_FC } from './modules/circsnake/matrixtwo/main'
include { CIRCSNAKE_MATRIXTWO as CIRCSNAKE_MATRIXTWO_CX } from './modules/circsnake/matrixtwo/main'
include { CIRCSNAKE_VOTE } from './modules/circsnake/vote/main'
include { CIRCSNAKE_READCOUNTS } from './modules/circsnake/readcounts/main'
include { CIRCSNAKE_NORM } from './modules/circsnake/norm/main'

workflow {
    input_aligned_reads = channel.fromPath(params.aligned_reads.class == String ? params.aligned_reads.split(',') as List : params.aligned_reads).map{ [["id": it.baseName, "is_paired_end": params.is_paired_end.toBoolean()], it] }
    cram_reference = channel.fromPath(params.cram_reference).first()
    star_genome = channel.fromPath(params.star_genome).first()
    reference = channel.fromPath(params.reference).first()
    reference_index = channel.fromPath(params.reference_index).first()
    refseq_gtf = channel.fromPath(params.refseq_gtf).first()
    refseq_bed = channel.fromPath(params.refseq_bed).first()
    refseq_annot = channel.fromPath(params.refseq_annot).first()
    bowtie_index = channel.fromPath(params.bowtie_index).first()
    mm1_refseq = channel.fromPath(params.mm1_refseq).first()
    mm1_circ = channel.fromPath(params.mm1_circ).first()
    micrornas = channel.fromPath(params.micrornas).first()
    coding_circs = channel.fromPath(params.coding_circs).first()
    hallmarks = channel.fromPath(params.hallmarks).first()
    ensembl_gene_descriptions = channel.fromPath(params.ensembl_gene_descriptions).first()

    SAMTOOLS_FASTQ(input_aligned_reads, cram_reference)

    CIRCSNAKE_READCOUNTS(SAMTOOLS_FASTQ.out.fastq.map { meta, files -> "${meta.id}\t${files[0].countFastq()}" }.collect())

    if (params.run_dcc) {
        log.info "Running DCC"

        paired_fq_channel = SAMTOOLS_FASTQ.out.fastq
        r1_fq_channel = SAMTOOLS_FASTQ.out.fastq.map { meta, files -> [meta + ["is_paired_end": false], files[0]] }
        r2_fq_channel = SAMTOOLS_FASTQ.out.fastq.map { meta, files -> [meta + ["is_paired_end": false], files[1]] }

        STAR_ALIGN_DCC_PAIR(paired_fq_channel, star_genome)
        STAR_ALIGN_DCC_R1(r1_fq_channel, star_genome)
        STAR_ALIGN_DCC_R2(r2_fq_channel, star_genome)

        x = STAR_ALIGN_DCC_PAIR.out.junctions_out.map { meta, file -> [meta.subMap("id"), file] }
        y = STAR_ALIGN_DCC_R1.out.junctions_out.map { meta, file -> [meta.subMap("id"), file] }
        z = STAR_ALIGN_DCC_R2.out.junctions_out.map { meta, file -> [meta.subMap("id"), file] }
        tabs = x.mix(y,z).groupTuple()
        a = STAR_ALIGN_DCC_PAIR.out.chimeric_junctions.map { meta, file -> [meta.subMap("id"), file] }
        b = STAR_ALIGN_DCC_R1.out.chimeric_junctions.map { meta, file -> [meta.subMap("id"), file] }
        c = STAR_ALIGN_DCC_R2.out.chimeric_junctions.map { meta, file -> [meta.subMap("id"), file] }
        dcc_main_ch = a.join(b).join(c).join(tabs)

        DCC_MAIN(dcc_main_ch, reference, refseq_gtf)
        BEDTOOLS_WINDOW_DCC(DCC_MAIN.out.counts, refseq_bed)
        outreader_ch = BEDTOOLS_WINDOW_DCC.out.windows.join(DCC_MAIN.out.coordinates)
        DCC_OUTREADER(outreader_ch)
        matrixmaker_dcc_channel = DCC_OUTREADER.out.processed_circs.map { meta, file -> file }.collect().map{ files -> [["id":params.output_basename], files] }
        CIRCSNAKE_MATRIXMAKER_DCC(matrixmaker_dcc_channel, mm1_circ, mm1_refseq)
        CIRCSNAKE_MATRIXTWO_DCC(CIRCSNAKE_MATRIXMAKER_DCC.out.matrix, micrornas, coding_circs, hallmarks, ensembl_gene_descriptions)
    }

    if (params.run_cx) {
        log.info "Running CX"

        STAR_ALIGN_CX(SAMTOOLS_FASTQ.out.fastq, star_genome)
        chimeric_junc_channel = STAR_ALIGN_CX.out.chimeric_junctions.map { meta, file -> [meta.subMap("id"), file] }

        CIRCEXPLORER_STARPARSE(chimeric_junc_channel)

        CIRCEXPLORER_MAIN(CIRCEXPLORER_STARPARSE.out.fusion_junctions, reference, refseq_annot)
        CIRCEXPLORER_OUTREADER(CIRCEXPLORER_MAIN.out.circs)
        matrixmaker_cx_channel = CIRCEXPLORER_OUTREADER.out.processed_circs.map { meta, file -> file }.collect().map{ files -> [["id":params.output_basename], files] }
        CIRCSNAKE_MATRIXMAKER_CX(matrixmaker_cx_channel, mm1_circ, mm1_refseq)
        CIRCSNAKE_MATRIXTWO_CX(CIRCSNAKE_MATRIXMAKER_CX.out.matrix, micrornas, coding_circs, hallmarks, ensembl_gene_descriptions)
    }

    if (params.run_fc) {
        log.info "Running FC"

        BOWTIE2_ALIGN_PRIMARY(SAMTOOLS_FASTQ.out.fastq, bowtie_index, cram_reference, channel.value(false), channel.value(true))
        primary_alignment_channel = BOWTIE2_ALIGN_PRIMARY.out.bam.join(BOWTIE2_ALIGN_PRIMARY.out.csi)
        SAMTOOLS_VIEW(primary_alignment_channel, channel.value([[],[]]), channel.value([]), channel.value([]))
        FINDCIRC_ANCHORS(SAMTOOLS_VIEW.out.bam)
        SAMTOOLS_SPLITFASTA(reference.merge(reference_index))
        FINDCIRC_MAIN(FINDCIRC_ANCHORS.out.anchors.map { meta, file -> [meta + ["single_end": true], file]}, bowtie_index, SAMTOOLS_SPLITFASTA.out.split_fasta, channel.value(false))
        FINDCIRC_FILTER(FINDCIRC_MAIN.out.bed_ci.map { meta, file -> [meta.subMap('id'), file]})
        BEDTOOLS_WINDOW_FC(FINDCIRC_FILTER.out.candidates, refseq_bed)
        FINDCIRC_OUTREADER(BEDTOOLS_WINDOW_FC.out.windows)
        matrixmaker_fc_channel = FINDCIRC_OUTREADER.out.processed_circs.map { meta, file -> file }.collect().map{ files -> [["id":params.output_basename], files] }
        CIRCSNAKE_MATRIXMAKER_FC(matrixmaker_fc_channel, mm1_circ, mm1_refseq)
        CIRCSNAKE_MATRIXTWO_FC(CIRCSNAKE_MATRIXMAKER_FC.out.matrix, micrornas, coding_circs, hallmarks, ensembl_gene_descriptions)
    }

    CIRCSNAKE_VOTE(CIRCSNAKE_MATRIXTWO_FC.out.matrix.join(CIRCSNAKE_MATRIXTWO_CX.out.matrix).join(CIRCSNAKE_MATRIXTWO_DCC.out.matrix))
    votes = CIRCSNAKE_VOTE.out.circex_voted.mix(CIRCSNAKE_VOTE.out.find_circ_voted, CIRCSNAKE_VOTE.out.dcc_voted)
    norm_channel = votes.combine(CIRCSNAKE_READCOUNTS.out.readcounts)
    CIRCSNAKE_NORM(norm_channel)
}
