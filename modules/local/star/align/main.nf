process STAR_ALIGN {
    label 'process_star'
    container "pgc-images.sbgenomics.com/d3b-bixu/star:2.7.10a"

    input:
    path(genomeDir)
    val(readFilesCommand)
    val(readFilesManifest)
    val(outFileNamePrefix)

    output:
    path('*Log.progress.out'), emit: log_progress_out
    path('*Log.out'), emit: log_out
    path('*Log.final.out'), emit: log_final_out
    path('*Aligned.out.bam'), emit: genomic_bam_out
    path('*SJ.out.tab.gz'), emit: junctions_out
    path('*Aligned.toTranscriptome.out.bam'), emit: transcriptome_bam_out
    path('*Chimeric.out.sam'), optional: true, emit: chimeric_sam_out
    path('*Chimeric.out.junction'), emit: chimeric_junctions
    path('*ReadsPerGene.out.tab.gz'), emit: gene_counts

    script:
    def star_ext_args = task.ext.args ?: ''
    """
    tar -I pigz -xvf $genomeDir \\
    && STAR \\
    --genomeDir ./${genomeDir.getBaseName().replace(".tar", "")} \\
    --readFilesCommand $readFilesCommand \\
    --readFilesManifest $readFilesManifest \\
    --outFileNamePrefix "${outFileNamePrefix}." \\
    $star_ext_args \\
    && pigz .*.tab
    """

}