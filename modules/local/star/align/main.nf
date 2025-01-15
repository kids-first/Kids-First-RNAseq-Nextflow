process STAR_ALIGN {
    label 'process_star'
    container "pgc-images.sbgenomics.com/d3b-bixu/star:2.7.10a"

    input:
    path(genomeDir)
    val(readFilesCommand)
    val(readFilesManifest)
    val(outFileNamePrefix)

    output:
    path('*Aligned.out.bam'), emit: genomic_bam_out

    script:
    """
    tar -I pigz -xvf $genomeDir \\
    && STAR \\
    --genomeDir ./${genomeDir.getName().replace(".tar", "")} \\
    --readFilesCommand $readFilesCommand \\
    --readFilesManifest $readFilesManifest \\
    --outFileNamePrefix $outFileNamePrefix

    """

}