process ARRIBA_DRAW {
    label 'process_arriba_draw'
    container "pgc-images.sbgenomics.com/d3b-bixu/arriba:2.2.1"

    input:
    tuple path(genome_aligned_bam), path(genome_aligned_bai)
    path(fusions)
    path(gtf_anno)

    output:
    path('*pdf'), emit: arriba_pdf

    script:
    def arriba_ext_args = task.ext.args ?: ''
    """
    /arriba_v2.2.1/draw_fusions.R \\
    --fusions=$fusions \\
    --alignments=$genome_aligned_bam \\
    --annotation=$gtf_anno \\
    --output=${fusions.getBaseName()}.pdf \\
    $arriba_ext_args
    """

}