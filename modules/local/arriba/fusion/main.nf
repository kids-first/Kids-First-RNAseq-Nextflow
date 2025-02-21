process ARRIBA_FUSION {
    label 'process_arriba_fusion'
    container "pgc-images.sbgenomics.com/d3b-bixu/arriba:2.2.1"

    input:
    path(genome_aligned_bam)
    path(reference_fasta)
    path(gtf_anno)
    val(outFileNamePrefix)
    val(arriba_strand_flag)
    val(assembly)

    output:
    path('*arriba_2.2.1.fusions.tsv'), emit: arriba_fusions

    script:
    def arriba_ext_args = task.ext.args ?: ''
    """
    /arriba_v2.2.1/arriba \\
    -x $genome_aligned_bam \\
    -a $reference_fasta \\
    -g $gtf_anno \\
    -s $arriba_strand_flag \\
    -o ${outFileNamePrefix}.arriba_2.2.1.fusions.tsv \\
    -O ${outFileNamePrefix}.arriba_2.2.1.discarded_fusions.tsv \\
    $arriba_ext_args \\
    """

}