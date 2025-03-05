process RSEM {
    label 'C16'
    container "images.sbgenomics.com/uros_sipetic/rsem:1.3.1"

    input:
    path(rsem_genome_tar)
    path(transcriptome_bam)
    val(outFileNamePrefix)
    val(pairedness)
    val(strandedness)

    output:
    path('*genes.results.gz'), emit: gene_out
    path('*isoforms.results.gz'), emit: isoform_out

    script:
    def rsem_ext_args = task.ext.args ?: ''
    """
    tar xzf $rsem_genome_tar \\
    && rsem-calculate-expression \\
    --no-bam-output \\
    --alignments \\
    --strandedness $strandedness \\
    ${pairedness ? "--paired-end" : ""} \\
    $rsem_ext_args \\
    $transcriptome_bam \\
    ./${rsem_genome_tar.getBaseName().replace(".tar", "")}/${rsem_genome_tar.getBaseName().replace(".tar", "")} \\
    ${outFileNamePrefix}.rsem \\
    && gzip *results
    """

}