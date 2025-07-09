process RNASEQC {
    label 'C4'
    container 'pgc-images.sbgenomics.com/d3b-bixu/rnaseqc:v2.4.2'

    input:
    path(collapsed_gtf)
    tuple path(aligned_sorted_reads), path(bai)
    val(stranded)
    val(is_paired)

    output:
    path('output/*metrics.tsv'), emit: Metrics
    path('output/*gene_tpm.gct'), emit: Gene_TPM
    path('output/*gene_reads.gct'), emit: Gene_count
    path('output/*exon_reads.gct'), emit: Exon_count

    script:
    """
    rnaseqc \\
    $collapsed_gtf \\
    $aligned_sorted_reads \\
    output \\
    --legacy \\
    ${stranded != "" ? "--stranded=$stranded" : ""} \\
    ${!is_paired ? "--unpaired" : ""}
    """

}