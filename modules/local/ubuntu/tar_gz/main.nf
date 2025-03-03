process TAR_GZ {
    label 'process_annofuse'
    container 'ubuntu:22.04'

    input:
    path(Gene_TPM)
    path(Gene_count)
    path(Exon_count)
    val(outFileNamePrefix)

    output:
    path('*tar.gz'), emit: RNASeQC_counts

    script:
    """
    mkdir ${outFileNamePrefix}_RNASeQC_counts

    cp $Gene_TPM \\
    $Gene_count \\
    $Exon_count \\
    ${outFileNamePrefix}_RNASeQC_counts

    tar -czf \\
    ${outFileNamePrefix}.RNASeQC.counts.tar.gz \\
    ${outFileNamePrefix}_RNASeQC_counts

    """

}