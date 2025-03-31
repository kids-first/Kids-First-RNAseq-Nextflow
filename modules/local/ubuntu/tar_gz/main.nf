process TAR_GZ {
    label 'C2'
    container 'ubuntu:22.04'

    input:
    path(Gene_TPM)
    path(Gene_count)
    path(Exon_count)

    output:
    path('*tar.gz'), emit: RNASeQC_counts

    script:
    """
    mkdir ${task.ext.prefix}_RNASeQC_counts

    cp $Gene_TPM \\
    $Gene_count \\
    $Exon_count \\
    ${task.ext.prefix}_RNASeQC_counts

    tar -czf \\
    ${task.ext.prefix}.RNASeQC.counts.tar.gz \\
    ${task.ext.prefix}_RNASeQC_counts

    """

}