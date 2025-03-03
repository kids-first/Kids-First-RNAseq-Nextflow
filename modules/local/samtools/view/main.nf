process SAMTOOLS_VIEW{
    label 'process_rsem'
    container 'pgc-images.sbgenomics.com/d3b-bixu/samtools:1.9'

    input:
    tuple path(reference), path(fai)
    tuple path(input_bam), path(bai)

    output:
    tuple path('*cram'), path('*crai'), emit: cram

    script:
    """
    samtools view \\
    -C \\
    -o ${input_bam.getBaseName()}.cram \\
    -T $reference \\
    -@ 16 \\
    && samtools \\
    index \\
    ${input_bam.getBaseName()}.cram
    """

}