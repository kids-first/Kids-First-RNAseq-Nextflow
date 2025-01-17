process SAMTOOLS_SORT {
    container "pgc-images.sbgenomics.com/d3b-bixu/samtools:1.9"

    input:
    path(unsorted_bam)
    val(output_basename)
    val(threads)

    output:
    path('*.bam'), emit: sorted_bam
    path('*.bai'), emit: sorted_bai

    script:
    """
    samtools sort \\
    $unsorted_bam \\
    --threads $threads \\
    -m 1G \\
    -O bam \\
    > ${output_basename}.sorted.bam \\
    && samtools index \\
    -@ $threads \\
    ${output_basename}.sorted.bam \\
    ${output_basename}.sorted.bai
    """

    stub:
    """
    touch test.bam
    touch test.bai
    """
}