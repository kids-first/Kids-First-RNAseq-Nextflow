process SAMTOOLS_SORT {
    label 'C8'
    container "pgc-images.sbgenomics.com/d3b-bixu/samtools:1.9"

    input:
    path(unsorted_bam)
    val(output_basename)

    output:
    path('*.bam'), emit: sorted_bam
    path('*.bai'), emit: sorted_bai

    script:
    """
    samtools sort \\
    $unsorted_bam \\
    --threads $task.cpus \\
    -m 1G \\
    -O bam \\
    > ${output_basename}.sorted.bam \\
    && samtools index \\
    -@ $task.cpus \\
    ${output_basename}.sorted.bam \\
    ${output_basename}.sorted.bai
    """

    stub:
    """
    touch test.bam
    touch test.bai
    """
}