process SAMTOOLS_SORT {
    label 'C8'
    container "pgc-images.sbgenomics.com/d3b-bixu/samtools:1.9"

    input:
    path(unsorted_bam)

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
    > ${task.prefix}.sorted.bam \\
    && samtools index \\
    -@ $task.cpus \\
    ${task.prefix}.sorted.bam \\
    ${task.prefix}.sorted.bai
    """

    stub:
    """
    touch test.bam
    touch test.bai
    """
}