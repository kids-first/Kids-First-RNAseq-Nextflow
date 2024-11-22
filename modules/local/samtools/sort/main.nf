process SAMTOOLS_SORT {
    tag "$meta.id"
    label 'process_extra_medium'
    container "pgc-images.sbgenomics.com/d3b-bixu/samtools:1.9"

    input:
    tuple val(meta), path(unsorted_bam)


    output:
    tuple val(meta), path('*.bam'), emit: sorted_bam
    tuple val(meta), path('*.bai'), emit: sorted_bai

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools sort \\
    $unsorted_bam \\
    --threads task.cpus \\
    -m 1G \\
    -O bam
    > ${prefix}.sorted.bam
    && samtools index
    -@ task.cpus
    ${prefix}.sorted.bam
    ${prefix}.sorted.bai
    """

    stub:
    """
    touch test.bam
    touch test.bai
    """
}