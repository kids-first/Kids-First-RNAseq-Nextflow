process SAMTOOLS_SPLIT {
    label 'C8'
    container "staphb/samtools:1.20"

    input:
    tuple val(meta), path(input_reads)
    path(reference)

    output:
    tuple val(meta), path('*.bam'), emit: bams

    script:
    """
    samtools split -f '%*_%#.bam' -@ $task.cpus --reference $reference $input_reads
    """

    stub:
    """
    touch test.bam
    """
}
