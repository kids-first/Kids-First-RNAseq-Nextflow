process SAMTOOLS_SPLIT {
    label 'process_medium'
    container "staphb/samtools:1.20"

    input:
    tuple val(meta), path(input_reads)
    path(reference)
    val(threads)

    output:
    tuple val(meta), path('*.bam'), emit: bams

    script:
    """
    samtools split -f '%*_%#.bam' -@ $threads --reference $reference $input_reads
    """

    stub:
    """
    touch test.bam
    """
}
