process SAMTOOLS_SPLIT {
    label 'process_medium'
    container "staphb/samtools:1.20"

    input:
    path(input_reads)
    path(reference)
    val(threads)

    output:
    tuple path('*.bam'), optional: true, emit: bam_files

    script:
    """
    RG_NUM=`samtools head $input_reads | grep -c ^@RG`
    if [ \$RG_NUM != 1 ]; then
      samtools split -f '%*_%#.bam' -@ $threads --reference $reference $input_reads
    fi
    """

    stub:
    """
    touch test.bam
    """
}


workflow samtools_split {
    take:
    input_reads
    reference
    threads
    main:
    SAMTOOLS_SPLIT(input_reads, reference, threads)
    emit:
    SAMTOOLS_SPLIT.out
}

workflow {
    input_reads = Channel.fromPath(params.unsorted_bam)
    reference = Channel.fromPath(params.reference)
    threads = Channel.from(params.threads)
    samtools_split(input_reads, reference, threads)

}