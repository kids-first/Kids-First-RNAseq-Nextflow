process SAMTOOLS_SPLIT {
    label 'process_medium'
    container "staphb/samtools:1.20"

    input:
    path(input_reads)
    path(reference)
    val(threads)

    output:
    path('*.bam'), optional: true, emit: bam_files

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
