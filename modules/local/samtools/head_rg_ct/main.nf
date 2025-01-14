process SAMTOOLS_HEAD_RG_CT {
    label 'process_low'
    container "staphb/samtools:1.20"

    input:
    tuple val(meta), path(input_align)

    output:
    tuple val(meta), env('RG_NUM'), path(input_align), emit: reads

    script:
    """
    RG_NUM=`samtools head $input_align | grep -c ^@RG`
    """

}
