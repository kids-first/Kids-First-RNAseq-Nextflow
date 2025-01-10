process SAMTOOLS_HEAD_RG_CT {
    label 'process_low'
    container "staphb/samtools:1.20"

    input:
    path(input_align)
    val line_filter

    output:
    tuple(stdout, path(input_align))

    script:
    def grep_line =
        line_filter != "" ? "| grep $line_filter"
        : ''

    """
    samtools \\
    head \\
    $input_align \\
    $grep_line \\
    | wc -l \\
    | tr -d "\n"
    """

}
