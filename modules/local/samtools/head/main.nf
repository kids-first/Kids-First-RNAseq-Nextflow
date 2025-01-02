process SAMTOOLS_HEAD {
    label 'process_low'
    container "staphb/samtools:1.20"

    input:
    path(input_align)
    val line_filter

    output:
    path '*.txt', emit: header_file

    script:
    def grep_line =
        line_filter != "" ? "| grep $line_filter"
        : ''

    """
    samtools \\
    head \\
    $input_align \\
    $grep_line \\
    > ${input_align.getBaseName()}.txt
    """

    stub:
    """
    touch header.txt
    """
}
