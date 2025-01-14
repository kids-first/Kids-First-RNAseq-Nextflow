process SAMTOOLS_HEAD {
    label 'process_low'
    container "staphb/samtools:1.20"

    input:
    tuple val(meta), path(input_align)
    val line_filter

    output:
    tuple val(meta), stdout, path(input_align), emit: header_info

    script:
    def grep_line =
        line_filter != "" ? "| grep $line_filter"
        : ''

    """
    samtools \\
    head \\
    $input_align \\
    $grep_line
    """
}
