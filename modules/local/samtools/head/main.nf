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

workflow samtools_head {
    take:
    input_align
    line_filter
    
    main:
    SAMTOOLS_HEAD(input_align, line_filter)
    emit:
    SAMTOOLS_HEAD.out
}

workflow  {
    input_align = Channel.fromPath(params.input_alignment_reads)
    line_filter = params.line_filter
    samtools_head(input_align, line_filter)
    
}