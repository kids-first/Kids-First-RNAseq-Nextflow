process KALLISTO {
    label 'C8'
    container 'images.sbgenomics.com/uros_sipetic/kallisto:0.43.1'

    input:
    path(transcript_idx)
    val(strand)
    path(reads)
    val(std_dev)
    val(avg_frag_len)
    val(is_paired_end)

    output:
    path('*kallisto.abundance.tsv.gz'), emit: abundance_out

    script:
    def read_line = ""
    read_line += is_paired_end ? "" : "--single -l $avg_frag_len -s $std_dev "
    read_line += reads instanceof List ? reads.join(" ") : reads
    """
    kallisto quant \\
    -i $transcript_idx \\
    -o output \\
    -b 10 \\
    -t 8 \\
    ${strand != "" ? "--" + strand : ""} \\
    $read_line \\
    && mv output/abundance.tsv ${task.ext.prefix}.kallisto.abundance.tsv \\
    && gzip ${task.ext.prefix}.kallisto.abundance.tsv

    """

}