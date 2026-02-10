process CIRCSNAKE_READCOUNTS {
    label 'process_single'
    container "ubuntu:22.04"

    input:
    val(sample_count_strings)

    output:
    path('*_sample_readcounts.tsv'), emit: readcounts

    script:
    def prefix = task.ext.prefix ?: "all"
    def body = sample_count_strings.join('\\n')
    """
    printf 'sample_short\\treads_total\\n' > ${prefix}_sample_readcounts.tsv
    printf '$body' >> ${prefix}_sample_readcounts.tsv
    """
}
