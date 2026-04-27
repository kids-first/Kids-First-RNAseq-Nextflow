process CIRCSNAKE_READCOUNTS {
    label 'process_single'
    container "ubuntu:22.04"

    input:
    tuple val(meta), path(reads_fastq)

    output:
    path('*_sample_readcounts.tsv'), emit: readcounts

    script:
    def prefix = task.ext.prefix ?: "all"
    """
    total_lines="\$(zcat $reads_fastq | wc -l | tr -d '[:space:]')"
    reads_total="\$(( total_lines / 4 ))"

    printf 'sample_short\\treads_total\\n' > ${prefix}_sample_readcounts.tsv
    printf "%s\\t%s\\n" $meta.sample_name \$reads_total >> ${prefix}_sample_readcounts.tsv
    """
}
