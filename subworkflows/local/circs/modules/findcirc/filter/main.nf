process FINDCIRC_FILTER {
    label 'process_single'
    container "pgc-images.sbgenomics.com/danmiller/circs-findcirc:0.1.0"

    input:
    tuple val(meta), path(sites)

    output:
    path('*.circ_candidates_auto.bed'), emit: candidates

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    grep circ $sites \\ 
    | grep -v chrM \\
    | python scripts/pipelines/sum.py -2,3 \\
    | python scripts/pipelines/scorethresh.py -16 1 \\
    | python scripts/pipelines/scorethresh.py -15 2 \\
    | python scripts/pipelines/scorethresh.py -14 2 \\
    | python scripts/pipelines/scorethresh.py 7 2 \\
    | python scripts/pipelines/scorethresh.py 8,9 35 \\
    | python scripts/pipelines/scorethresh.py -17 100000 \\
    > ${prefix}.circ_candidates_auto.bed
    """
}
