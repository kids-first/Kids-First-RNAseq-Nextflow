process FINDCIRC_FILTER {
    label 'process_single'
    container "pgc-images.sbgenomics.com/danmiller/circs-dcc:0.1.0"

    input:
    tuple val(meta), path(sites)

    output:
    tuple val(meta), path('*.circ_candidates_auto.bed'), emit: candidates

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    grep circ $sites \\
    | grep -v chrM \\
    | python /opt/circs_snake/scripts/pipelines/sum.py -2,3 \\
    | python /opt/circs_snake/scripts/pipelines/scorethresh.py -15 1 \\
    | python /opt/circs_snake/scripts/pipelines/scorethresh.py -14 2 \\
    | python /opt/circs_snake/scripts/pipelines/scorethresh.py -13 2 \\
    | python /opt/circs_snake/scripts/pipelines/scorethresh.py 7 2 \\
    | python /opt/circs_snake/scripts/pipelines/scorethresh.py 9,10 35 \\
    | python /opt/circs_snake/scripts/pipelines/scorethresh.py -19 100000 \\
    > ${prefix}.circ_candidates_auto.bed
    """
}
