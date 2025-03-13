process ALIGNMENT_PAIREDNESS {
    label 'C8'
    container "quay.io/biocontainers/pysam:0.22.0--py310h41dec4a_0"

    input:
    tuple val(meta), path(input_reads)
    path input_reference
    val max_reads

    output:
    env(RESULT), emit: result

    script:
    """
    RESULT=`alignmentfile_pairedness.py \\
    --input_reads $input_reads \\
    --input_reference $input_reference \\
    --max_reads $max_reads \\
    --threads $task.cpus`

    if [ \$RESULT == 'ReadType:MIXED' ]
    then
        echo "Result was mixed, could not determine pairedness"
        exit 1
    elif [ \$RESULT == 'ReadType:PAIRED' ]
    then
        RESULT=true
    else
        RESULT=false
    fi

    """
}