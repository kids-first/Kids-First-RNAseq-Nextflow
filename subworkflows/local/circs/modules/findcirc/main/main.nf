process FINDCIRC_MAIN {
    label 'C16'
    container "pgc-images.sbgenomics.com/danmiller/findcirc_mp:0.1.0"

    input:
    tuple val(meta) , path(reads)
    path(index_tar)
    path(reference)
    val   save_unaligned

    output:
    tuple val(meta), path("*.sam")      , emit: sam     , optional:true
    tuple val(meta), path("*.bam")      , emit: bam     , optional:true
    tuple val(meta), path("*.cram")     , emit: cram    , optional:true
    tuple val(meta), path("*.csi")      , emit: csi     , optional:true
    tuple val(meta), path("*.crai")     , emit: crai    , optional:true
    tuple val(meta), path("*fastq.gz")  , emit: fastq   , optional:true
    path("*.secondpass.log")            , emit: bowtie_log
    path("*.stats.txt")         , emit: stats
    tuple val(meta), path("*.sites.bed")                 , emit: bed_ci
    tuple val(meta), path("*.sites.reads")               , emit: reads_ci
    path("worker_error*.log")  , emit: fc_worker_error   , optional:true
    path  "versions.yml"                , emit: versions

    script:
    def args = task.ext.args ?: ""
    def args2 = task.ext.args2 ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"

    def unaligned = ""
    def reads_args = ""
    if (meta.single_end) {
        unaligned = save_unaligned ? "--un-gz ${prefix}.unmapped.fastq.gz" : ""
        reads_args = "-U ${reads}"
    } else {
        unaligned = save_unaligned ? "--un-conc-gz ${prefix}.unmapped.fastq.gz" : ""
        reads_args = "-1 ${reads[0]} -2 ${reads[1]}"
    }
    def sample_name = meta.sample_name ?: "${meta.id}"
    """
    tar --use-compress-program=pigz -xvf $index_tar

    INDEX=`find -L ./ -name "*.rev.1.bt2" | sed "s/\\.rev.1.bt2\$//"`
    [ -z "\$INDEX" ] && INDEX=`find -L ./ -name "*.rev.1.bt2l" | sed "s/\\.rev.1.bt2l\$//"`
    [ -z "\$INDEX" ] && echo "Bowtie2 index files not found" 1>&2 && exit 1

    bowtie2 \\
        -x \$INDEX \\
        --reorder \\
        $reads_args \\
        --threads 6 \\
        $unaligned \\
        $args \\
        2>| >(tee ${prefix}.secondpass.log >&2) \\
    | find_circ_mp.py \\
        -G $reference \\
        -p $sample_name \\
        --output_prefix ${prefix} \\
        --procs ${task.cpus - 6} \\
        $args2

    if [ -f ${prefix}.unmapped.fastq.1.gz ]; then
        mv ${prefix}.unmapped.fastq.1.gz ${prefix}.unmapped_1.fastq.gz
    fi

    if [ -f ${prefix}.unmapped.fastq.2.gz ]; then
        mv ${prefix}.unmapped.fastq.2.gz ${prefix}.unmapped_2.fastq.gz
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """
}
