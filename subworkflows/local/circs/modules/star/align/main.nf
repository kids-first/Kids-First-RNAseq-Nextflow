process STAR_ALIGN {
    label 'M16'
    container "pgc-images.sbgenomics.com/d3b-bixu/star:2.7.10a"

    input:
    tuple val(meta), path(reads)
    path(genomeDir)

    output:
    tuple val(meta), path('*Log.progress.out'), emit: log_progress_out
    tuple val(meta), path('*Log.out'), emit: log_out
    tuple val(meta), path('*Log.final.out'), emit: log_final_out
    tuple val(meta), path('*Aligned.out.bam'), optional:true, emit: genomic_bam_out
    tuple val(meta), path('*Aligned.out.sam'), optional:true, emit: genomic_sam_out
    tuple val(meta), path('*SJ.out.tab.gz'), emit: junctions_out
    tuple val(meta), path('*Aligned.toTranscriptome.out.bam'), optional:true, emit: transcriptome_bam_out
    tuple val(meta), path('*Chimeric.out.sam'), optional: true, emit: chimeric_sam_out
    tuple val(meta), path('*Chimeric.out.junction'), emit: chimeric_junctions
    tuple val(meta), path('*ReadsPerGene.out.tab.gz'), optional:true, emit: gene_counts

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def star_ext_args = task.ext.args ?: ''
    def reads1 = []
    def reads2 = []
    meta.is_paired_end ? reads.eachWithIndex{ v, ix -> ( ix & 1 ? reads2 : reads1) << v } : [reads].flatten().each{ read -> reads1 << read}
    def read_command = reads1[0].name.endsWith(".gz") ? "zcat" : "cat"
    """
    tar -I pigz -xvf $genomeDir \\
    && STAR \\
    --runThreadN $task.cpus \\
    --genomeDir ./${genomeDir.getBaseName().replace(".tar", "")} \\
    --readFilesIn ${reads1.join(",")} ${reads2.join(",")} \\
    --readFilesCommand $read_command \\
    --outFileNamePrefix $prefix. \\
    $star_ext_args \\
    && pigz *tab \\
    && rm -rf ./${genomeDir.getBaseName().replace(".tar", "")}
    """

}
