process ANNOTATE_ARRIBA {
    label 'C4'
    container "pgc-images.sbgenomics.com/d3b-bixu/fusionannotator:0.1.1"

    input:
    path(input_fusion_file)
    path(fusion_annotator_tar)

    output:
    path('*annotated.tsv'), emit: annotated_tsv

    script:
    def ext_args = task.ext.args ?: ''
    """
    tar -I pigz -xf $fusion_annotator_tar \\
    && /opt/FusionAnnotator/FusionAnnotator \\
    --annotate $input_fusion_file \\
    $ext_args \\
    > ${task.prefix}.annotated.tsv \\

    """

}