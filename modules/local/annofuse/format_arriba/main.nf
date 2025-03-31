process FORMAT_ARRIBA {
    label 'C2'
    container "pgc-images.sbgenomics.com/d3b-bixu/annofuse:0.92.0"

    input:
    path(arriba_fusion_file)

    output:
    path('*arriba_formatted.tsv'), emit: formatted_fusion_tsv

    script:
    """
    formatArribaFusionCalls.R \\
    --fusionfile $arriba_fusion_file \\
    --outputfile ${task.ext.prefix}.arriba_formatted.tsv
    """

}