process ANNOFUSE {
    label 'C4'
    container "pgc-images.sbgenomics.com/d3b-bixu/annofuse:0.92.0"

    input:
    path(arriba_formatted_fusions)
    path(starfusion_formatted_fusions)
    path(rsem_expr_file)
    val(sample_name)
    val(output_basename)

    output:
    path('*annoFuse_filter.tsv'), optional: true, emit: filtered_fusions_tsv

    script:
    """
    A_CT=`wc -l $arriba_formatted_fusions | cut -f 1 -d " "`

    S_CT=`wc -l $starfusion_formatted_fusions | cut -f 1 -d " "`

    if [ \$A_CT -eq 1 ] && [ \$S_CT -eq 1 ]; then
    echo "Both inputs are empty, will skip processing as there no fusions." >&2;
    exit 0;
    fi

    annoFusePerSample.R \\
    --fusionfileArriba $arriba_formatted_fusions \\
    --fusionfileStarFusion $starfusion_formatted_fusions \\
    --expressionFile $rsem_expr_file \\
    --tumorID $sample_name \\
    --outputfile ${output_basename}.annoFuse_filter.tsv

    """

}