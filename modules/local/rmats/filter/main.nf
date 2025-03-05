process AWK_JC_FILTER {
    label 'C2'
    container "ubuntu:22.04"

    input:
    path(input_jc_file)

    output:
    path('*filtered.*MATS.JC.txt'), emit: output

    script:
    def output_jc_file = input_jc_file.name.tokenize(".")
    def ins_pos = output_jc_file.size() - 4
    output_jc_file.add(ins_pos, "filtered")
    output_jc_file = output_jc_file.join(".")
    """
    awk -F'\t' 'NR==1 { for (i=1; i<=NF; i++) { f[\$i] = i } print \$line } \$f["IJC_SAMPLE_1"]+\$f["IJC_SAMPLE_2"]+\$f["SJC_SAMPLE_1"]+\$f["SJC_SAMPLE_2"] >= 10 {print \$line}' $input_jc_file > $output_jc_file;
    """

}