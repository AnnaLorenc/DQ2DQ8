process MERGE_SAMPLES {

    publishDir "${params.outdir}/merged/", mode: 'copy'

    input:
    tuple val(files), val(merge_columns), val(count_column), val(type)
    
    output: 
    path "merged_${merge_columns.join('_')}_${type}.tsv.gz", emit: merged    

    script:
    def merge_cols_arg = merge_columns.join(' ')
    def merge_cols_filename = merge_columns.join('_')
    def input_files_arg = files.join(' ')
    
    """
    python ${projectDir}/bin/merge_trb_data.py \
        --input_files ${input_files_arg} \
        --merge_columns ${merge_cols_arg} \
        --count_column "${count_column}" \
        --output_file merged_${merge_cols_filename}_${type}.tsv.gz
    """
}