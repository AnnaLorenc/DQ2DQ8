process merge_amino_vfam {

    input:
    tuple file_info from input_channel // Expects a tuple: (file_names, merge_columns, count_column)

    output:
    path "merged_output.tsv.gz" into output_channel

    script:
    """
    python ./bin/merge_trb_data.py \
        --input_files ${file_info[0].join(' ')} \
        --merge_columns ${file_info[1].join(' ')} \
        --count_column "${file_info[2]}" \
        --output_file merged_output.tsv.gz
    """
}