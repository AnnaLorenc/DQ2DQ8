include { COMPUTE_OVERLAPS } from '../modules/compute_overlaps.nf'
include { COMPUTE_OVERLAPS_SUBSAMPLED } from '../modules/compute_overlaps_subsampled.nf'

workflow OVERLAP_ANALYSIS {
    take:
    merged_ch
    subsampled_data_ch  // This could be sample_names + subsampled_dir
    index_columns

    main:
    def first_column_with_sample_in_merged = index_columns.size()

    // Compute overlaps on merged data
    compute_overlaps_input_ch = merged_ch
        .map { merged_file -> tuple(merged_file, "productive", first_column_with_sample_in_merged) }
        
    COMPUTE_OVERLAPS(compute_overlaps_input_ch)

    // Compute overlaps on subsampled data
    COMPUTE_OVERLAPS_SUBSAMPLED(subsampled_data_ch)

    emit:
    merged_overlaps = COMPUTE_OVERLAPS.out
    subsampled_overlaps = COMPUTE_OVERLAPS_SUBSAMPLED.out
}