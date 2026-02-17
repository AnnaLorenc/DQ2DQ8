include { MERGE_SAMPLES as MERGE_AMINO_VFAM } from '../modules/merge_samples.nf'
include { SUBSAMPLE_FROM_MERGED } from '../modules/subsample_from_merged.nf'

workflow MERGING_AND_SUBSAMPLING {
    take:
    productive_ch
    index_columns
    N
    M

    main:
    // Prepare merge channel
    prod_merge_ch = productive_ch
        .map { sample_name, file_path -> file_path }
        .collect()
        .map { files -> 
            tuple(files, index_columns, "count (templates/reads)", "productive") 
        }

    // Merge sequences
    MERGE_AMINO_VFAM(prod_merge_ch)

    // Prepare subsample channel
    subsample_input_ch = MERGE_AMINO_VFAM.out.merged
        .map { merged_file -> tuple(merged_file, N, M, 3, "productive") }

    // Subsample sequences
    SUBSAMPLE_FROM_MERGED(subsample_input_ch)

    emit:
    merged = MERGE_AMINO_VFAM.out.merged
    subsampled = SUBSAMPLE_FROM_MERGED.out.subsampled
    sample_names = SUBSAMPLE_FROM_MERGED.out.sample_names  // if available
    subsampled_dir = SUBSAMPLE_FROM_MERGED.out.subsampled_dir  // if available
}