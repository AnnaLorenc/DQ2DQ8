include { COMPUTE_FREQS_SUBSAMPLED } from '../modules/compute_freqs_subsampled.nf'
include { COMBINE_FREQS_SUBSAMPLED } from '../modules/combine_freqs_subsampled.nf'

workflow SUBSAMPLED_ANALYSIS {
    take:
    subsampled_ch

    main:
    // Prepare input for frequency computation
    compute_freqs_subsampled_input_ch = subsampled_ch
        .flatten()
        .map { subsampled_file -> 
            def sample_name = subsampled_file.baseName.replaceAll(/_subsampled\.tsv$/, '')
            tuple(subsampled_file, sample_name, "productive", "aminoAcid", ["vFamilyName","jGeneName"])
        }
    
    // Compute frequencies on subsampled data
    COMPUTE_FREQS_SUBSAMPLED(compute_freqs_subsampled_input_ch)

    // Combine subsampled frequencies
    combine_freqs_subsampled_input_ch = COMPUTE_FREQS_SUBSAMPLED.out.subsampled
        .map { big_file, summ_file -> summ_file }
        .collect()
        .map { files ->
            tuple(files, "productive_subs_freqs.tsv")
        }

    COMBINE_FREQS_SUBSAMPLED(combine_freqs_subsampled_input_ch)

    emit:
    combined_freqs = COMBINE_FREQS_SUBSAMPLED.out
    individual_freqs = COMPUTE_FREQS_SUBSAMPLED.out.subsampled
}