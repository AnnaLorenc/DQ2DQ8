include { INITIAL_CLEANUP_SPLIT } from '../modules/read_sheet.nf'

workflow DATA_PREPROCESSING {
    take:
    sample_filepath_ch

    main:
    INITIAL_CLEANUP_SPLIT(sample_filepath_ch)

    emit:
    productive = INITIAL_CLEANUP_SPLIT.out.productive
    nonproductive = INITIAL_CLEANUP_SPLIT.out.nonproductive
}