include { COMPUTE_FREQS as COMPUTE_FREQS_PROD } from '../modules/compute_freqs_metrics.nf'
include { COMPUTE_FREQS as COMPUTE_FREQS_NONPROD } from '../modules/compute_freqs_metrics.nf'
include { COMBINE_FREQS as COMBINE_FREQS_PROD } from '../modules/combine_freqs.nf'
include { COMBINE_FREQS as COMBINE_FREQS_NONPROD } from '../modules/combine_freqs.nf'

workflow FREQUENCY_ANALYSIS {
    take:
    productive_ch
    nonproductive_ch

    main:
    // Compute frequencies
    COMPUTE_FREQS_PROD(productive_ch)
    COMPUTE_FREQS_NONPROD(nonproductive_ch)

    // Collect productive frequencies
    productive_freqs = COMPUTE_FREQS_PROD.out.freqs
        .unique()
        .toList()
        .map { items -> 
            def samples = items.collect { it[0] }
            def files = items.collect { it[1] }
            tuple('productive', samples, files)
        }

    // Collect nonproductive frequencies
    nonproductive_freqs = COMPUTE_FREQS_NONPROD.out.freqs
        .unique()
        .toList()
        .map { items -> 
            def samples = items.collect { it[0] }
            def files = items.collect { it[1] }
            tuple('nonproductive', samples, files)
        }

    // Combine frequencies
    COMBINE_FREQS_PROD(productive_freqs)
    COMBINE_FREQS_NONPROD(nonproductive_freqs)

    emit:
    productive_combined = COMBINE_FREQS_PROD.out
    nonproductive_combined = COMBINE_FREQS_NONPROD.out
}