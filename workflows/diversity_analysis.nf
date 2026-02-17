include { EXTRACT_DIVERSITY_METRICS as EXTRACT_DIVERSITY_METRICS_PROD } from '../modules/extract_diversity_metrics.nf'
include { EXTRACT_DIVERSITY_METRICS as EXTRACT_DIVERSITY_METRICS_NONPROD } from '../modules/extract_diversity_metrics.nf'
include { COMBINE_DIVERSITY_METRICS as COMBINE_DIVERSITY_METRICS_PROD } from '../modules/combine_diversity_metrics.nf'
include { COMBINE_DIVERSITY_METRICS as COMBINE_DIVERSITY_METRICS_NONPROD } from '../modules/combine_diversity_metrics.nf'

workflow DIVERSITY_ANALYSIS {
    take:
    productive_ch
    nonproductive_ch

    main:
    // Extract diversity metrics
    EXTRACT_DIVERSITY_METRICS_PROD(productive_ch)
    EXTRACT_DIVERSITY_METRICS_NONPROD(nonproductive_ch)

    // Collect and combine productive data
    productive_data = EXTRACT_DIVERSITY_METRICS_PROD.out.diversity_metrics
        .unique()
        .toList()
        .map { items -> 
            def samples = items.collect { it[0] }
            def files = items.collect { it[1] }
            tuple('productive', samples, files)
        }
    
    // Collect and combine nonproductive data
    nonproductive_data = EXTRACT_DIVERSITY_METRICS_NONPROD.out.diversity_metrics
        .unique()
        .toList()
        .map { items -> 
            def samples = items.collect { it[0] }
            def files = items.collect { it[1] }
            tuple('nonproductive', samples, files)
        }
    
    // Combine diversity metrics
    COMBINE_DIVERSITY_METRICS_PROD(productive_data)
    COMBINE_DIVERSITY_METRICS_NONPROD(nonproductive_data)

    emit:
    productive_combined = COMBINE_DIVERSITY_METRICS_PROD.out
    nonproductive_combined = COMBINE_DIVERSITY_METRICS_NONPROD.out
}