#!/usr/bin/env nextflow


nextflow.enable.dsl=2

// Parameters
params.input = './data/collated_info.csv'
params.outdir = './results'
params.seqs_to_remove = './data/D20210208D_1-overrepresented-sequences.txt'


include {INITIAL_CLEANUP_SPLIT } from './modules/read_sheet.nf'
include { EXTRACT_DIVERSITY_METRICS as EXTRACT_DIVERSITY_METRICS_PROD } from './modules/extract_diversity_metrics.nf'
include { EXTRACT_DIVERSITY_METRICS as EXTRACT_DIVERSITY_METRICS_NONPROD } from './modules/extract_diversity_metrics.nf'
include { COMPUTE_FREQS as COMPUTE_FREQS_PROD} from './modules/compute_freqs_metrics.nf'
include { COMPUTE_FREQS as COMPUTE_FREQS_NONPROD} from './modules/compute_freqs_metrics.nf'
include{ COMBINE_DIVERSITY_METRICS as COMBINE_DIVERSITY_METRICS_PROD } from './modules/combine_diversity_metrics.nf'
include{ COMBINE_DIVERSITY_METRICS as COMBINE_DIVERSITY_METRICS_NONPROD } from './modules/combine_diversity_metrics.nf'
include{ COMBINE_FREQS as COMBINE_FREQS_PROD } from './modules/combine_freqs.nf'
include{ COMBINE_FREQS as COMBINE_FREQS_NONPROD } from './modules/combine_freqs.nf'


// Print a header for your pipeline 
log.info """\

=======================================================================================
DQ2DQ8 
=======================================================================================

Created by Ania


=======================================================================================
Workflow run parameters 
=======================================================================================
input       : ${params.input}
results     : ${params.outdir}
workDir     : ${workflow.workDir}
=======================================================================================

"""


workflow {

// Show help message if --help is run or (||) a required parameter (input) is not provided
		
		log.info "seqs_to_remove: ${params.seqs_to_remove}"
		log.info "DEBUG: params.input = '${params.input}'"  // Add this debug line
    	log.info "DEBUG: params.input type = ${params.input?.getClass()}"  // Check type
		// Channel 1: sample_short + genotype_short
		channel.fromPath(params.input)
			.splitCsv(header: true, sep: ',')
			.map { row -> tuple(row.sample_short, row.genotype_short) }
			.set { sample_genotype_ch }

		// Channel 2: sample_short + filepath
		channel.fromPath(params.input)
			.splitCsv(header: true, sep: ',')
			.map { row -> tuple(row.sample_short, "${row.LOC}/${row.SAMPLE}", params.seqs_to_remove, params.outdir) }
			.set { sample_filepath_ch}

		// Use channels
		
		INITIAL_CLEANUP_SPLIT (sample_filepath_ch)

		// INITIAL_CLEANUP_SPLIT.out.productive.view { sample, file -> "Sample ${sample}: ${file}"    }
		EXTRACT_DIVERSITY_METRICS_PROD(INITIAL_CLEANUP_SPLIT.out.productive)
		EXTRACT_DIVERSITY_METRICS_NONPROD(INITIAL_CLEANUP_SPLIT.out.nonproductive)
		COMPUTE_FREQS_PROD(INITIAL_CLEANUP_SPLIT.out.productive)
		COMPUTE_FREQS_NONPROD(INITIAL_CLEANUP_SPLIT.out.nonproductive)

		// Collect productive CSV files
		productive_data = EXTRACT_DIVERSITY_METRICS_PROD.out.diversity_metrics
			.unique()  // Remove duplicates
			.toList()
			.map { items -> 
				def samples = items.collect { it[0] }
				def files = items.collect { it[1] }
				tuple('productive', samples, files)
			}
		
		// Collect nonproductive with deduplication
		nonproductive_data = EXTRACT_DIVERSITY_METRICS_NONPROD.out.diversity_metrics
			.unique()  // Remove duplicates
			.toList()
			.map { items -> 
				def samples = items.collect { it[0] }
				def files = items.collect { it[1] }
				tuple('nonproductive', samples, files)
			}
		
		// Combine each type
		COMBINE_DIVERSITY_METRICS_PROD(productive_data)
		COMBINE_DIVERSITY_METRICS_NONPROD(nonproductive_data)


		// Collect productive CSV files for freqs
		productive_freqs = COMPUTE_FREQS_PROD.out.freqs
		.unique()  // Remove duplicates
		.toList()
		.map { items -> 
			def samples = items.collect { it[0] }
			def files = items.collect { it[1] }
			tuple('productive', samples, files)
		}
	
	// Collect nonproductive with deduplication
		nonproductive_freqs = COMPUTE_FREQS_NONPROD.out.freqs
			.unique()  // Remove duplicates
			.toList()
			.map { items -> 
				def samples = items.collect { it[0] }
				def files = items.collect { it[1] }
				tuple('nonproductive', samples, files)
			}

	COMBINE_FREQS_PROD(productive_freqs)
	COMBINE_FREQS_NONPROD(nonproductive_freqs)

}

// Print workflow execution summary 
workflow.onComplete {
summary = """
=======================================================================================
Workflow execution summary
=======================================================================================

Duration    : ${workflow.duration}
Success     : ${workflow.success}
workDir     : ${workflow.workDir}
Exit status : ${workflow.exitStatus}
results     : ${params.outdir}

=======================================================================================
  """
println summary

}
