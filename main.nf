#!/usr/bin/env nextflow


nextflow.enable.dsl=2

// Parameters
params.input = './data/collated_info.csv'
params.outdir = './results'
params.seqs_to_remove = './data/D20210208D_1-overrepresented-sequences.txt'
params.N = 10000 
params.M = 25
params.index_columns = ["aminoAcid", "vFamilyName", "jGeneName"] //to join samples in merging
params.min_len = 8
params.max_len = 25	// CDR3 length filters for IMGT freqs

include {INITIAL_CLEANUP_SPLIT } from './modules/read_sheet.nf'
include { EXTRACT_DIVERSITY_METRICS as EXTRACT_DIVERSITY_METRICS_PROD } from './modules/extract_diversity_metrics.nf'
include { EXTRACT_DIVERSITY_METRICS as EXTRACT_DIVERSITY_METRICS_NONPROD } from './modules/extract_diversity_metrics.nf'
include { COMPUTE_FREQS as COMPUTE_FREQS_PROD} from './modules/compute_freqs_metrics.nf'
include { COMPUTE_FREQS as COMPUTE_FREQS_NONPROD} from './modules/compute_freqs_metrics.nf'
include{ COMBINE_DIVERSITY_METRICS as COMBINE_DIVERSITY_METRICS_PROD } from './modules/combine_diversity_metrics.nf'
include{ COMBINE_DIVERSITY_METRICS as COMBINE_DIVERSITY_METRICS_NONPROD } from './modules/combine_diversity_metrics.nf'
include{ COMBINE_FREQS as COMBINE_FREQS_PROD } from './modules/combine_freqs.nf'
include{ COMBINE_FREQS as COMBINE_FREQS_NONPROD } from './modules/combine_freqs.nf'
include { MERGE_SAMPLES as MERGE_AMINO_VFAM } from './modules/merge_samples.nf'
include { SUBSAMPLE_FROM_MERGED } from './modules/subsample_from_merged.nf'
include {COMPUTE_FREQS_SUBSAMPLED} from './modules/compute_freqs_subsampled.nf'
include { COMBINE_FREQS_SUBSAMPLED } from './modules/combine_freqs_subsampled.nf'
include {COMPUTE_OVERLAPS} from './modules/compute_overlaps.nf'
include {COMPUTE_OVERLAPS_SUBSAMPLED} from './modules/compute_overlaps_subsampled.nf'
include {COMPUTE_IMGT_AA_FREQS_SUBS} from './modules/imgt_aa_freqs_subsampled.nf'
include { COMBINE_IMGT_AA_FREQS_MED } from './modules/combine_IMGT_AA_freqs.nf'
include { COMBINE_IMGT_AA_FREQS_MED as COMBINE_IMGT_AA_FREQS_MED_WL } from './modules/combine_IMGT_AA_freqs.nf'

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

	INITIAL_CLEANUP_SPLIT.out.productive
        .map { sample_name, file_path -> file_path }  // Extract just the files
        .collect()  // Collect all files into a single list
        .map { files -> 
            tuple(files, params.index_columns, "count (templates/reads)","productive") 
        }
        .set { prod_merge_ch }

	// prepare one database of productive sequences
    MERGE_AMINO_VFAM(prod_merge_ch)

	//subsample Vfam, Jgene, CDR3
	MERGE_AMINO_VFAM.out.merged
		.map { merged_file -> tuple(merged_file, params.N, params.M, 3, "productive") } // Example values for N, M, num_index_columns .set { subsample_input_ch }
		.set { subsample_input_ch }
	//sample M times N sequences from each sample
	SUBSAMPLE_FROM_MERGED( subsample_input_ch)



	SUBSAMPLE_FROM_MERGED.out.subsampled_files
		.flatten()
    	.map { subsampled_file -> 
        def sample_name = subsampled_file.baseName.replaceAll(/_subsampled\.tsv$/, '')
        tuple(subsampled_file, sample_name, "productive", "aminoAcid", [ "vFamilyName","jGeneName"])
    } 
    .set { compute_freqs_subsampled_input_ch }
	
	COMPUTE_FREQS_SUBSAMPLED(compute_freqs_subsampled_input_ch)

	COMPUTE_FREQS_SUBSAMPLED.out.subsampled
		.map { big_file, summ_file -> summ_file}
		.collect()
		.map{ files ->
		tuple(files,  "productive_subs_freqs.tsv") }
		.set { combine_freqs_subsampled_input_ch  }	

	COMBINE_FREQS_SUBSAMPLED(combine_freqs_subsampled_input_ch)

	def first_column_with_sample_in_merged = params.index_columns.size()

	MERGE_AMINO_VFAM.out.merged
		.map { merged_file -> tuple(merged_file,  "productive", first_column_with_sample_in_merged) }
		.set { compute_overlaps_input_ch }
		
	COMPUTE_OVERLAPS(compute_overlaps_input_ch) 

	// Coompute overlaps on subsampled
	SUBSAMPLE_FROM_MERGED.out.sample_names
		.collect()
		.combine(SUBSAMPLE_FROM_MERGED.out.subsampled_dir.collect())
		.map { sample_names, subsampled_dirs -> 
			tuple(sample_names, subsampled_dirs,  "productive")  // Assuming single files
		}
	.set { compute_overlaps_subsampled_input_ch }
	COMPUTE_OVERLAPS_SUBSAMPLED(compute_overlaps_subsampled_input_ch  )

	//----aminoacid usage in CDR3
	SUBSAMPLE_FROM_MERGED.out.subsampled_files
		.flatten() 
		.map { subsampled_file -> 
		// println "Processing file: ${subsampled_file}"
        // println "File name: ${subsampled_file.name}"
        // println "Base name: ${subsampled_file.baseName}"
		def sample_name = subsampled_file.baseName.replaceAll(/_subsampled\.tsv$/, '')
		// println "Extracted sample name: ${sample_name}"
		tuple(subsampled_file, sample_name, "productive", params.min_len, params.max_len)
	} 
	.set { compute_imgt_freqs_subsampled_input_ch }	

	COMPUTE_IMGT_AA_FREQS_SUBS(compute_imgt_freqs_subsampled_input_ch)

	// Collect all aa_imgt_freq_full_med files
	COMPUTE_IMGT_AA_FREQS_SUBS.out.aa_imgt_freq_full_med
		.map { sample_name, file -> file }  // Extract just the files
		.collect()  // Collect all files into a single list
		.set { freq_med_files_ch }

	// Combine the frequency median files
	def index_cols = params.imgt_index_columns ?: ["IMGT_position", "AA", "aminoAcid_length"]
	COMBINE_IMGT_AA_FREQS_MED(freq_med_files_ch, index_cols,"_full")

	COMPUTE_IMGT_AA_FREQS_SUBS.out.aa_imgt_freq_WL_med
	.map { sample_name, file -> file }  // Extract just the files
	.collect()  // Collect all files into a single list
	.set { freq_med_WL_files_ch }

	def index_cols_WL = params.imgt_index_columns_WL ?: ["IMGT_position", "AA"]
	COMBINE_IMGT_AA_FREQS_MED_WL(freq_med_WL_files_ch, index_cols_WL, "_WL")

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
