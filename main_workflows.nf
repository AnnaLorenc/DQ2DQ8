#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Parameters
params.input = './data/collated_info.csv'
params.outdir = './results'
params.seqs_to_remove = './data/D20210208D_1-overrepresented-sequences.txt'
params.N = 10000 
params.M = 25
params.index_columns = ["aminoAcid", "vFamilyName", "jGeneName"]

// Include subworkflows
include { DATA_PREPROCESSING } from './workflows/data_preprocessing.nf'
include { DIVERSITY_ANALYSIS } from './workflows/diversity_analysis.nf'
include { FREQUENCY_ANALYSIS } from './workflows/frequency_analysis.nf'
include { MERGING_AND_SUBSAMPLING } from './workflows/merging_and_subsampling.nf'
include { SUBSAMPLED_ANALYSIS } from './workflows/subsampled_freqs.nf'
include { OVERLAP_ANALYSIS } from './workflows/overlap_analysis.nf'

// Print header
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
    // Create input channels
    sample_filepath_ch = Channel.fromPath(params.input)
        .splitCsv(header: true, sep: ',')
        .map { row -> tuple(row.sample_short, "${row.LOC}/${row.SAMPLE}", params.seqs_to_remove, params.outdir) }

    // Run subworkflows
    DATA_PREPROCESSING(sample_filepath_ch)
    
    DIVERSITY_ANALYSIS(
        DATA_PREPROCESSING.out.productive,
        DATA_PREPROCESSING.out.nonproductive
    )
    
    FREQUENCY_ANALYSIS(
        DATA_PREPROCESSING.out.productive,
        DATA_PREPROCESSING.out.nonproductive
    )
    
    MERGING_AND_SUBSAMPLING(
        DATA_PREPROCESSING.out.productive,
        params.index_columns,
        params.N,
        params.M
    )
    
    SUBSAMPLED_ANALYSIS(
        MERGING_AND_SUBSAMPLING.out.subsampled
    )
    
    // Prepare data for overlap analysis
    subsampled_data_ch = MERGING_AND_SUBSAMPLING.out.sample_names
        .combine(MERGING_AND_SUBSAMPLING.out.subsampled_dir)
        .map { sample_names_file, subsampled_dir -> 
            tuple(sample_names_file, subsampled_dir, params.index_columns)
        }
    
    OVERLAP_ANALYSIS(
        MERGING_AND_SUBSAMPLING.out.merged,
        subsampled_data_ch,
        params.index_columns
    )
}