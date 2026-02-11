process COMBINE_DIVERSITY_METRICS {
    tag "Combining ${type} diversity metrics"
    publishDir "${params.outdir}/combined_metrics", mode: 'copy'

    input:
    tuple val(type), val(sample_names), path(csv_files, stageAs: 'input_?.csv')

    output:
    path "combined_${type}_diversity_metrics.csv", emit: combined_diversity

    script:
    def samples_list = sample_names.join(',')
    """
    #!/usr/bin/env python3
    import pandas as pd
    import glob

    # Sample names passed from Nextflow
    sample_names = '${samples_list}'.split(',')
    
    # Get all CSV files in order
    csv_files = sorted(glob.glob("input_*.csv"))
    
    # List to store dataframes
    dfs = []
    
    for sample_name, file in zip(sample_names, csv_files):
        # Read CSV
        df = pd.read_csv(file)
        
        # Add sample_name column at the beginning
        df.insert(0, 'sample_name', sample_name)
        
        dfs.append(df)
    
    # Combine all dataframes
    if dfs:
        combined_df = pd.concat(dfs, ignore_index=True)
        combined_df.to_csv('combined_${type}_diversity_metrics.csv', index=False)
    else:
        pd.DataFrame().to_csv('combined_${type}_diversity_metrics.csv', index=False)
    """
}