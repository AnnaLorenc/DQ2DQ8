# Pipeline name
Analyse bulk TCR repertoires from Immunoseq

## Workflow description 
### Cleanup: 
- Removal of 2 sequences we identified as contamination
- Only sequences with VFamily are kept
- Split into productive and nonproductive

### Diversities - separately for productive and nonnproductive
- Gini_coefficient
- Shannon
- total_sequences, unique_sequences
- total_count
- top_sequence
- convergence

### Frequencies
- vFamilyName_jFamilyName_cdr3Length - each computed independently
- jGeneName
- vFamilyName_jFamilyName_cdr3Length- each combination

### Subsampling
-from each original sample we sample 10K cells, with probability equal to the frequencies of VfamilyName+JgeneName +CDR3. Each sample is subsampled 25 times

### Frequencies on subsamples
All frequencies are computed also on subsamples. results/freqs_subsampled/productive/XXXX.tsv_summ.tsv contains the original frequencies along with means for subsamples. This is collected across all samples in results/combined_freqs/productive_subs_freqs.tsv.tsv

### Repertoire overlap
For each pair of samples within E and withn N samples, also for E-N samples from the same individual. Computed as number of sequences shared, number of cells with shared sequences, also Jaccard coefficient (numbers of sequences and cells available) *overlaps_raw*.
This was also repeated for subsamples - each combination of 25 subsamples was used.*overlaps_subs*
This is to be grouped by samples,short_id and medians taken:
    sample1	subsample1	sample2	subsample2	unique_overlap	unique1	unique2	jaccard	comparison_type	short_id
    F15625E	F15625E_subsample_1	F16018E	F16018E_subsample_1	12	6874	8684	0.0007719027402547279	within_E	AMIN_VFAM_JGEN  

### Aminoacids inn CDR3 
CDR3 are filtered, to eliminate ones without starting C, ending YV,W, F, too short, too long.
Then each aminoacid is described in terms of AA, IMGT position and CR3_length.
The following computations are perforemd in 2 versions:
-grouping by AA, IMGT position and CR3_length;
-grouping by AA, IMGT position ony (VL)
For each AA-pos-(length), counts and frequencies per sample are computed. This is on the level of sequences (rows) and cells (counts). This is performed for the original sequences as well as for subsamples. Finaly, for subsampled samples, median from subsamples is computed.
The resulting data is collected in the files corresponding to the source sequences and way it was computed:

- `comb_aa_imgt_full_orig_counts.tsv`
- `comb_aa_imgt_full_orig_counts_freq.tsv`
- `comb_aa_imgt_full_orig_rows.tsv`
- `comb_aa_imgt_full_orig_rows_freq.tsv`
- `comb_aa_imgt_full_subs_counts.tsv`
- `comb_aa_imgt_full_subs_counts_freq.tsv`
- `comb_aa_imgt_full_subs_rows.tsv`
- `comb_aa_imgt_full_subs_rows_freq.tsv`
- `comb_aa_imgt_WL_orig_counts.tsv`
- `comb_aa_imgt_WL_orig_counts_freq.tsv`
- `comb_aa_imgt_WL_orig_rows.tsv`
- `comb_aa_imgt_WL_orig_rows_freq.tsv`
- `comb_aa_imgt_WL_subs_counts.tsv`
- `comb_aa_imgt_WL_subs_counts_freq.tsv`
- `comb_aa_imgt_WL_subs_rows.tsv`
- `comb_aa_imgt_WL_subs_rows_freq.tsv`

The main  files to use for modelling are: comb_aa_imgt_full_subs_rows_freq.tsv, comb_aa_imgt_WL_subs_rows_freq.tsv.


## User guide 

## Component tools 

## Additional notes

## Help / FAQ / Troubleshooting

## License(s)

## Acknowledgements/citations/credits
Ania