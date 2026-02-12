import polars as pl
import numpy as np
import gzip
import os

def subsample_merged_sequences(input_file, output_dir, N, M, num_index_columns):
    """
    Subsample merged sequences based on weights from the input file.

    Parameters:
        input_file (str): Path to the input merged file (tsv.gz).
        output_dir (str): Directory to save the output files.
        N (int): Number of indexes to sample per column.
        M (int): Number of repetitions for subsampling.
        num_index_columns (int): Number of columns to treat as indexes.
    """
    # Read the input file
    with gzip.open(input_file, "rt") as f:
        df = pl.read_csv(f, separator="\t")

    # Treat the specified number of columns as indexes
    index_columns = df.columns[:num_index_columns]
    sample_columns = df.columns[num_index_columns:]

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    for sample in sample_columns[0:1]:
        # Filter rows where the sample column > 0
        filtered_df = df.filter(pl.col(sample) > 0)

        # Extract weights and indexes
        weights = filtered_df[sample].to_numpy()
        indexes = filtered_df[index_columns].to_numpy()

        # Prepare the output DataFrame
        output_df = filtered_df.select(index_columns + [pl.col(sample).alias(f"{sample}_orig")])

        # Perform M subsampling iterations
        for i in range(M):
            sampled_indexes = np.random.choice(
                len(indexes), size=N, replace=True, p=weights / weights.sum()
            )

            # Create a column initialized with zeros
            sampled_column = np.zeros(len(filtered_df), dtype=int)

            # Increment counts for sampled rows
            for idx in sampled_indexes:
                sampled_column[idx] += 1

            # Add the sampled column to the output DataFrame
            output_df = output_df.with_columns(
                pl.Series(f"{sample}_subsample_{i + 1}", sampled_column)
            )

        # Write the output to a gzipped TSV file
        output_file = os.path.join(output_dir, f"{sample}_subsampled.tsv.gz")
        output_df.write_csv(output_file, separator="\t", compression="gzip")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Subsample merged sequences based on weights.")
    parser.add_argument("--input_file", required=True, help="Path to the input merged file (tsv.gz).")
    parser.add_argument("--output_dir", required=True, help="Directory to save the output files.")
    parser.add_argument("--N", type=int, required=True, help="Number of indexes to sample per column.")
    parser.add_argument("--M", type=int, required=True, help="Number of repetitions for subsampling.")
    parser.add_argument("--num_index_columns", type=int, required=True, help="Number of columns to treat as indexes.")

    args = parser.parse_args()

    subsample_merged_sequences(args.input_file, args.output_dir, args.N, args.M, args.num_index_columns)