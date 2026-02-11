#!/bin/bash

# This script is used to create a samplesheet

# Exit immediately if a command exits with a non-zero status
set -e

# Define the directories containing the samples
SAMPLES_DIRS=("/Users/ania/Documents/DQ2DQ8/data_link/D20210208A_1/samples/" "/Users/ania/Documents/DQ2DQ8/data_link/D20210208A_1/TN_samples/")

# Create the samplesheet
OUTPUT_FILE="samplesheet.csv"
echo -e "Sample\tloc" > "$OUTPUT_FILE"

# Iterate over each directory and process the files
for SAMPLES_DIR in "${SAMPLES_DIRS[@]}"; do
  # Check if the directory exists
  if [ ! -d "$SAMPLES_DIR" ]; then
    echo "Warning: Directory $SAMPLES_DIR does not exist. Skipping." >&2
    continue
  fi

  for file in $(ls "$SAMPLES_DIR" | grep gz); do
    echo -e "$file\t$SAMPLES_DIR" >> "$OUTPUT_FILE"
  done

done

echo "Samplesheet created: $OUTPUT_FILE"