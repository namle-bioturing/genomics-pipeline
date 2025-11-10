#!/bin/bash

SAMPLES_TSV=$1
SAMPLE_TYPE=$2
OUTPUT_DIR=$3

while IFS=$'\t' read -r sample_id fq1 fq2 bench_vcf bench_bed target_bed; do
    [[ "$sample_id" == "sample_id" ]] && continue  # skip header
    [[ "$sample_id" =~ ^#.* ]] && continue  # skip comments

    # Build nextflow command
    NF_CMD="nextflow run parabricks.nf \
        --sample_id \"$sample_id\" \
        --input \"$fq1,$fq2\" \
        --benchmark_vcf \"$bench_vcf\" \
        --benchmark_bed \"$bench_bed\" \
        --sample_type \"$SAMPLE_TYPE\" \
        --output_dir \"$OUTPUT_DIR\""

    # Add target_bed if provided
    if [[ -n "$target_bed" ]]; then
        NF_CMD="$NF_CMD --target_bed \"$target_bed\""
    fi

    NF_CMD="$NF_CMD -resume --profile standard"

    # Execute the command
    eval $NF_CMD
done < "$SAMPLES_TSV"
