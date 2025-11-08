#!/bin/bash

pod5=${1}

reference_dir=""

bam_file=""
sample_name=""

# Step 1: Run the basecalling & alignment
nextflow run namle-bioturing/wf-basecalling \
    --basecaller_cfg 'dna_r10.4.1_e8.2_400bps_hac@v5.0.0' \
    --dorado_ext 'pod5' \
    --input ${pod5} \
    --out_dir 
    --ref '/mnt/disk3/namle/data/references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta' \
    --remora_cfg 'dna_r10.4.1_e8.2_400bps_hac@v5.0.0_5mCG_5hmCG@v2' \
    --cuda_device 'cuda:0,2' \
    --sample_name ${sample_name}
    --output_fmt 'bam' \
    --use_parabricks \
    -profile standard \
    -resume


# Step 2: Run the variant calling


# Step 3: Run annotation
