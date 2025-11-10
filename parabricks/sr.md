```bash
# BAM input
nextflow run parabricks_sr.nf \
  --input_type bam \
  --input "sample.bam,sample.bam.bai" \
  --sample_id HG001 \
  --output_dir output

# FASTQ input (single pair)
nextflow run parabricks_sr.nf \
  --input_type fastq \
  --input "sample_R1.fastq.gz,sample_R2.fastq.gz" \
  --sample_id HG001 \
  --output_dir output

# FASTQ input (multiple pairs)
nextflow run parabricks_sr.nf \
  --input_type fastq \
  --input "lane1_R1.fq.gz,lane1_R2.fq.gz,lane2_R1.fq.gz,lane2_R2.fq.gz" \
  --sample_id HG001 \
  --output_dir output
```

```bash
pbrun bamsort \
    --ref /workdir/${REFERENCE_FILE} \
    --in-bam /workdir/${INPUT_BAM} \
    --out-bam /outputdir/${OUTPUT_BAM} \
    --sort-order coordinate

pbrun markdup \
  --ref /workdir/${REFERENCE_FILE} \
  --in-bam /workdir/${INPUT_BAM} \
  --out-bam /outputdir/${OUTPUT_BAM}


pbrun bqsr \
     --ref /workdir/${REFERENCE_FILE} \
     --in-bam /workdir/${INPUT_BAM} \
     --knownSites /workdir/${KNOWN_SITES_FILE} \
     --out-recal-file /outputdir/${INPUT_RECAL_FILE}

  pbrun applybqsr \
    --ref /workdir/${REFERENCE_FILE} \
    --in-bam /workdir/${INPUT_BAM} \
    --in-recal-file /workdir/${INPUT_RECAL_FILE}  \
    --out-bam /outputdir/${OUTPUT_BAM}

 pbrun haplotypecaller \
    --ref /workdir/${REFERENCE_FILE} \
    --in-bam /workdir/${INPUT_BAM_FILE} \
    --in-recal-file /workdir/${INPUT_RECAL_FILE} \
    --out-variants /outputdir/${OUTPUT_VCF_FILE}
```
