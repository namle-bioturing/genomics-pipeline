pbrun germline \
    --ref /workdir/${REFERENCE_FILE} \
    --in-fq /workdir/${INPUT_FASTQ_1} /workdir/${INPUT_FASTQ_2} \
    --knownSites /workdir/${KNOWN_SITES_FILE} \
    --out-bam /outputdir/${OUTPUT_BAM} \
    --out-variants /outputdir/${OUTPUT_VCF} \
    --out-recal-file /outputdir/${OUT_RECAL_FILE}


 pbrun deepvariant \
    --ref /workdir/${REFERENCE_FILE} \
    --in-bam /workdir/${INPUT_BAM} \
    --out-variants /outputdir/${OUTPUT_VCF}

# Example
export HGREF=/mnt/disk1/namle/data/references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta
cd /workspace && /opt/hap.py/bin/hap.py \
    /mnt/disk1/namle/data/benchmarks/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
    ${somevcf} \
    -f /mnt/disk1/namle/data/benchmarks/HG001_GRCh38_1_22_v4.2.1_benchmark.bed \
    --threads 20 \
    --engine vcfeval \
    -o PAW79146-hac-filtered-vcfeval