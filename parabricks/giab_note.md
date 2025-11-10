```bash
/bin/bash /mnt/disk1/namle/run/giab/run_samples.sh /mnt/disk1/namle/run/giab/wgs_samples_test.tsv wgs /mnt/disk1/namle/run/giab/output/wgs_pcr_free

/bin/bash /mnt/disk1/namle/run/giab/run_samples.sh /mnt/disk1/namle/run/giab/wgs_samples.tsv wgs /mnt/disk1/namle/run/giab/output/wgs_pcr_free
```

```bash
# Start a tmux session
tmux new-session -s giab-wes-idt

# Attach later
tmux attach -t giab-wes-idt

# Remove 
tmux kill-session -t giab-wes-idt

# View logs
tmux capture-pane -pt giab-wes-idt -S - | less
```

```bash
/bin/bash /mnt/disk1/namle/run/giab/run_samples.sh /mnt/disk1/namle/run/giab/wes_samples_idt_50x.tsv wes /mnt/disk1/namle/run/giab/output/wes_idt_50x

/bin/bash /mnt/disk1/namle/run/giab/run_samples.sh /mnt/disk1/namle/run/giab/wes_samples_idt_75x.tsv wes /mnt/disk1/namle/run/giab/output/wes_idt_75x

/bin/bash /mnt/disk1/namle/run/giab/run_samples.sh /mnt/disk1/namle/run/giab/wes_samples_idt_100x.tsv wes /mnt/disk1/namle/run/giab/output/wes_idt_100x
```

```bash
# Start a tmux session
tmux new-session -s deepvariant-g400

# Attach later
tmux attach -t deepvariant-g400

# Remove 
tmux kill-session -t deepvariant-g400

# View logs
tmux capture-pane -pt deepvariant-g400 -S -
```

```bash
# Run bam
bam=/mnt/disk1/namle/data/complete_g400/HG002.complete_g400.V350151728.grch38.bam
bai=/mnt/disk1/namle/data/complete_g400/HG002.complete_g400.V350151728.grch38.bam.bai
nextflow run main_sr.nf \
    -resume \
    --sample_id HG002 \
    --input "$bam,$bai" \
    --benchmark_vcf "/mnt/disk1/namle/data/benchmarks/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" \
    --benchmark_bed "/mnt/disk1/namle/data/benchmarks/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed" \
    --sample_type "wgs" \
    --input_type "bam" \
    --gg_deepvariant_model "/mnt/disk1/namle/data/complete_g400/weights-60-0.993753.ckpt" \
    --output_dir "/mnt/disk1/namle/run/giab/output/complete_g400"
```

```bash
# Run bam
bam=/mnt/disk1/namle/data/tests/HG002.bam
bai=/mnt/disk1/namle/data/tests/HG002.bam.bai
nextflow run main_sr.nf \
    -resume \
    --sample_id HG002 \
    --input "$bam,$bai" \
    --benchmark_vcf "/mnt/disk1/namle/data/benchmarks/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" \
    --benchmark_bed "/mnt/disk1/namle/data/benchmarks/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed" \
    --sample_type "wgs" \
    --input_type "bam" \
    --output_dir "/mnt/disk1/namle/run/giab/output/tests"
```