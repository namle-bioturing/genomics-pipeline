### Check nvidia-docker2 installation

```bash
docker run --rm --gpus 1,2 nvidia/cuda:12.9-base-ubuntu22.04 nvidia-smi
```

### Deploy NVIDIA Parabricks

```bash
docker pull nvcr.io/nvidia/clara/clara-parabricks:4.5.0-1
```

### Getting The Sample Data

```bash
wget -O parabricks_sample.tar.gz "https://s3.amazonaws.com/parabricks.sample/parabricks_sample.tar.gz"
```


```bash
docker run -ti --rm --gpus '"device=0,1"' -v .:/workspace --volume $(pwd):/workdir --volume $(pwd):/outputdir nvcr.io/nvidia/clara/clara-parabricks:4.5.0-1 bash
```

```bash
pbrun minimap2 \
    --ref /workdir/epi2me/wf-basecalling-demo/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
    --in-fq /workdir/parabricks/out.fq \
    --out-bam /outputdir/out.bam
```

## Variant calling

### Variant calling PAW79146

```bash
docker run -ti --rm --gpus '"device=2,3"' -v .:/workspace -v /mnt/disk1/namle:/mnt/disk1/namle --workdir /workspace nvcr.io/nvidia/clara/clara-parabricks:4.5.0-1 bash

/usr/bin/time -v \
pbrun deepvariant \
  --ref /mnt/disk1/namle/data/references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
  --in-bam /workspace/output/SAMPLE.modified.bam \
  --mode ont \
  --out-variants /workspace/out.vcf
```

### Variant calling PAW81754

```bash
samtools addreplacerg -w  -r  "ID:150d2e3f-5093-4f3e-bf40-714b9c053121_dna_r10.4.1_e8.2_400bps_hac@v5.0.0\tSM:SAMPLE\tPU:PAW81754\tPM:PCAMP156     DT:2024-09-11T11:29:15.831+00:00\tPL:ONT\tDS:basecall_model=dna_r10.4.1_e8.2_400bps_hac@v5.0.0 modbase_models=dna_r10.4.1_e8.2_400bps_hac@v5.0.0_5mCG_5hmCG@v2 runid=150d2e3f-5093-4f3e-bf40-714b9c053121\tLB:hg001" -@ 50 -o SAMPLE.modified.bam SAMPLE.pass.bam

docker run -ti --rm --gpus '"device=2,3"' -v .:/workspace -v /mnt/disk3/namle:/mnt/disk3/namle --workdir /workspace nvcr.io/nvidia/clara/clara-parabricks:4.5.0-1 bash

/usr/bin/time -v \
pbrun deepvariant \
  --ref /mnt/disk3/namle/data/references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
  --in-bam /workspace/output/SAMPLE.modified.bam \
  --mode ont \
  --out-variants /workspace/out.vcf

bcftools view -f "PASS,LowQual" out.vcf -Oz -o filtered.vcf.gz
```


## Benchmarking

### Benchmarking for PAW79146

```bash
docker run -itv .:/workspace -v /mnt/disk1/namle:/mnt/disk1/namle hap.py:latest bash

export HGREF=/mnt/disk1/namle/data/references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta
cd /workspace && /opt/hap.py/bin/hap.py \
    /mnt/disk1/namle/data/benchmarks/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
    /mnt/disk1/namle/run/PAW79146/parabricks/wf-human-variation/filtered.vcf.gz \
    -f /mnt/disk1/namle/data/benchmarks/HG001_GRCh38_1_22_v4.2.1_benchmark.bed \
    --threads 80 \
    --engine vcfeval \
    --pass-only \
    -o PAW79146-hac-filtered-vcfeval
```

### Benchmarking for PAW81754

```bash
docker run -itv .:/workspace -v /mnt/disk3/namle:/mnt/disk3/namle hap.py:latest bash

export HGREF=/mnt/disk3/namle/data/references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta
cd /workspace && /opt/hap.py/bin/hap.py \
    /mnt/disk3/namle/data/benchmarks/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
    /mnt/disk3/namle/run/PAW81754/parabricks/wf-human-variation/filtered.vcf.gz \
    -f /mnt/disk3/namle/data/benchmarks/HG001_GRCh38_1_22_v4.2.1_benchmark.bed \
    --threads 80 \
    --engine vcfeval \
    -o PAW81754-hac-filtered-vcfeval
```
