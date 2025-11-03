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

```bash
docker run -ti --rm --gpus '"device=2,3"' -v .:/workspace -v /mnt/disk1/namle:/mnt/disk1/namle --workdir /workspace nvcr.io/nvidia/clara/clara-parabricks:4.5.0-1 bash

/usr/bin/time -v \
pbrun deepvariant \
  --ref /mnt/disk1/namle/data/references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
  --in-bam /workspace/output/SAMPLE.modified.bam \
  --mode ont \
  --out-variants /workspace/out.vcf
```

### Benchmarking

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
    -o bm-filtered-vcfeval
```