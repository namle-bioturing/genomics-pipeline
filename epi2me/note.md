### Test

#### Moving
```bash
nohup rsync -ah --info=progress2 --stats \
  /nfsshared/personal/nam.le_bioturing.com/data/giab/epi2me/HG001/PAW79146 \
  /mnt/disk1/namle/data/ &
```

#### Containers
```bash
docker run -tiv .:/workspace -v /mnt/disk1/namle:/mnt/disk1/namle --user $(id -u):$(id -g) --group-add 100 ontresearch/dorado:shae423e761540b9d08b526a1eb32faf498f32e8f22 bash
```

#### Base calling
```bash
wget https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-basecalling/wf-basecalling-demo.tar.gz
tar -xzvf wf-basecalling-demo.tar.gz
```

```bash
nextflow run epi2me-labs/wf-basecalling \
    --basecaller_cfg 'dna_r10.4.1_e8.2_400bps_hac@v5.0.0' \
    --dorado_ext 'pod5' \
    --input 'wf-basecalling-demo/input' \
    --ref 'wf-basecalling-demo/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta' \
    --remora_cfg 'dna_r10.4.1_e8.2_400bps_hac@v5.0.0_5mCG_5hmCG@v2' \
    --cuda_device 'cuda:0,4' \
    --output_fmt 'bam' \
    -profile standard

nextflow pull namle-bioturing/wf-basecalling && \
nextflow run namle-bioturing/wf-basecalling \
    --basecaller_cfg 'dna_r10.4.1_e8.2_400bps_hac@v5.0.0' \
    --dorado_ext 'pod5' \
    --input 'wf-basecalling-demo/input' \
    --ref 'wf-basecalling-demo/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta' \
    --remora_cfg 'dna_r10.4.1_e8.2_400bps_hac@v5.0.0_5mCG_5hmCG@v2' \
    --cuda_device 'cuda:0,4' \
    --output_fmt 'bam' \
    --use_parabricks \
    -profile standard
```

#### Human variation
```bash
# For testing the base calling result
nextflow run epi2me-labs/wf-human-variation \
    --bam 'basecalling-output/SAMPLE.pass.bam' \
    --ref 'wf-basecalling-demo/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta' \
    --sample_name 'SAMPLE' \
    --snp \
    --annotation false \
    -profile standard

# For testing the human variation workflow
nextflow run epi2me-labs/wf-human-variation -r master \
    --bam 'wf-human-variation-demo/demo.bam' \
    --ref 'wf-human-variation-demo/demo.fasta' \
    --bed 'wf-human-variation-demo/demo.bed' \
    --sample_name 'DEMO' \
    --snp \
    --annotation false \
    -profile standard
```

### Run for HG001 - 1st FlowCell


#### Base calling
```bash
nextflow run epi2me-labs/wf-basecalling \
    -bg \
    -resume \
    --basecaller_args '--sample-sheet /mnt/disk1/namle/data/PAW79146/sample_sheet_PAW79146_20240911_1129_c6420e1f.csv' \
    --basecaller_cfg 'dna_r10.4.1_e8.2_400bps_hac@v5.0.0' \
    --dorado_ext 'pod5' \
    --input '/mnt/disk1/namle/data/PAW79146/pod5' \
    --ref '/mnt/disk1/namle/data/references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta' \
    --remora_cfg 'dna_r10.4.1_e8.2_400bps_hac@v5.0.0_5mCG_5hmCG@v2' \
    --cuda_device 'cuda:0,4' \
    --output_fmt 'bam' \
    --ubam_map_threads 25 \
    --ubam_sort_threads 10 \
    --ubam_bam2fq_threads 10 \
    -profile standard

nextflow pull namle-bioturing/wf-basecalling && \
nextflow run namle-bioturing/wf-basecalling \
    -bg \
    --basecaller_cfg 'dna_r10.4.1_e8.2_400bps_hac@v5.0.0' \
    --dorado_ext 'pod5' \
    --input '/mnt/disk1/namle/data/PAW79146/pod5' \
    --ref '/mnt/disk1/namle/data/references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta' \
    --remora_cfg 'dna_r10.4.1_e8.2_400bps_hac@v5.0.0_5mCG_5hmCG@v2' \
    --cuda_device 'cuda:0,2' \
    --output_fmt 'bam' \
    --use_parabricks \
    -profile standard \
    -resume
```


#### Human variation

```bash
nextflow run epi2me-labs/wf-human-variation \
    --bam 'wf-basecalling-output/SAMPLE.pass.bam' \
    --ref '/mnt/disk1/namle/data/references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta' \
    --sample_name 'PAW79146' \
    --snp \
    --annotation false \
    -profile standard

nextflow run epi2me-labs/wf-human-variation \
    --bam 'basecalling-output/SAMPLE.pass.bam' \
    --ref '/mnt/disk1/namle/data/references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta' \
    --sample_name 'PAW79146' \
    --snp \
    --annotation false \
    -profile standard
```

### Run for HG001 - 2nd FlowCell


#### Base calling
```bash
nextflow run epi2me-labs/wf-basecalling \
    -bg \
    -resume \
    --basecaller_cfg 'dna_r10.4.1_e8.2_400bps_hac@v5.0.0' \
    --dorado_ext 'pod5' \
    --input '/mnt/disk3/namle/data/PAW81754/pod5' \
    --ref '/mnt/disk3/namle/data/references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta' \
    --remora_cfg 'dna_r10.4.1_e8.2_400bps_hac@v5.0.0_5mCG_5hmCG@v2' \
    --cuda_device 'cuda:0,2' \
    --output_fmt 'bam' \
    --ubam_map_threads 25 \
    --ubam_sort_threads 10 \
    --ubam_bam2fq_threads 10 \
    -profile standard

nextflow pull namle-bioturing/wf-basecalling && \
nextflow run namle-bioturing/wf-basecalling \
    -bg \
    --basecaller_cfg 'dna_r10.4.1_e8.2_400bps_hac@v5.0.0' \
    --dorado_ext 'pod5' \
    --input '/mnt/disk3/namle/data/PAW81754/pod5' \
    --ref '/mnt/disk3/namle/data/references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta' \
    --remora_cfg 'dna_r10.4.1_e8.2_400bps_hac@v5.0.0_5mCG_5hmCG@v2' \
    --cuda_device 'cuda:0,2' \
    --output_fmt 'bam' \
    --use_parabricks \
    -profile standard \
    -resume
```


#### Human variation

```bash
nextflow run epi2me-labs/wf-human-variation \
    --bam 'basecalling-output/SAMPLE.pass.bam' \
    --ref '/mnt/disk3/namle/data/references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta' \
    --sample_name 'SAMPLE' \
    --snp \
    -profile standard

nextflow run epi2me-labs/wf-human-variation \
    --bam 'basecalling-output/SAMPLE.pass.bam' \
    --ref '/mnt/disk3/namle/data/references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta' \
    --sample_name 'SAMPLE' \
    --snp \
    --annotation false \
    -profile standard

```

### Benchmarking

```bash
docker run -itv .:/workspace -v /mnt/disk1/namle:/mnt/disk1/namle hap.py:latest bash

export HGREF=/mnt/disk1/namle/data/references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta
cd /workspace && /opt/hap.py/bin/hap.py \
    /mnt/disk1/namle/data/benchmarks/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
    /mnt/disk1/namle/run/PAW79146/epi2me/wf-human-variation/PAW79146.wf_snp.vcf.gz \
    -f /mnt/disk1/namle/data/benchmarks/HG001_GRCh38_1_22_v4.2.1_benchmark.bed \
    --threads 80 \
    -o test
```