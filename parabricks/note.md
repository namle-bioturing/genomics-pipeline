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
docker run -ti --rm --gpus '"device=0,1"' --volume $(pwd):/workdir --volume $(pwd):/outputdir nvcr.io/nvidia/clara/clara-parabricks:4.5.0-1 bash
```

```bash
pbrun minimap2 \
    --ref /workdir/epi2me/wf-basecalling-demo/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
    --in-fq /workdir/parabricks/out.fq \
    --out-bam /outputdir/out.bam
```