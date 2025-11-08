

# Run VEP

## Build docker
```bash
docker build -t nam.le_bioturing.com/ensembl-vep:release_115.2 .
```

## Run docker
```bash
docker run -u $(id -u):$(id -g) -it --rm --name test_vep \
    -v /mnt/nasdev2/namle/references:/references \
    -v /mnt/nasdev2/namle/run:/workspace \
    -v /mnt/nasdev2/namle/.vep:/opt/vep/.vep \
    --workdir /workspace \
    nam.le_bioturing.com/ensembl-vep:release_115.2 \
    bash
    /workspace/run_vep.sh PAW81754 PAW81754.vcf.gz
```