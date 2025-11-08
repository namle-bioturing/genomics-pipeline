```bash
nohup rsync -aAXH --numeric-ids -W --inplace --partial \
  --info=progress2 --human-readable \
  /nfsshared/namle/brain-genomic-public/research/sequencing/fastq/novaseq/wgs_pcr_free/30x \
  /mnt/nasdev2/namle/giab/data/wgs_pcr_free/ &

```

```bash
nohup rsync -aAXH --numeric-ids -W --inplace --partial \
  --info=progress2 --human-readable \
  /nfsshared/namle/brain-genomic-public/research/sequencing/fastq/novaseq/wgs_pcr_free/30x \
  /mnt/disk1/namle/data/wgs_pcr_free/ &

```