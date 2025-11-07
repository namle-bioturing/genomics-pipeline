```bash
python3 process_vep_parallel.py \
    --vcf PAW81754.out.vep.vcf.gz \
    --hgnc-mapping /mnt/nasdev2/namle/references/hgnc.txt \
    --inheritance-genes /mnt/nasdev2/namle/references/inheritance_genes.tsv \
    --omim-info /mnt/nasdev2/namle/references/omim_info.txt \
    --output parquet/vep_result_v3 \
    --processes 20
```