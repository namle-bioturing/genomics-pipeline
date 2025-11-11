#!/bin/bash -ue
echo "[INTERVAR] Starting InterVar annotation at $(date)"
echo "[INTERVAR] Working directory: $(pwd)"

# Run InterVar
Intervar.py         -c /apps/config.ini         -b hg38         -i PAW81754.out.normalized.vcf.gz         --input_type=VCF         -o PAW81754.intervar         --threads 20         -t /apps/intervardb

echo "[INTERVAR] InterVar annotation completed at $(date)"
