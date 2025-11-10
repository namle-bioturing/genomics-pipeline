# Genomics Annotation Pipeline

This directory contains the annotation pipeline that combines VEP, InterVar, and custom ACMG annotations.

## Overview

The pipeline performs the following steps:
1. **VCF Normalization** - Normalizes VCF variants
2. **VCF Filtering** - Filters variants using BED file
3. **VEP Annotation** - Adds comprehensive variant effect predictions
4. **InterVar Annotation** - Adds ACMG classification evidence
5. **Combined Annotation** - Merges all annotations into Parquet format

## Files

- `annotation.nf` - Main Nextflow workflow
- `annotate.py` - Python script for combining annotations
- `examples/intervar.example.txt` - Example InterVar output format
- `INTERVAR_VERIFICATION.md` - Technical documentation

## Docker Images

### 1. VEP Container
- Image: `nam.le_bioturing.com/ensembl-vep:release_115.2`
- Contains: Ensembl VEP with plugins (REVEL, AlphaMissense, CADD, SpliceAI)

### 2. InterVar Container
- Image: `nam.le_bioturing.com/intervar:latest`
- Build: `docker build -f docker/intervar.docker -t nam.le_bioturing.com/intervar:latest .`
- Contains: InterVar with ANNOVAR

### 3. Annotate Container
- Image: `nam.le_bioturing.com/annotate:latest`
- Build: `./docker/build_annotate.sh`
- Contains: Python 3.11 with polars, cyvcf2, tabix

## Building Docker Images

### Build InterVar Image
```bash
cd /path/to/genomics-pipeline
docker build -f docker/intervar.docker -t nam.le_bioturing.com/intervar:latest .
docker push nam.le_bioturing.com/intervar:latest
```

### Build Annotate Image
```bash
cd /path/to/genomics-pipeline
./docker/build_annotate.sh
docker push nam.le_bioturing.com/annotate:latest
```

## Reference Files

The pipeline requires the following reference files mounted at `/mnt/nasdev2/namle/references/`:

- `hgnc.txt` - HGNC ID to gene symbol mapping
- `inheritance_genes.tsv` - Gene inheritance patterns
- `omim_info.txt` - OMIM phenotype information
- `GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta` - Reference genome
- `synonyms.txt` - Chromosome synonyms
- `hg38_Twist_Bioscience_for_Illumina_Exome_2_5_Mito.bed` - Target regions

## Running the Pipeline

### Basic Usage
```bash
nextflow run annotation.nf \
  -resume \
  --sample SAMPLE_ID \
  --vcf input.vcf.gz \
  -with-timeline timeline.html \
  -with-report report.html \
  -with-trace trace.txt \
  -profile standard
```

### Parameters
- `--sample` - Sample identifier
- `--vcf` - Input VCF file (bgzipped)
- `--threads` - Number of threads (default: 20)
- `--outdir` - Output directory (default: current directory)

### Example
```bash
nextflow run annotation.nf \
  -resume \
  --sample PAW81754 \
  --vcf PAW81754.vcf.gz \
  --threads 20 \
  -with-timeline timeline.html \
  -with-report report.html \
  -with-trace trace.txt \
  -profile standard
```

## Output Files

### Intermediate Files
- `{sample}.out.normalized.vcf.gz` - Normalized VCF
- `{sample}.out.filtered.vcf.gz` - Filtered VCF (BED regions)
- `{sample}.out.vep.vcf.gz` - VEP-annotated VCF
- `{sample}.intervar.hg38_multianno.txt.intervar` - InterVar classifications

### Final Output
- `{sample}.annotated.parquet` - Combined annotations in Parquet format

## Output Schema

The final Parquet file contains the following columns:

### Variant Information
- `chrom`, `pos`, `ref`, `alt`, `qual`
- `variant_id` - MD5 hash of chrom_pos_ref_alt
- `genes` - Array of all genes associated with variant

### Sample Genotype (FORMAT fields)
- `gt` - Genotype (0=HOM_REF, 1=HET, 2=HOM_ALT)
- `dp` - Total depth
- `gq` - Genotype quality
- `ad` - Alternate allele depth

### VEP Annotations (from PICK=1 transcript)
- `consequence` - Variant consequence
- `feature` - Transcript ID
- `exon`, `intron` - Exon/intron numbers
- `hgvsc`, `hgvsp` - HGVS nomenclature
- `rsid` - dbSNP ID
- `hgnc_id` - HGNC ID
- `mane` - MANE transcript

### Prediction Scores
- `sift`, `polyPhen` - Deleteriousness predictions
- `revel` - REVEL score
- `cadd_phred`, `cadd_raw` - CADD scores
- `spliceai_pred_*` - SpliceAI predictions
- `clinsig` - ClinVar significance
- `am_class`, `am_pathogenicity` - AlphaMissense predictions

### Population Frequencies
- `af` - Global allele frequency
- `afr_af`, `amr_af`, `eas_af`, `eur_af`, `sas_af` - 1000 Genomes
- `gnomade_*` - gnomAD Exomes frequencies
- `gnomadg_*` - gnomAD Genomes frequencies
- `max_af`, `max_af_pops` - Maximum frequency

### OMIM Annotations
- `omim_inheritance` - Inheritance patterns
- `omim_phenotype` - Associated phenotypes

### ACMG Classifications (from InterVar)
- `ACMG_classification` - Most significant classification
- `ACMG_evidences` - Evidence codes (e.g., "PVS1, PS1")
- `ACMG_data` - Array of per-gene ACMG data

## ACMG Classification Significance Order
1. Pathogenic
2. Likely pathogenic
3. Uncertain significance
4. Likely benign
5. Benign

## Querying Parquet Output

### Using Python (Polars)
```python
import polars as pl

# Load annotations
df = pl.read_parquet("SAMPLE.annotated.parquet")

# Filter pathogenic variants
pathogenic = df.filter(
    pl.col("ACMG_classification").is_in(["Pathogenic", "Likely pathogenic"])
)

# Get variants in specific gene
brca1_variants = df.filter(
    pl.col("genes").list.contains("BRCA1")
)

# High impact variants
high_impact = df.filter(
    (pl.col("consequence").str.contains("stop_gained|frameshift")) |
    (pl.col("cadd_phred") > 20)
)
```

### Using DuckDB
```python
import duckdb

con = duckdb.connect()
result = con.execute("""
    SELECT chrom, pos, ref, alt, genes, ACMG_classification, ACMG_evidences
    FROM 'SAMPLE.annotated.parquet'
    WHERE ACMG_classification IN ('Pathogenic', 'Likely pathogenic')
    ORDER BY chrom, pos
""").df()
```

## Performance

- **VEP**: ~5-10 minutes for WES (~100K variants)
- **InterVar**: ~10-15 minutes for WES
- **Annotate**: ~2-5 minutes for WES
- **Total**: ~20-30 minutes for typical WES sample

Memory requirements:
- VEP: 20 GB
- InterVar: 20 GB
- Annotate: 32 GB

## Troubleshooting

### Docker Issues
If you see "Command 'ps' required by nextflow":
```bash
# Rebuild with procps included
docker build -f docker/intervar.docker -t nam.le_bioturing.com/intervar:latest .
```

### InterVar Column Mismatch
If InterVar output format changes, update column name in `annotate.py`:
```python
columns=["#Chr", "Start", "Ref", "Alt", "Ref.Gene", "InterVar: InterVar and Evidence"]
```

### Memory Issues
Increase memory allocation in `annotation.nf`:
```groovy
memory '64 GB'  // Increase as needed
```

## References

- **VEP**: https://www.ensembl.org/info/docs/tools/vep/
- **InterVar**: https://github.com/WGLab/InterVar
- **ACMG Guidelines**: https://www.acmg.net/
- **Polars**: https://pola.rs/

## Citation

If you use this pipeline, please cite:
- Ensembl VEP
- InterVar
- ANNOVAR
- Relevant VEP plugins (REVEL, AlphaMissense, CADD, SpliceAI)


```bash
# Start a tmux session
tmux new-session -s namle-annotation

# Attach later
tmux attach -t namle-annotation

# Remove 
tmux kill-session -t namle-annotation

# View logs
tmux capture-pane -pt namle-annotation -S -
```
