# Setup Summary - Annotation Pipeline

## ✅ Files Created

### 1. Docker Configuration
- **File**: [docker/annotate.docker](../docker/annotate.docker)
- **Purpose**: Dockerfile for running annotate.py script
- **Contents**: Python 3.11 + polars + cyvcf2 + tabix + procps

### 2. Build Script
- **File**: [docker/build_annotate.sh](../docker/build_annotate.sh)
- **Purpose**: Automated build script for Docker image
- **Usage**: `./docker/build_annotate.sh [tag]`

### 3. Updated Workflow
- **File**: [annotation.nf](annotation.nf)
- **Changes**:
  - Added annotate container configuration (lines 33-40)
  - Added ANNOTATE_VCF process (lines 248-285)
  - Updated workflow to call ANNOTATE_VCF (lines 59-65)

### 4. Documentation
- **File**: [README.md](README.md)
- **Contents**: Complete pipeline documentation

## 🚀 Quick Start

### Step 1: Build Docker Image
```bash
cd /home/bioturing/Documents/GitHub/bioturing/genomics-pipeline
./docker/build_annotate.sh
```

### Step 2: Push to Registry
```bash
docker push nam.le_bioturing.com/annotate:latest
```

### Step 3: Run Pipeline
```bash
cd annotations
nextflow run annotation.nf \
  -resume \
  --sample PAW81754 \
  --vcf PAW81754.vcf.gz \
  -with-timeline timeline.html \
  -with-report report.html \
  -with-trace trace.txt \
  -profile standard
```

## 📊 Pipeline Flow

```
Input VCF
    ↓
NORMALIZE_VCF (bcftools norm)
    ↓
    ├─→ RUN_VEP ────────────────┐
    │                           │
    └─→ RUN_INTERVAR ──────────┤
                                ↓
                         ANNOTATE_VCF
                                ↓
                    {sample}.annotated.parquet
```

## 🔧 Reference Files Required

All mounted from `/mnt/nasdev2/namle/references/`:
- ✅ `hgnc.txt`
- ✅ `inheritance_genes.tsv`
- ✅ `omim_info.txt`
- ✅ VEP cache and plugins
- ✅ ANNOVAR databases

## 📦 Output

**Final File**: `{sample}.annotated.parquet`

**Contains**:
- VEP annotations (consequence, predictions, frequencies)
- OMIM annotations (inheritance, phenotypes)
- ACMG classifications (InterVar-based)
- All genes per variant
- Sample genotype information

## ⚠️ Important Notes

1. **VEP Output Required**: ANNOTATE_VCF depends on RUN_VEP completing successfully
2. **InterVar Format**: Must match column structure in `intervar.example.txt`
3. **Memory**: ANNOTATE_VCF requires 32GB memory for typical WES samples
4. **Tabix Index**: VCF files must be bgzipped and indexed (.tbi)

## 🐛 Troubleshooting

### Missing `ps` command
Already fixed in `docker/intervar.docker` with `procps` package.

### Docker build fails
Make sure you're in the repository root when building:
```bash
cd /home/bioturing/Documents/GitHub/bioturing/genomics-pipeline
./docker/build_annotate.sh
```

### InterVar column mismatch
Column name is now correct in `annotate.py` (line 184):
```python
"InterVar: InterVar and Evidence"
```

### Process hangs or crashes
Check Nextflow work directory:
```bash
tail -f work/*/*/.command.log
```

## 📝 Next Steps

1. ✅ Build annotate Docker image
2. ✅ Push to registry
3. ✅ Test with sample data
4. ✅ Verify Parquet output schema
5. ✅ Query results with Polars/DuckDB

## 🔍 Verification

After pipeline completes, verify output:
```python
import polars as pl

df = pl.read_parquet("PAW81754.annotated.parquet")
print(df.shape)  # Should show (n_variants, ~90 columns)
print(df.columns)  # Check all expected columns present
print(df.select(["chrom", "pos", "genes", "ACMG_classification"]).head())
```

## 📚 Documentation

- [Full README](README.md) - Complete documentation
- [InterVar Verification](INTERVAR_VERIFICATION.md) - Technical details
- [annotate.py](annotate.py) - Source code with docstrings
