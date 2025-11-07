# VEP VCF Annotation Processor

High-performance Python script for processing VEP-annotated VCF files using cyvcf2, polars, and multiprocessing.

**Script:**
- `process_vep_parallel.py`: Multi-core parallel version for large datasets (processes by chromosome)

## Features

- **Single-Sample, Single-ALT VCF**: Designed for VCF files with one sample and one ALT allele per variant
- **Standard Chromosome Filtering**: Only processes standard chromosomes (1-22, X, Y, MT or chr1-chr22, chrX, chrY, chrM)
- **Parallel Processing**: Multi-core chromosome-level parallelization for faster processing of large VCFs
- **Intelligent Gene Extraction**: Prioritizes HGNC_ID over SYMBOL with automatic deduplication
- **Selective Field Extraction**: Captures all VCF standard fields, INFO fields, FORMAT fields, and 45 selected CSQ fields
- **Quality Filtering**: Automatically excludes CSQ records where BAM_EDIT = FAILED
- **High Performance**: Uses cyvcf2 for fast VCF parsing and Polars for efficient DataFrame operations
- **Optimized Output**: Writes to Parquet format with ZSTD compression and statistics enabled
- **Progress Tracking**: Logs progress every 10,000 variants (with chromosome prefixes in parallel mode)
- **Robust Error Handling**: Gracefully handles malformed CSQ entries and missing values

## Installation

```bash
pip install -r requirements.txt
```

## Requirements

- Python 3.8+
- cyvcf2 >= 0.30.0
- polars >= 0.19.0

## Usage

```bash
python process_vep_parallel.py \
  --vcf annotated.vcf.gz \
  --hgnc-mapping hgnc_map.tsv \
  --output results.parquet \
  --processes 8
```

**Requirements:**
- VCF file must be bgzipped (`.vcf.gz`)
- Tabix index must exist (`.vcf.gz.tbi`)
- Create index with: `tabix -p vcf annotated.vcf.gz`

### Arguments

- `--vcf`: Path to VEP-annotated VCF file (must be .vcf.gz with .tbi index)
- `--hgnc-mapping`: Path to TSV file with columns [hgnc_id, hgnc_symbol]
- `--output`: Path for output Parquet file
- `--processes`: Number of parallel processes (default: number of CPU cores)

### HGNC Mapping File Format

The HGNC mapping file should be a TSV with two columns:

```
hgnc_id	hgnc_symbol
1	A1BG
2	A2M
5	NAT1
```

## Gene Extraction Logic

For each variant, the script extracts variant-gene pairs using this priority:

1. **Quality Filter**: Skip CSQ records where BAM_EDIT = "FAILED"
2. **HGNC Priority**: If HGNC_ID exists in CSQ record, map it to HGNC symbol using provided mapping
3. **Symbol Fallback**: If no HGNC_ID, use SYMBOL field from CSQ
4. **Deduplication**: Track genes per variant - only output each gene once per variant
5. **Exclusion**: Skip records with no gene OR if gene was already seen for this variant

### Example

For a variant with 5 CSQ records:
- Record 1: BAM_EDIT="FAILED" → Skip (quality filter)
- Record 2: HGNC_ID=1 (maps to "GENE_E") → Extract: Variant_A - GENE_E
- Record 3: No HGNC_ID, SYMBOL="GENE_F" → Extract: Variant_A - GENE_F
- Record 4: No HGNC_ID, No SYMBOL → Skip (no gene)
- Record 5: No HGNC_ID, SYMBOL="GENE_E" → Skip (already extracted)

Result: 2 records output

## Parallel Processing Architecture

The parallel version (`process_vep_parallel.py`) distributes work by chromosome:

1. **Main Process**:
   - Parses VCF header and extracts chromosome list (filtered to standard chromosomes only)
   - Loads HGNC mapping (shared across workers)
   - Spawns worker processes using `multiprocessing.Pool`

2. **Worker Processes**:
   - Each worker processes one chromosome independently
   - Opens its own VCF handle (cyvcf2 is not thread-safe)
   - Uses region queries to iterate only assigned chromosome
   - Returns list of record dictionaries (avoids DataFrame pickling overhead)

3. **Result Combination**:
   - Main process collects all worker results
   - Flattens into single list
   - Converts to Polars DataFrame once
   - Writes single Parquet file

**Performance Benefits**:
- Linear scaling with number of CPU cores (up to number of chromosomes)
- Efficient memory usage (workers return lists, not DataFrames)
- Progress logging from each worker with chromosome prefixes

## Output Format

The output Parquet file contains one row per variant-gene pair with:

**Empty Field Handling**: All empty or missing fields are stored as `null` values in Parquet for optimal type preservation and query performance

**Benefits of null handling:**
- Numeric fields (qual, dp, gq, ad, allele frequencies) maintain their native types (Int64, Float64)
- Easy filtering with `.is_null()` and `.is_not_null()` methods
- Better compression and smaller file sizes
- Direct numeric operations without type casting

### Standard VCF Columns
- `chrom`, `pos`, `ref`, `alt`, `qual`
- **Note**: Each variant contains only a single ALT allele (no multi-allelic variants)

### Computed Columns
- `variant_id`: MD5 hash of "{CHROM}_{POS}_{REF}_{ALT}" (32 characters)
- `pair_id`: MD5 hash of "{CHROM}_{POS}_{REF}_{ALT}_{GENE}" (32 characters)
- `gene`: Extracted gene symbol

### FORMAT Fields
- `gt`: Genotype
- `dp`: Read depth
- `gq`: Genotype quality
- `ad`: Allelic depth

### CSQ Fields
Selected CSQ fields from VEP annotation (see below for full list)

### CSQ Fields Extracted

The script extracts only the following 45 CSQ fields (to reduce output size):

**Functional Annotations:**
- Consequence, Feature, EXON, INTRON, HGVSc, HGVSp, Existing_variation, HGNC_ID, MANE, BAM_EDIT

**Prediction Scores:**
- SIFT, PolyPhen, REVEL, CLIN_SIG, am_class, am_pathogenicity

**Population Frequencies:**
- AF, AFR_AF, AMR_AF, EAS_AF, EUR_AF, SAS_AF
- gnomADe_AF, gnomADe_AFR_AF, gnomADe_AMR_AF, gnomADe_ASJ_AF, gnomADe_EAS_AF, gnomADe_FIN_AF, gnomADe_MID_AF, gnomADe_NFE_AF, gnomADe_REMAINING_AF, gnomADe_SAS_AF
- gnomADg_AF, gnomADg_AFR_AF, gnomADg_AMI_AF, gnomADg_AMR_AF, gnomADg_ASJ_AF, gnomADg_EAS_AF, gnomADg_FIN_AF, gnomADg_MID_AF, gnomADg_NFE_AF, gnomADg_REMAINING_AF, gnomADg_SAS_AF
- MAX_AF, MAX_AF_POPS

### Customizing CSQ Field Names

The script uses a field mapping dictionary (`CSQ_FIELD_MAPPING`) to control which CSQ fields are extracted and what they're named in the output.

**To customize output field names:**
1. Edit the `CSQ_FIELD_MAPPING` dictionary at the top of the script
2. Change the value (output name) for any field you want to rename
3. Example:

```python
CSQ_FIELD_MAPPING = {
    'Consequence': 'variant_consequence',  # Renamed
    'HGVSc': 'hgvs_coding',                # Renamed
    'SIFT': 'SIFT',                        # Unchanged
    # ... etc
}
```

**To add/remove fields:**
- Add new field: Add a new key-value pair to `CSQ_FIELD_MAPPING`
- Remove field: Delete the key-value pair from `CSQ_FIELD_MAPPING`

CSQ fields will be output with the custom names you specify in the mapping. For example, if you map `'Consequence': 'variant_consequence'`, the output column will be `variant_consequence`.

## Performance Optimizations

The script implements several optimizations for high performance:

1. Uses Polars instead of pandas for DataFrame operations
2. Parses CSQ header once at startup
3. Pre-compiles regex patterns
4. Uses efficient string operations and list comprehensions
5. Minimizes object creation in loops
6. Enables Parquet statistics for faster downstream querying
7. Uses native Polars Parquet writer for better performance
8. Pre-computes CSQ field indices before worker loop (O(1) field lookup)
9. Uses set for `seen_genes` tracking (O(1) membership test)
10. Workers return `list[dict]`, not DataFrame (avoids pickling overhead)
11. Each worker opens its own VCF handle (cyvcf2 is not thread-safe)
12. Chromosome-level parallelization for optimal load balancing

## Output Statistics

**Single-threaded version** prints:
- Total variants processed
- Total variant-gene pairs extracted
- Average pairs per variant
- Output file size

**Parallel version** prints:
- Total variants processed (sum across all workers)
- Total variant-gene pairs extracted
- Processing time and rate (pairs/second)
- Output file size
- Per-chromosome statistics (variants and pairs)

### Example Parallel Output

```
[INFO] Loading HGNC mapping: 19,234 entries
[INFO] Parsing VCF header...
[INFO] Found 84 CSQ fields
[INFO] Found 24 standard chromosomes: chr1, chr2, chr3, ..., chrM
[INFO] Starting parallel processing with 8 workers...
[INFO] [chr1] Processing variants...
[INFO] [chr2] Processing variants...
[INFO] [chr1] Processed 10,000 variants...
[INFO] [chr1] Processed 20,000 variants...
[INFO] [chr22] Completed: 12,443 variants → 3,892 pairs extracted
[INFO] [chr1] Completed: 152,443 variants → 47,892 pairs extracted
[INFO] [chr2] Completed: 148,731 variants → 45,123 pairs extracted
...
[INFO] All workers completed
[INFO] Total variant-gene pairs: 945,678
[INFO] Creating Polars DataFrame with 945,678 records...
[INFO] Writing Parquet file...
[INFO] Output written: results.parquet (142.7 MB)
[INFO] Total time: 12.3 seconds
[INFO] Processing rate: 76,910 pairs/second

[INFO] ✓ Processing complete!
```

## Querying the Parquet File

Use the provided `query_vep_parquet.py` script to explore and analyze the output Parquet file.

### Usage

```bash
python query_vep_parquet.py results.parquet
```

### Example Queries Included

The script demonstrates 12 common query patterns:

1. **Filter by chromosome and position**: Get variants from chr1, sorted by position
2. **Gene-specific queries**: Get all variants for a specific gene (e.g., BRCA1)
3. **High-impact variants**: Filter by consequence type (stop_gained, frameshift)
4. **Rare variants**: Filter by allele frequency (gnomAD AF < 1%)
5. **Pathogenic variants**: Filter by ClinVar classification
6. **Genomic regions**: Get variants in specific coordinate ranges
7. **Prediction scores**: Filter by SIFT/PolyPhen predictions
8. **Summary statistics**: Count variants per chromosome
9. **Top genes**: Find genes with most variants
10. **Complex filters**: Combine multiple conditions (rare + pathogenic + protein-altering)
11. **Export subsets**: Save filtered results to new Parquet files
12. **Null handling**: Work with missing values efficiently

### Quick Query Examples

```python
import polars as pl

# Load the Parquet file
df = pl.read_parquet("results.parquet")

# Example 1: Chromosome 1, sorted by position, first 30 records
result = (
    df
    .filter(pl.col("chrom") == "1")
    .sort("pos")
    .head(30)
)

# Example 2: Rare variants (AF < 0.01)
rare = (
    df
    .with_columns(pl.col("gnomadg_af").cast(pl.Float64).alias("af_float"))
    .filter((pl.col("af_float") < 0.01) | pl.col("af_float").is_null())
)

# Example 3: Specific gene variants
brca1 = df.filter(pl.col("gene") == "BRCA1").sort(["chrom", "pos"])

# Example 4: Summary statistics
stats = (
    df
    .group_by("chrom")
    .agg([
        pl.col("variant_id").n_unique().alias("unique_variants"),
        pl.col("gene").n_unique().alias("unique_genes")
    ])
    .sort("chrom")
)
```

## Error Handling

- Gracefully handles malformed CSQ entries with warning messages
- Handles missing/null values throughout
- Validates input files before processing
- Provides informative error messages
