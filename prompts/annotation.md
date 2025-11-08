# Clean Prompt for Claude Code Generation

Create a high-performance Python script using cyvcf2, polars, and multiprocessing to process VEP-annotated VCF files in parallel by chromosome.

## Requirements:

### Input:
- Single-sample VCF file annotated with VEP (.vcf.gz format with tabix index required)
- VCF contains CSQ field in INFO with pipe-delimited consequences

### Gene Extraction Logic (Critical):
For each variant, extract variant-gene pairs using this priority:
1. **HGNC Priority**: If HGNC_ID exists in CSQ record, map it to HGNC symbol using provided mapping dict
2. **Symbol Fallback**: If no HGNC_ID, use SYMBOL field from CSQ
3. **Deduplication**: Track genes per variant - only output each gene once per variant
4. **Exclusion**: Skip records with no gene OR if gene was already seen for this variant

Example:
- Variant A has 4 CSQ records:
  - Record 1: HGNC_ID=1 (maps to "GENE_E") → Extract: Variant_A - GENE_E
  - Record 2: No HGNC_ID, SYMBOL="GENE_F" → Extract: Variant_A - GENE_F  
  - Record 3: No HGNC_ID, No SYMBOL → Skip
  - Record 4: No HGNC_ID, SYMBOL="GENE_E" → Skip (already extracted)
- Result: 2 records output

### Fields to Extract from VCF:
**Standard columns:**
- CHROM, POS, REF, ALT, QUAL, FILTER
- Only 1 ALT exist, we don't have multiple ALT

**INFO fields:**
- All INFO fields available (flexible - may include DP, AC, AF, etc.)

**FORMAT fields (single sample):**
- GT (genotype)
- DP (read depth)
- GQ (genotype quality)
- AD (allelic depth)
- Any other FORMAT fields present

**CSQ fields to extract (from pipe-delimited string):**
Dynamically from the VCF header

### Filtering Requirements: No filter

### Parallelization Strategy:

**Architecture:**
```
Main Process
├── Check VCF has tabix index (.tbi file must exist)
├── Parse VCF header (extract chromosomes list)
├── Load HGNC mapping (shared across workers)
├── Parse CSQ field structure from header (shared across workers)
└── Spawn worker processes via multiprocessing.Pool
    ├── Worker 1: Process chr1
    ├── Worker 2: Process chr2
    ├── Worker 3: Process chr3
    └── ... (distribute chromosomes across workers)
        ↓
    Each worker:
    - Opens its own VCF handle with cyvcf2
    - Iterates only variants in assigned chromosome
    - Returns list of record dicts (NOT DataFrame)
        ↓
Main Process:
├── Collect all worker results
├── Flatten into single list
├── Convert to Polars DataFrame once
├── Add computed columns
└── Write single Parquet file
```

**Implementation Requirements:**

1. **Chromosome extraction:**
   - Read VCF header to get list of chromosomes from contigs
   - Support both "chr1" and "1" naming formats
   - Include all chromosomes present in VCF (autosomes, sex chromosomes, MT, etc.)

2. **Worker function:**
   - Takes: chromosome name, vcf_path, hgnc_map, csq_field_names
   - Opens VCF independently in each worker (cyvcf2 is NOT thread-safe)
   - Uses region query: `for variant in vcf(chrom):` to iterate only assigned chromosome
   - Returns plain list of dicts (avoids DataFrame pickling overhead)
   - Logs progress every 10k variants within worker

3. **Shared data handling:**
   - HGNC mapping dict is passed to workers (will be pickled)
   - CSQ field names list is passed to workers
   - Each worker opens its own VCF file handle from the path

4. **Result combination:**
   - Use `pool.map()` to distribute chromosomes to workers
   - Collect all worker results: `results = pool.map(worker_func, chromosomes)`
   - Flatten: `all_records = [rec for chrom_results in results for rec in chrom_results]`
   - Build single Polars DataFrame from flattened list
   - Add computed columns in main process
   - Write single Parquet file

5. **Process management:**
   - Use `multiprocessing.Pool` with context manager or explicit close/join
   - Handle KeyboardInterrupt gracefully (terminate workers)
   - Number of processes: command-line argument (default: cpu_count())
   - Adjust if processes > number of chromosomes

### Command Line Interface:

```bash
python process_vep_parallel.py \
  --vcf annotated.vcf.gz \
  --hgnc-mapping hgnc_map.tsv \
  --output results.parquet \
  --processes 8
```

**Arguments:**
- `--vcf`: Path to VEP-annotated VCF (must be bgzipped with .tbi index)
- `--hgnc-mapping`: Path to HGNC mapping TSV (columns: hgnc_id, hgnc_symbol)
- `--output`: Path for output Parquet file
- `--processes`: Number of parallel processes (default: cpu_count())

### Output Format:
- Parquet file with ZSTD compression
- Each row represents one variant-gene pair
- The field GENE can be extracted from mapped from HGNC_ID if exist, or else using the SYMBOL
- Include all extracted VCF fields + selected CSQ fields
- Add computed columns:
  - variant_id: hash("{CHROM}_{POS}_{REF}_{ALT}")
  - pair_id: hash("{CHROM}_{POS}_{REF}_{ALT}_{GENE}")

Use hashlib.md5 for consistent hashing:
```python
import hashlib
def generate_hash(value: str) -> str:
    return hashlib.md5(value.encode()).hexdigest()
```

### Performance Optimizations Required:
1. Use Polars (not pandas) for DataFrame operations
2. Parse CSQ header once at start (don't repeat)
3. Pre-compile any regex patterns
4. Use list comprehension where possible
5. Minimize object creation in loops
6. Use efficient string operations
7. Enable Parquet statistics for faster querying
8. **Pre-compute CSQ field indices before worker loop** (use dict for O(1) lookup)
9. **Use set for seen_genes tracking** (O(1) membership test)
10. **Workers return list[dict], not DataFrame** (avoid pickling overhead)

### Code Structure:
```
1. Parse command-line arguments
2. Check VCF has tabix index
3. Load HGNC mapping (dict: hgnc_id -> hgnc_symbol)
4. Parse CSQ field structure from VCF header
5. Extract chromosome list from VCF header
6. Create worker function that:
   - Takes chromosome, vcf_path, hgnc_map, csq_field_names
   - Opens VCF with cyvcf2 for that chromosome
   - Iterates variants in assigned chromosome
   - Extracts all VCF fields
   - Parses CSQ string and applies gene extraction logic
   - Returns list of record dicts
7. Create multiprocessing.Pool
8. Map worker function across chromosomes
9. Collect and flatten results from all workers
10. Convert to Polars DataFrame
11. Add computed columns (variant_id, pair_id hashes)
12. Write to Parquet with compression
```

### Input Parameters:
- vcf_file: Path to VEP-annotated VCF (can be .vcf.gz)
- hgnc_mapping_file: TSV with columns [hgnc_id, hgnc_symbol]
- output_parquet: Path for output Parquet file
- processes: Number of parallel processes

### Example Usage:
```python
python process_vep.py \
  --vcf annotated.vcf.gz \
  --hgnc-mapping hgnc_map.tsv \
  --output results.parquet \
  --processes 8
```

### Additional Requirements:
- Add progress logging (print stats every 10k variants within workers)
- Handle missing/null values gracefully
- Add error handling for malformed CSQ fields
- Print summary statistics at end (total variants, total pairs extracted)
- Use type hints
- Add docstrings
- Verify tabix index exists before starting
- Log from each worker with chromosome prefix: `[chr1] Processing...`
- Log total execution time and processing rate (variants/second)

### Logging Example:
```
[INFO] Loading HGNC mapping: 19,234 entries
[INFO] Parsing VCF header...
[INFO] Found 24 chromosomes: chr1, chr2, ..., chrM
[INFO] Starting parallel processing with 8 workers...
[INFO] [chr1] Processing variants...
[INFO] [chr2] Processing variants...
[INFO] [chr1] Processed 10,000 variants...
[INFO] [chr1] Completed: 152,443 variants → 47,892 pairs extracted
[INFO] [chr2] Completed: 148,731 variants → 45,123 pairs extracted
...
[INFO] All workers completed
[INFO] Total variants processed: 3,088,269
[INFO] Total variant-gene pairs: 945,678
[INFO] Writing Parquet file...
[INFO] Output written: results.parquet (142.7 MB)
[INFO] Total time: 12.3 seconds
[INFO] Processing rate: 251,042 variants/second
```

Generate clean, production-ready code optimized for speed with parallel processing by chromosome.