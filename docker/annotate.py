#!/usr/bin/env python3
"""
High-performance parallel VEP-annotated VCF processor using cyvcf2, polars, and multiprocessing.

This script processes VEP-annotated VCF files in parallel by chromosome. For each variant, it:
- Selects the canonical transcript (PICK=1) for detailed annotations
- Extracts the gene from the PICK=1 transcript and stores it in a 'gene' field
- Annotates with OMIM inheritance patterns and phenotype information
Output is in Parquet format with one row per variant.
"""

import argparse
import hashlib
import multiprocessing as mp
import re
import signal
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import polars as pl
from cyvcf2 import VCF

# Compiled regex patterns for performance
INTERVAR_CLASS_PATTERN = re.compile(r'^(Likely pathogenic|Likely benign|Uncertain significance|Pathogenic|Benign)\s+')
INTERVAR_EVIDENCE_PATTERN = re.compile(r'(\w+)=(\[[\d,\s]+\]|\d+)')
INTERVAR_ARRAY_PATTERN = re.compile(r'\d+')


# CSQ field mapping: {original_field_name: output_column_name}
# Edit the output column names as needed
CSQ_FIELD_MAPPING = {
    # Functional annotations
    'Consequence': 'consequence',
    'Feature': 'feature',
    'EXON': 'exon',
    'INTRON': 'intron',
    'HGVSc': 'hgvsc',
    'HGVSp': 'hgvsp',
    'existing_variation': 'rsid',
    'HGNC_ID': 'hgnc_id',
    'MANE': 'mane',

    # Prediction scores
    'SIFT': 'sift',
    'PolyPhen': 'polyPhen',
    'REVEL': 'revel',
    'CLIN_SIG': 'clinsig',
    'am_class': 'am_class',
    'am_pathogenicity': 'am_pathogenicity',
    'CADD_PHRED': 'cadd_phred',
    'CADD_RAW': 'cadd_raw',
    'SpliceAI_pred_DP_AG': 'spliceai_pred_dp_ag',
    'SpliceAI_pred_DP_AL': 'spliceai_pred_dp_al',
    'SpliceAI_pred_DP_DG': 'spliceai_pred_dp_dg',
    'SpliceAI_pred_DP_DL': 'spliceai_pred_dp_dl',
    'SpliceAI_pred_DS_AG': 'spliceai_pred_ds_ag',
    'SpliceAI_pred_DS_AL': 'spliceai_pred_ds_al',
    'SpliceAI_pred_DS_DG': 'spliceai_pred_ds_dg',
    'SpliceAI_pred_DS_DL': 'spliceai_pred_ds_dl',
    'SpliceAI_pred_SYMBOL': 'spliceai_pred_symbol',


    # Population frequencies - 1000 Genomes
    'AF': 'af',
    'AFR_AF': 'afr_af',
    'AMR_AF': 'amr_af',
    'EAS_AF': 'eas_af',
    'EUR_AF': 'eur_af',
    'SAS_AF': 'sas_af',

    # Population frequencies - gnomAD Exomes
    'gnomADe_AF': 'gnomade_af',
    'gnomADe_AFR_AF': 'gnomade_afr_af',
    'gnomADe_AMR_AF': 'gnomade_amr_af',
    'gnomADe_ASJ_AF': 'gnomade_asj_af',
    'gnomADe_EAS_AF': 'gnomade_eaf_af',
    'gnomADe_FIN_AF': 'gnomade_fin_af',
    'gnomADe_MID_AF': 'gnomade_mid_af',
    'gnomADe_NFE_AF': 'gnomade_nfe_af',
    'gnomADe_REMAINING_AF': 'gnomade_remaining_af',
    'gnomADe_SAS_AF': 'gnomade_sas_af',

    # Population frequencies - gnomAD Genomes
    'gnomADg_AF': 'gnomadg_af',
    'gnomADg_AFR_AF': 'gnomadg_afr_af',
    'gnomADg_AMI_AF': 'gnomadg_ami_af',
    'gnomADg_AMR_AF': 'gnomadg_amr_af',
    'gnomADg_ASJ_AF': 'gnomadg_asj_af',
    'gnomADg_EAS_AF': 'gnomadg_eass_af',
    'gnomADg_FIN_AF': 'gnomadg_fin_af',
    'gnomADg_MID_AF': 'gnomadg_mid_af',
    'gnomADg_NFE_AF': 'gnomadg_nfe_af',
    'gnomADg_REMAINING_AF': 'gnomadg_remaining_af',
    'gnomADg_SAS_AF': 'gnomadg_sas_af',

    # Maximum frequencies
    'MAX_AF': 'max_af',
    'MAX_AF_POPS': 'max_af_pops',

    # Transcript selection
    # 'PICK': 'pick',
}


def compute_hash(value: str) -> str:
    """
    Compute MD5 hash of a string value.

    Args:
        value: String to hash

    Returns:
        Hexadecimal hash string (32 characters)
    """
    return hashlib.md5(value.encode()).hexdigest()


def parse_intervar_evidence(evidence_string: str) -> Tuple[Optional[str], Optional[str]]:
    """
    Parse InterVar evidence string to extract ACMG classification and evidences (optimized).

    Example input:
    "InterVar: Benign PVS1=0 PS=[0, 0, 0, 0, 0] PM=[0, 0, 0, 0, 0, 0, 0] PP=[0, 0, 0, 0, 0, 0] BA1=1 BS=[1, 0, 0, 0, 0] BP=[0, 0, 0, 0, 0, 0, 0, 0]"

    Returns:
        (classification, evidences) tuple, e.g., ("Benign", "BA1, BS1")
    """
    if not evidence_string or "InterVar:" not in evidence_string:
        return None, None

    # Remove "InterVar: " prefix - handle potential leading whitespace
    idx = evidence_string.find("InterVar:")
    content = evidence_string[idx + 9:].strip()  # Skip "InterVar:" (9 chars)

    # Extract classification using compiled regex
    match = INTERVAR_CLASS_PATTERN.match(content)
    if not match:
        return None, None

    classification = match.group(1)
    content = content[match.end():]

    # Parse evidence codes using compiled regex
    evidences = []

    for match in INTERVAR_EVIDENCE_PATTERN.finditer(content):
        code = match.group(1)
        value = match.group(2)

        if value.startswith('['):
            # Array format - extract numbers efficiently
            nums = INTERVAR_ARRAY_PATTERN.findall(value)
            for idx, val in enumerate(nums):
                if val != '0':
                    evidences.append(f"{code}{idx+1}")
        elif value != '0':
            # Single non-zero value
            evidences.append(code)

    return classification, ", ".join(evidences) if evidences else None


def load_intervar_data(intervar_file: Path) -> Dict[str, Tuple[str, Optional[str]]]:
    """
    Load InterVar annotation data from TSV file (optimized for large files).

    Args:
        intervar_file: Path to InterVar TSV file

    Returns:
        Dictionary mapping "chrom_pos_ref_alt_gene" to (classification, evidences)
    """
    try:
        print(f"[INFO] Loading InterVar data from {intervar_file}...")

        # Read only required columns by name for efficiency
        # Column names from InterVar: #Chr, Start, Ref, Alt, Ref.Gene, InterVar: InterVar and Evidence
        df = pl.read_csv(
            intervar_file,
            separator="\t",
            has_header=True,
            columns=["#Chr", "Start", "Ref", "Alt", "Ref.Gene", "InterVar: InterVar and Evidence"],
            schema_overrides={
                "#Chr": pl.Utf8,
                "Start": pl.Utf8,
                "Ref": pl.Utf8,
                "Alt": pl.Utf8,
                "Ref.Gene": pl.Utf8,
                "InterVar: InterVar and Evidence": pl.Utf8
            },
            ignore_errors=True
        )

        # Rename columns for easier access
        df = df.rename({
            "#Chr": "Chr",
            "Ref.Gene": "Gene",
            "InterVar: InterVar and Evidence": "InterVar"
        })

        # Filter out rows with missing/invalid genes
        df = df.filter(
            (pl.col("Gene").is_not_null()) &
            (pl.col("Gene") != "NONE") &
            (pl.col("Gene") != ".") &
            (pl.col("Gene") != "")
        )

        # Build mapping dictionary
        mapping = {}
        skipped_count = 0

        for row in df.iter_rows(named=False):
            chrom, start, ref, alt, gene, intervar = row

            if not all([chrom, start, ref, alt, gene, intervar]):
                skipped_count += 1
                continue

            # Normalize chromosome name (remove 'chr' prefix if present for consistent matching)
            chrom_normalized = chrom.replace('chr', '') if chrom.startswith('chr') else chrom

            # Create key with normalized chromosome
            key = f"{chrom_normalized}_{start}_{ref}_{alt}_{gene}"

            # Parse InterVar evidence
            classification, evidences = parse_intervar_evidence(intervar)

            if classification:
                mapping[key] = (classification, evidences)

        print(f"[INFO] Loaded InterVar mapping: {len(mapping):,} entries (skipped {skipped_count:,} invalid rows)")
        if len(mapping) > 0:
            # Print first 3 keys as examples for debugging
            example_keys = list(mapping.keys())[:3]
            print(f"[INFO] Example InterVar keys: {example_keys}")
        return mapping

    except Exception as e:
        print(f"[ERROR] Error loading InterVar data: {e}")
        sys.exit(1)


def load_hgnc_mapping(hgnc_file: Path) -> Dict[str, str]:
    """
    Load HGNC ID to symbol mapping from TSV file.

    Args:
        hgnc_file: Path to TSV file with columns [hgnc_id, hgnc_symbol]

    Returns:
        Dictionary mapping HGNC IDs to symbols
    """
    try:
        # Read HGNC mapping file - only the first two columns
        # First, get column names from header
        with open(hgnc_file, 'r') as f:
            header = f.readline().strip().split('\t')
            col_names = header[:2]

        # Read only the first two columns as strings
        df = pl.read_csv(
            hgnc_file,
            separator="\t",
            columns=col_names,  # Only read first two columns
            schema_overrides={col_names[0]: pl.Utf8, col_names[1]: pl.Utf8}
        )

        # Create mapping from first two columns
        mapping = {
            row[0]: row[1]
            for row in df.iter_rows()
        }
        print(f"[INFO] Loading HGNC mapping: {len(mapping):,} entries")
        return mapping
    except Exception as e:
        print(f"[ERROR] Error loading HGNC mapping: {e}")
        sys.exit(1)


def load_inheritance_mapping(inheritance_file: Path) -> Dict[str, str]:
    """
    Load gene to inheritance mode mapping from TSV file.

    Args:
        inheritance_file: Path to TSV file with columns [gene, inheritance_mode]

    Returns:
        Dictionary mapping gene symbols to inheritance modes
    """
    try:
        df = pl.read_csv(
            inheritance_file,
            separator="\t",
            has_header=False,
            new_columns=["gene", "inheritance"],
            schema_overrides={"gene": pl.Utf8, "inheritance": pl.Utf8}
        )
        mapping = {row[0]: row[1] for row in df.iter_rows()}
        print(f"[INFO] Loading inheritance mapping: {len(mapping):,} entries")
        return mapping
    except Exception as e:
        print(f"[ERROR] Error loading inheritance mapping: {e}")
        sys.exit(1)


def load_omim_phenotype_mapping(omim_file: Path) -> Dict[str, str]:
    """
    Load gene to OMIM phenotype mapping from TSV file.
    Handles multiple phenotypes per gene by concatenating unique values with ' | '.

    Args:
        omim_file: Path to TSV file where column 1 is gene, column 4 is phenotype

    Returns:
        Dictionary mapping gene symbols to concatenated unique phenotypes
    """
    try:
        df = pl.read_csv(
            omim_file,
            separator="\t",
            has_header=False,
            schema_overrides={f"column_{i}": pl.Utf8 for i in range(1, 15)}
        )

        # Group by gene (column_1) and aggregate unique phenotypes (column_4)
        grouped = (
            df.group_by("column_1")
            .agg(pl.col("column_4").unique().str.concat(delimiter=" | "))
        )

        mapping = {row[0]: row[1] for row in grouped.iter_rows()}
        print(f"[INFO] Loading OMIM phenotype mapping: {len(mapping):,} genes")
        return mapping
    except Exception as e:
        print(f"[ERROR] Error loading OMIM phenotype mapping: {e}")
        sys.exit(1)


def parse_csq_header(vcf: VCF) -> List[str]:
    """
    Parse CSQ field structure from VCF header.

    Args:
        vcf: cyvcf2 VCF object

    Returns:
        List of CSQ field names in order
    """
    for header_line in vcf.header_iter():
        info = header_line.info()
        if info.get('ID') == 'CSQ':
            description = info.get('Description', '')
            match = re.search(r'Format:\s*([^\s]+)', description)
            if match:
                fields = match.group(1).split('|')
                return fields

    print("[ERROR] Could not find CSQ field in VCF header")
    sys.exit(1)


def extract_chromosomes(vcf_path: Path) -> List[str]:
    """
    Extract list of standard chromosomes from VCF header.

    Only includes standard chromosomes: 1-22, X, Y, MT (or chr1-chr22, chrX, chrY, chrM).

    Args:
        vcf_path: Path to VCF file

    Returns:
        List of standard chromosome names
    """
    vcf = VCF(str(vcf_path))
    chromosomes = []

    # Define standard chromosomes
    standard_autosomes = set(str(i) for i in range(1, 23))  # 1-22
    standard_autosomes_chr = set(f"chr{i}" for i in range(1, 23))  # chr1-chr22
    standard_sex = {'X', 'Y', 'chrX', 'chrY'}
    standard_mt = {'MT', 'M', 'chrM'}

    standard_chromosomes = standard_autosomes | standard_autosomes_chr | standard_sex | standard_mt

    for header_line in vcf.header_iter():
        info = header_line.info()
        if info.get('HeaderType') == 'CONTIG':
            chrom = info.get('ID')
            if chrom and chrom in standard_chromosomes:
                chromosomes.append(chrom)

    vcf.close()
    return chromosomes


def check_tabix_index(vcf_path: Path) -> None:
    """
    Check if tabix index exists for VCF file.

    Args:
        vcf_path: Path to VCF file

    Raises:
        SystemExit if index doesn't exist
    """
    tbi_path = Path(str(vcf_path) + '.tbi')
    if not tbi_path.exists():
        print(f"[ERROR] Tabix index not found: {tbi_path}")
        print(f"[ERROR] Please create index with: tabix -p vcf {vcf_path}")
        sys.exit(1)


def extract_gene_from_csq(
    csq_values: List[str],
    csq_field_indices: Dict[str, int],
    hgnc_mapping: Dict[str, str]
) -> Optional[str]:
    """
    Extract gene name from CSQ values using priority logic.

    Priority:
    1. HGNC_ID (mapped to symbol)
    2. SYMBOL field

    Args:
        csq_values: List of CSQ field values for one consequence
        csq_field_indices: Dictionary mapping field names to indices (for O(1) lookup)
        hgnc_mapping: HGNC ID to symbol mapping

    Returns:
        Gene symbol or None if no gene found
    """
    # Priority 1: HGNC_ID
    hgnc_idx = csq_field_indices.get('HGNC_ID')
    if hgnc_idx is not None and hgnc_idx < len(csq_values):
        hgnc_id = csq_values[hgnc_idx].strip()
        if hgnc_id:
            hgnc_id = hgnc_id.replace('HGNC:', '')
            if hgnc_id in hgnc_mapping:
                return hgnc_mapping[hgnc_id]

    # Priority 2: SYMBOL
    symbol_idx = csq_field_indices.get('SYMBOL')
    if symbol_idx is not None and symbol_idx < len(csq_values):
        symbol = csq_values[symbol_idx].strip()
        if symbol:
            return symbol

    return None


def process_chromosome(
    chromosome: str,
    vcf_path: str,
    hgnc_mapping: Dict[str, str],
    inheritance_mapping: Dict[str, str],
    omim_phenotype_mapping: Dict[str, str],
    intervar_mapping: Dict[str, Tuple[str, Optional[str]]],
    csq_field_names: List[str]
) -> List[Dict]:
    """
    Process all variants in a single chromosome (worker function).

    Args:
        chromosome: Chromosome name to process
        vcf_path: Path to VCF file
        hgnc_mapping: HGNC ID to symbol mapping
        inheritance_mapping: Gene to inheritance mode mapping
        omim_phenotype_mapping: Gene to OMIM phenotypes mapping
        intervar_mapping: InterVar data mapping (chrom_pos_ref_alt_gene -> (classification, evidences))
        csq_field_names: List of CSQ field names

    Returns:
        List of record dictionaries (one per variant with PICK=1 transcript and its gene)
    """
    # Pre-compute CSQ field indices for O(1) lookup
    csq_field_indices = {name: idx for idx, name in enumerate(csq_field_names)}

    records = []
    variant_count = 0
    acmg_match_count = 0
    acmg_lookup_count = 0

    # Open VCF for this worker (cyvcf2 is not thread-safe, each worker needs its own handle)
    vcf = VCF(vcf_path)

    print(f"[INFO] [{chromosome}] Processing variants...")

    try:
        # Iterate only variants in this chromosome using region query
        for variant in vcf(chromosome):
            variant_count += 1

            # Progress logging every 10k variants
            if variant_count % 10000 == 0:
                print(f"[INFO] [{chromosome}] Processed {variant_count:,} variants...")

            # Extract standard VCF fields
            chrom = variant.CHROM
            pos = variant.POS
            ref = variant.REF
            qual = variant.QUAL if variant.QUAL is not None else None

            # Handle single ALT (VCF contains only one ALT per variant)
            alt = variant.ALT[0] if variant.ALT and len(variant.ALT) > 0 else None
            if not alt:
                continue


            # Extract FORMAT fields (single sample, index 0)
            format_fields = {}
            if len(variant.gt_types) > 0:
                format_fields['gt'] = variant.gt_types[0] if variant.gt_types is not None else None
                format_fields['dp'] = variant.gt_depths[0] if (variant.gt_depths is not None and len(variant.gt_depths) > 0) else None
                format_fields['gq'] = variant.gt_quals[0] if (variant.gt_quals is not None and len(variant.gt_quals) > 0) else None

                if hasattr(variant, 'gt_alt_depths') and variant.gt_alt_depths is not None and len(variant.gt_alt_depths) > 0:
                    format_fields['ad'] = variant.gt_alt_depths[0]
                else:
                    format_fields['ad'] = None
            else:
                # Set all FORMAT fields to None if no genotype data
                format_fields['gt'] = None
                format_fields['dp'] = None
                format_fields['gq'] = None
                format_fields['ad'] = None

            # Create variant identifier
            variant_key = f"{chrom}_{pos}_{ref}_{alt}"
            variant_id = compute_hash(variant_key)

            # Parse CSQ field
            csq_string = variant.INFO.get('CSQ', '')
            if not csq_string:
                continue

            # Split CSQ by comma (multiple consequences)
            csq_entries = csq_string.split(',')

            # Find the PICK=1 transcript
            picked_csq_values = None

            for csq_entry in csq_entries:
                try:
                    # Split by pipe to get individual fields
                    csq_values = csq_entry.split('|')

                    # Filter out CSQ records where BAM_EDIT = FAILED
                    bam_edit_idx = csq_field_indices.get('BAM_EDIT')
                    if bam_edit_idx is not None and bam_edit_idx < len(csq_values):
                        if csq_values[bam_edit_idx] == 'FAILED':
                            continue

                    # Check if this is the picked transcript
                    pick_idx = csq_field_indices.get('PICK')
                    if pick_idx is not None and pick_idx < len(csq_values):
                        if csq_values[pick_idx] == '1':
                            picked_csq_values = csq_values
                            break  # Found PICK=1, stop searching

                except Exception as e:
                    print(f"[WARNING] [{chromosome}] Error parsing CSQ entry at {chrom}:{pos}: {e}")
                    continue

            # Skip if no picked transcript found
            if not picked_csq_values:
                continue

            # Extract gene from the picked transcript only
            gene = extract_gene_from_csq(picked_csq_values, csq_field_indices, hgnc_mapping)

            # Build record from the picked transcript
            record = {
                'chrom': chrom,
                'pos': pos,
                'ref': ref,
                'alt': alt,
                'qual': qual,
                'gene': gene,  # Store single gene from PICK=1 transcript
                'variant_id': variant_id
            }

            # Add FORMAT fields
            record.update(format_fields)

            # Add CSQ fields from picked transcript (only mapped fields from CSQ_FIELD_MAPPING)
            for idx, field_name in enumerate(csq_field_names):
                if field_name in CSQ_FIELD_MAPPING:
                    output_name = CSQ_FIELD_MAPPING[field_name]
                    if idx < len(picked_csq_values):
                        # Convert empty strings to None, keep non-empty values as strings
                        value = picked_csq_values[idx]
                        record[output_name] = value if value and value.strip() else None
                    else:
                        record[output_name] = None

            # Add OMIM annotations for the single gene
            record['omim_inheritance'] = inheritance_mapping.get(gene) if gene else None
            record['omim_phenotype'] = omim_phenotype_mapping.get(gene) if gene else None

            # Add InterVar ACMG annotations for the single gene
            record['ACMG_classification'] = None
            record['ACMG_evidences'] = None

            if gene:
                acmg_lookup_count += 1
                # Normalize chromosome name (remove 'chr' prefix if present for consistent matching)
                chrom_normalized = chrom.replace('chr', '') if chrom.startswith('chr') else chrom

                # Create InterVar lookup key: chrom_pos_ref_alt_gene (with normalized chromosome)
                intervar_key = f"{chrom_normalized}_{pos}_{ref}_{alt}_{gene}"

                if intervar_key in intervar_mapping:
                    acmg_match_count += 1
                    classification, evidences = intervar_mapping[intervar_key]
                    record['ACMG_classification'] = classification
                    record['ACMG_evidences'] = evidences
                elif acmg_lookup_count <= 5:  # Print first 5 failed lookups for debugging
                    print(f"[DEBUG] [{chromosome}] No InterVar match for key: {intervar_key}")

            records.append(record)

    finally:
        vcf.close()

    variants_extracted = len(records)
    print(f"[INFO] [{chromosome}] Completed: {variant_count:,} variants → {variants_extracted:,} variants with PICK=1 extracted")

    # Report ACMG matching statistics
    if acmg_lookup_count > 0:
        match_rate = (acmg_match_count / acmg_lookup_count) * 100
        print(f"[INFO] [{chromosome}] ACMG matches: {acmg_match_count:,}/{acmg_lookup_count:,} ({match_rate:.1f}%)")
    else:
        print(f"[INFO] [{chromosome}] No ACMG lookups performed")

    return records


def process_vcf_parallel(
    vcf_file: Path,
    hgnc_mapping: Dict[str, str],
    inheritance_mapping: Dict[str, str],
    omim_phenotype_mapping: Dict[str, str],
    intervar_mapping: Dict[str, Tuple[str, Optional[str]]],
    output_file: Path,
    num_processes: int
) -> None:
    """
    Process VEP-annotated VCF file in parallel by chromosome.

    Creates one record per variant using the PICK=1 transcript for annotations,
    extracting the gene from the PICK=1 transcript into the 'gene' field.

    Args:
        vcf_file: Path to VEP-annotated VCF file
        hgnc_mapping: HGNC ID to symbol mapping
        inheritance_mapping: Gene to inheritance mode mapping
        omim_phenotype_mapping: Gene to OMIM phenotypes mapping
        intervar_mapping: InterVar data mapping
        output_file: Path for output Parquet file
        num_processes: Number of parallel processes
    """
    start_time = time.time()

    print(f"[INFO] Processing VCF: {vcf_file}")

    # Check tabix index exists
    check_tabix_index(vcf_file)

    # Parse VCF header
    print("[INFO] Parsing VCF header...")
    vcf = VCF(str(vcf_file))
    csq_field_names = parse_csq_header(vcf)
    print(f"[INFO] Found {len(csq_field_names)} CSQ fields")
    vcf.close()

    # Extract chromosomes from VCF header
    chromosomes = extract_chromosomes(vcf_file)
    if not chromosomes:
        print("[ERROR] No chromosomes found in VCF header")
        sys.exit(1)

    print(f"[INFO] Found {len(chromosomes)} standard chromosomes: {', '.join(chromosomes)}")

    # Adjust number of processes if needed (don't create more processes than chromosomes)
    actual_processes = min(num_processes, len(chromosomes))
    if actual_processes < num_processes:
        print(f"[INFO] Adjusting processes from {num_processes} to {actual_processes} (number of chromosomes)")

    print(f"[INFO] Starting parallel processing with {actual_processes} workers...")

    # Create worker arguments (distribute chromosomes across workers)
    worker_args = [
        (chrom, str(vcf_file), hgnc_mapping, inheritance_mapping, omim_phenotype_mapping, intervar_mapping, csq_field_names)
        for chrom in chromosomes
    ]

    # Process chromosomes in parallel using multiprocessing.Pool
    results = []
    try:
        with mp.Pool(processes=actual_processes) as pool:
            # Use pool.starmap to distribute work across workers
            results = pool.starmap(process_chromosome, worker_args)
    except KeyboardInterrupt:
        print("\n[INFO] Interrupted by user, terminating workers...")
        sys.exit(1)

    print("[INFO] All workers completed")

    # Flatten results from all workers into single list
    all_records = [record for chrom_results in results for record in chrom_results]

    total_variants = len(all_records)
    print(f"[INFO] Total variants with PICK=1: {total_variants:,}")

    if not all_records:
        print("[WARNING] No records extracted. Output file will not be created.")
        return

    # Convert to Polars DataFrame once (in main process)
    print(f"[INFO] Creating Polars DataFrame with {len(all_records):,} records...")
    # Use larger infer_schema_length or all records (whichever is smaller) to ensure consistent schema
    schema_length = min(len(all_records), 100000)
    df = pl.DataFrame(all_records, infer_schema_length=schema_length)

    # Write to Parquet with ZSTD compression
    print(f"[INFO] Writing Parquet file...")
    df.write_parquet(
        output_file,
        compression="zstd",
        statistics=True,  # Enable statistics for faster querying
        use_pyarrow=False  # Use native Polars writer for better performance
    )

    file_size_mb = output_file.stat().st_size / (1024 * 1024)
    print(f"[INFO] Output written: {output_file} ({file_size_mb:.1f} MB)")

    # Calculate and print statistics
    elapsed_time = time.time() - start_time
    print(f"[INFO] Total time: {elapsed_time:.1f} seconds")

    if elapsed_time > 0:
        # Calculate processing rate (variants/second)
        rate = total_variants / elapsed_time
        print(f"[INFO] Processing rate: {rate:,.0f} variants/second")


def main():
    """Main entry point for the script."""
    parser = argparse.ArgumentParser(
        description="Process VEP-annotated VCF files in parallel by chromosome",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "--vcf",
        type=Path,
        required=True,
        help="Path to VEP-annotated VCF file (must be bgzipped with .tbi index)"
    )

    parser.add_argument(
        "--hgnc-mapping",
        type=Path,
        required=True,
        help="Path to TSV file with columns [hgnc_id, hgnc_symbol]"
    )

    parser.add_argument(
        "--inheritance-genes",
        type=Path,
        required=True,
        help="Path to TSV file with columns [gene, inheritance_mode]"
    )

    parser.add_argument(
        "--omim-info",
        type=Path,
        required=True,
        help="Path to TSV file with gene and phenotype information"
    )

    parser.add_argument(
        "--intervar-file",
        type=Path,
        required=False,
        help="Path to InterVar TSV file with ACMG annotations (optional)"
    )

    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Path for output Parquet file"
    )

    parser.add_argument(
        "--processes",
        type=int,
        default=mp.cpu_count(),
        help=f"Number of parallel processes (default: {mp.cpu_count()})"
    )

    args = parser.parse_args()

    # Validate inputs
    if not args.vcf.exists():
        print(f"[ERROR] VCF file not found: {args.vcf}")
        sys.exit(1)

    if not args.hgnc_mapping.exists():
        print(f"[ERROR] HGNC mapping file not found: {args.hgnc_mapping}")
        sys.exit(1)

    if not args.inheritance_genes.exists():
        print(f"[ERROR] Inheritance genes file not found: {args.inheritance_genes}")
        sys.exit(1)

    if not args.omim_info.exists():
        print(f"[ERROR] OMIM info file not found: {args.omim_info}")
        sys.exit(1)

    if args.intervar_file and not args.intervar_file.exists():
        print(f"[ERROR] InterVar file not found: {args.intervar_file}")
        sys.exit(1)

    # Create output directory if needed
    args.output.parent.mkdir(parents=True, exist_ok=True)

    # Load mappings (shared across workers for O(1) lookups)
    hgnc_mapping = load_hgnc_mapping(args.hgnc_mapping)
    inheritance_mapping = load_inheritance_mapping(args.inheritance_genes)
    omim_phenotype_mapping = load_omim_phenotype_mapping(args.omim_info)

    # Load InterVar data if provided
    intervar_mapping = {}
    if args.intervar_file:
        intervar_mapping = load_intervar_data(args.intervar_file)
    else:
        print("[INFO] No InterVar file provided, ACMG annotations will be empty")

    # Process VCF in parallel
    process_vcf_parallel(
        args.vcf,
        hgnc_mapping,
        inheritance_mapping,
        omim_phenotype_mapping,
        intervar_mapping,
        args.output,
        args.processes
    )

    print("\n[INFO] ✓ Processing complete!")


if __name__ == "__main__":
    # Handle SIGINT gracefully
    signal.signal(signal.SIGINT, signal.default_int_handler)
    main()
