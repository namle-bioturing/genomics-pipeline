#!/usr/bin/env python3
"""
Query VEP Parquet File Examples

This script demonstrates various queries on the VEP-annotated variant Parquet file
using Polars for high-performance data analysis.

Usage:
    python query_vep_parquet.py results.parquet
"""

import sys
import time
from pathlib import Path
from typing import Optional

import polars as pl


# Standard fields to query from the VEP Parquet file
# This ensures consistent column selection across all examples
QUERY_FIELDS = [
    # Core variant identification
    "chrom",
    "pos",
    "ref",
    "alt",
    # "rsid",
    "qual",

    # VEP annotations
    "gene",
    "consequence",
    "hgvsc",
    "hgvsp",

    # # Population frequency
    # "gnomadg_af",
    # "gnomadg_af_afr",
    # "gnomadg_af_eas",
    # "gnomadg_af_nfe",

    # # Pathogenicity predictions
    # "revel",
    # "sift",
    # "polyphen",
    # "cadd_phred",

    # Clinical significance
    "clinsig",
    "am_class",
    "am_pathogenicity",

    # ACMG classification
    "ACMG_classification",
    "ACMG_evidences",
]


def time_query(query_name: str, query_func):
    """
    Time a query execution and print results.

    Args:
        query_name: Name of the query
        query_func: Function that executes the query and returns result

    Returns:
        Query result
    """
    start_time = time.time()
    result = query_func()
    elapsed = time.time() - start_time

    print(f"\n⏱️  Query time: {elapsed:.4f} seconds ({elapsed*1000:.2f} ms)")
    print(f"📊 Returned {len(result):,} records")

    return result


# Note: All examples use pl.scan_parquet() directly for optimal lazy evaluation
# No need for separate helper functions - each example builds its own query pipeline


def example_1_chromosome_position(parquet_file: Path) -> None:
    """
    Example 1: Lazy load + filter chromosome 1, sorted by position, first 30 records.

    Pattern: scan → filter → select → sort → limit → collect
    """
    print("\n" + "=" * 80)
    print("EXAMPLE 1: Lazy load chr1, sorted by position, first 30 records")
    print("=" * 80)

    def query():
        return (
            pl.scan_parquet(parquet_file)       # Lazy scan
            .filter(pl.col("chrom") == "chr1")  # Filter pushed to file reader
            .select(QUERY_FIELDS)             # Select standard fields
            .sort("pos")                      # Sort
            .head(30)                         # Limit to 30
            .collect()                        # Execute optimized plan
        )

    result = time_query("Lazy: chr1 filter + sort", query)
    print(result)


def example_3_high_impact_variants(parquet_file: Path) -> None:
    """
    Example 3: Lazy load + filter high-impact variants (stop gained, frameshift), first 30.

    Pattern: scan → filter consequence → select → sort → limit → collect
    """
    print("\n" + "=" * 80)
    print("EXAMPLE 3: Lazy load high-impact variants (stop_gained, frameshift), first 30")
    print("=" * 80)

    def query():
        return (
            pl.scan_parquet(parquet_file)
            .filter(
                pl.col("consequence").str.contains("stop_gained") |
                pl.col("consequence").str.contains("frameshift")
            )
            .select(QUERY_FIELDS)           # Select standard fields
            .sort(["chrom", "pos"])
            .head(30)
            .collect()
        )

    result = time_query("Lazy: high-impact filter", query)
    print(result)


def example_4_rare_variants(parquet_file: Path, max_af: float = 0.01) -> None:
    """
    Example 4: Lazy load + filter rare variants (gnomAD AF < 1%), first 30.

    Pattern: scan → filter AF → select columns → sort → limit → collect
    """
    print("\n" + "=" * 80)
    print(f"EXAMPLE 4: Lazy load rare variants (gnomADg_af < {max_af}), first 30")
    print("=" * 80)

    def query():
        return (
            pl.scan_parquet(parquet_file)
            .filter(
                # (pl.col("gnomadg_af").cast(pl.Float64) < max_af) |
                # pl.col("gnomadg_af").is_null()
                pl.col("am_class").is_not_null()
            )
            .select(QUERY_FIELDS)           # Select standard fields
            .sort(["chrom", "pos"])
            .head(30)
            .collect()
        )

    result = time_query("Lazy: rare variants filter", query)
    print(result)


def example_5_pathogenic_variants(parquet_file: Path) -> None:
    """
    Example 5: Lazy load + filter pathogenic/likely pathogenic variants from ClinVar, first 30.

    Pattern: scan → filter clinsig → select columns → sort → limit → collect
    """
    print("\n" + "=" * 80)
    print("EXAMPLE 5: Lazy load pathogenic/likely pathogenic variants (ClinVar), first 30")
    print("=" * 80)

    def query():
        return (
            pl.scan_parquet(parquet_file)
            .filter(
                pl.col("clinsig").str.contains("pathogenic") |
                pl.col("clinsig").str.contains("likely_pathogenic")
            )
            .select(QUERY_FIELDS)           # Select standard fields
            .sort(["chrom", "pos"])
            .head(30)
            .collect()
        )

    result = time_query("Lazy: pathogenic filter", query)
    print(result)


def example_6_specific_region(parquet_file: Path, chrom: str = "chr17",
                               start: int = 43000000, end: int = 43100000) -> None:
    """
    Example 6: Lazy load + filter variants in a specific genomic region, first 30.

    Pattern: scan → filter region → select → sort → limit → collect
    """
    print("\n" + "=" * 80)
    print(f"EXAMPLE 6: Lazy load variants in region {chrom}:{start:,}-{end:,}, first 30")
    print("=" * 80)

    def query():
        return (
            pl.scan_parquet(parquet_file)
            .filter(
                (pl.col("chrom") == chrom) &
                (pl.col("pos") >= start) &
                (pl.col("pos") <= end)
            )
            .select(QUERY_FIELDS)           # Select standard fields
            .sort("pos")
            .head(30)
            .collect()
        )

    result = time_query("Lazy: region filter", query)
    print(result)


def example_7_acmg_classification(parquet_file: Path) -> None:
    """
    Example 7: Lazy load + filter variants with ACMG classification, first 30.

    Pattern: scan → filter ACMG → select columns → sort → limit → collect
    """
    print("\n" + "=" * 80)
    print("EXAMPLE 7: Lazy load variants with ACMG classification, first 30")
    print("=" * 80)

    def query():
        return (
            pl.scan_parquet(parquet_file)
            .filter(pl.col("ACMG_classification").is_not_null())
            .select(QUERY_FIELDS)           # Select standard fields
            .sort(["chrom", "pos"])
            .head(30)
            .collect()
        )

    result = time_query("Lazy: ACMG classification filter", query)
    print(result)


def example_8_duplicate_variants(parquet_file: Path, min_occurrences: int = 2) -> None:
    """
    Example 8: Find variants that appear multiple times (default >=2 occurrences).

    Pattern: scan → group by variant keys → filter counts → join back → select → collect
    """
    print("\n" + "=" * 80)
    print(f"EXAMPLE 8: Variants appearing at least {min_occurrences} times")
    print("=" * 80)

    duplicates = (
        pl.scan_parquet(parquet_file)
        .group_by(["chrom", "pos", "ref", "alt"])
        .agg(pl.len().alias("variant_count"))
        .filter(pl.col("variant_count") >= min_occurrences)
    )

    def query():
        return (
            pl.scan_parquet(parquet_file)
            .join(duplicates, on=["chrom", "pos", "ref", "alt"], how="inner")
            .select(["chrom", "pos", "ref", "alt", "gene", "variant_count"])
            .sort(["variant_count", "chrom", "pos"], descending=[True, False, False])
            .head(30)
            .collect()
        )

    result = time_query("Lazy: duplicate variant finder", query)
    print(result)


def print_all_columns(parquet_file: Path) -> None:
    """
    Print all available columns in the Parquet file.
    """
    print("\n" + "=" * 80)
    print("Available Columns in Dataset")
    print("=" * 80)

    # Scan just the schema without loading data
    df = pl.scan_parquet(parquet_file)
    columns = df.collect_schema().names()

    print(f"\nTotal columns: {len(columns)}")
    print("\nColumn names:")
    for i, col in enumerate(columns, 1):
        print(f"  {i:2d}. {col}")


def export_to_tsv(
    parquet_file: Path,
    output_file: str = "filtered_variants.tsv",
    filter_func=None,
    columns: Optional[list[str]] = None
) -> None:
    """
    Export filtered results to TSV file.

    Args:
        parquet_file: Input Parquet file path
        output_file: Output TSV file path
        filter_func: Optional function that takes a LazyFrame and returns a filtered LazyFrame
        columns: Optional list of columns to export (exports all if None)

    Example usage:
        # Export all data
        export_to_tsv(parquet_file, "all_variants.tsv")

        # Export with custom filter
        def my_filter(lf):
            return lf.filter(pl.col("chrom") == "chr1")
        export_to_tsv(parquet_file, "chr1_variants.tsv", filter_func=my_filter)

        # Export specific columns
        export_to_tsv(parquet_file, "subset.tsv",
                     columns=["chrom", "pos", "ref", "alt", "gene"])
    """
    print("\n" + "=" * 80)
    print("Export to TSV")
    print("=" * 80)

    # Start with lazy scan
    lf = pl.scan_parquet(parquet_file)

    # Apply filter if provided
    if filter_func:
        lf = filter_func(lf)

    # Select columns if specified
    if columns:
        lf = lf.select(columns)

    # Collect and write to TSV
    print(f"\n[INFO] Processing data...")
    df = lf.collect()

    # Handle nested columns (lists, structs) by converting to JSON strings
    # CSV/TSV format doesn't support nested data structures
    print(f"[INFO] Converting nested columns to JSON strings for CSV export...")

    # Identify list/struct columns and convert them to JSON strings
    nested_cols = []
    for col_name in df.columns:
        dtype = df[col_name].dtype
        if isinstance(dtype, (pl.List, pl.Struct)):
            nested_cols.append(col_name)

    if nested_cols:
        print(f"[INFO] Found nested columns: {', '.join(nested_cols)}")
        # Convert nested columns to JSON strings
        # We need to convert to Python objects first, then serialize to JSON
        import json

        for col in nested_cols:
            dtype = df[col].dtype
            if isinstance(dtype, pl.List):
                # For list columns: convert to JSON array string
                # Extract as Python list, serialize to JSON string
                python_data = df[col].to_list()
                json_strings = [json.dumps(item) if item is not None else None for item in python_data]
                df = df.with_columns(
                    pl.Series(col, json_strings, dtype=pl.Utf8)
                )
            elif isinstance(dtype, pl.Struct):
                # For struct columns: use built-in json_encode
                df = df.with_columns(
                    pl.col(col).struct.json_encode().alias(col)
                )

    print(f"[INFO] Exporting {len(df):,} records to {output_file}")
    df.write_csv(output_file, separator="\t")

    print(f"[SUCCESS] Exported to {output_file}")

    # Show preview
    print(f"\nPreview (first 5 rows):")
    print(df.head(5))


def main():
    """
    Main function to run all query examples.

    All examples use the OPTIMAL pattern for large files:
    - Lazy scan (no upfront data loading)
    - Filter during file read (predicate pushdown)
    - Select only needed columns (projection pushdown)
    - Collect only matching records

    This is 5-10x faster than loading all data first!
    """
    if len(sys.argv) < 2:
        print("Usage: python query_vep_parquet.py <parquet_file>")
        print("\nExample:")
        print("  python query_vep_parquet.py results.parquet")
        print("\nOptimal Pattern (used by ALL examples):")
        print("  ✓ Lazy scan - no upfront data loading")
        print("  ✓ Filter during file read (predicate pushdown)")
        print("  ✓ Select only needed columns (projection pushdown)")
        print("  ✓ Collect only matching records")
        print("\nBenefits for 5M+ rows:")
        print("  • 5-10x faster than eager loading all data")
        print("  • 90% less memory usage")
        print("  • Automatic query plan optimization")
        sys.exit(1)

    parquet_file = Path(sys.argv[1])

    if not parquet_file.exists():
        print(f"[ERROR] File not found: {parquet_file}")
        sys.exit(1)

    # Check file size
    file_size_mb = parquet_file.stat().st_size / (1024 * 1024)
    print(f"\n{'='*80}")
    print(f"VEP Parquet Query Examples - Lazy Load + Filter Pattern")
    print(f"{'='*80}")
    print(f"[INFO] Parquet file: {parquet_file}")
    print(f"[INFO] File size: {file_size_mb:.1f} MB")
    print(f"\n[INFO] All examples use lazy evaluation for optimal performance")
    print(f"[INFO] No data loaded until each query executes")

    # Run examples - all use lazy load + filter pattern
    print("\n" + "=" * 80)
    print("🚀 Running Examples - All use LAZY LOAD + FILTER for optimal performance")
    print("=" * 80)

    try:
        # First, show all available columns
        print_all_columns(parquet_file)

        example_1_chromosome_position(parquet_file)
        example_3_high_impact_variants(parquet_file)
        example_4_rare_variants(parquet_file, max_af=0.01)
        example_5_pathogenic_variants(parquet_file)
        example_6_specific_region(parquet_file, chrom="chr17", start=43000000, end=43100000)
        example_7_acmg_classification(parquet_file)
        example_8_duplicate_variants(parquet_file, min_occurrences=2)

        # Export all chr1 variants to TSV
        export_to_tsv(
            parquet_file,
            "chr1_variants.tsv",
            filter_func=lambda lf: lf.filter(pl.col("chrom") == "chr1")
        )

        print("\n" + "=" * 80)
        print("✓ All examples completed successfully!")
        print("=" * 80)

    except Exception as e:
        print(f"\n[ERROR] An error occurred: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
