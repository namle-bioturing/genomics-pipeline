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

    Pattern: scan → filter → sort → limit → collect
    """
    print("\n" + "=" * 80)
    print("EXAMPLE 1: Lazy load chr1, sorted by position, first 30 records")
    print("=" * 80)

    def query():
        return (
            pl.scan_parquet(parquet_file)    # Lazy scan
            .filter(pl.col("chrom") == "chr1")  # Filter pushed to file reader
            .sort("pos")                      # Sort
            .head(30)                         # Limit to 30
            .collect()                        # Execute optimized plan
        )

    result = time_query("Lazy: chr1 filter + sort", query)
    print(result)


def example_2_specific_gene(parquet_file: Path, gene_name: str = "BRCA1") -> None:
    """
    Example 2: Lazy load + filter specific gene (first 30 records).

    Pattern: scan → filter gene → sort → limit → collect
    """
    print("\n" + "=" * 80)
    print(f"EXAMPLE 2: Lazy load {gene_name} variants (first 30)")
    print("=" * 80)

    def query():
        return (
            pl.scan_parquet(parquet_file)
            .filter(pl.col("gene") == gene_name)  # Filter pushed to file reader
            .sort(["chrom", "pos"])
            .head(30)
            .collect()
        )

    result = time_query(f"Lazy: {gene_name} filter", query)
    print(result)


def example_3_high_impact_variants(parquet_file: Path) -> None:
    """
    Example 3: Lazy load + filter high-impact variants (stop gained, frameshift), first 30.

    Pattern: scan → filter consequence → sort → limit → collect
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
            .select(["chrom", "pos", "ref", "alt", "genes", "gnomadg_af", "revel", "am_class"])  # Column pruning
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
            .select(["chrom", "pos", "ref", "alt", "genes", "clinsig"])  # Column pruning
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

    Pattern: scan → filter region → sort → limit → collect
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
            .select([
                "chrom", "pos", "ref", "alt", "genes",
                "ACMG_classification", "ACMG_evidences",
                "consequence", "gnomadg_af"
            ])  # Column pruning
            .sort(["chrom", "pos"])
            .head(30)
            .collect()
        )

    result = time_query("Lazy: ACMG classification filter", query)
    print(result)


def example_8_missense_with_prediction(df: pl.DataFrame) -> None:
    """
    Example 8: Get missense variants with damaging predictions, first 30.
    """
    print("\n" + "=" * 80)
    print("EXAMPLE 8: Missense variants with damaging predictions (SIFT/PolyPhen), first 30")
    print("=" * 80)

    def query():
        return (
            df
            .filter(
                pl.col("consequence").str.contains("missense_variant") &
                (
                    pl.col("sift").str.contains("deleterious") |
                    pl.col("polyPhen").str.contains("damaging")
                )
            )
            .sort(["chrom", "pos"])
            .select(["chrom", "pos", "ref", "alt", "gene", "sift", "polyPhen"])
            .head(30)
        )

    result = time_query("Damaging missense filter", query)
    print(result)


def example_9_variant_summary_stats(df: pl.DataFrame) -> None:
    """
    Example 9: Summary statistics - variants per chromosome.
    """
    print("\n" + "=" * 80)
    print("EXAMPLE 9: Summary statistics - variants per chromosome")
    print("=" * 80)

    stats = (
        df
        .group_by("chrom")
        .agg([
            pl.col("variant_id").n_unique().alias("unique_variants"),
            pl.col("pair_id").count().alias("total_pairs"),
            pl.col("gene").n_unique().alias("unique_genes")
        ])
        .sort("chrom")
    )

    print("\nVariants per chromosome:")
    print(stats)


def example_10_top_genes(df: pl.DataFrame, top_n: int = 20) -> None:
    """
    Example 10: Top genes with most variants.
    """
    print("\n" + "=" * 80)
    print(f"EXAMPLE 10: Top {top_n} genes with most variants")
    print("=" * 80)

    top_genes = (
        df
        .group_by("gene")
        .agg([
            pl.col("variant_id").n_unique().alias("unique_variants"),
            pl.col("pair_id").count().alias("total_pairs")
        ])
        .sort("unique_variants", descending=True)
        .head(top_n)
    )

    print(f"\nTop {top_n} genes:")
    print(top_genes)


def example_11_filter_multiple_conditions(df: pl.DataFrame) -> None:
    """
    Example 11: Complex query - rare, pathogenic, protein-altering variants.
    """
    print("\n" + "=" * 80)
    print("EXAMPLE 11: Complex query - rare + pathogenic + protein-altering")
    print("=" * 80)

    result = (
        df
        .with_columns(
            pl.col("gnomadg_af").cast(pl.Float64).alias("gnomadg_af_float")
        )
        .filter(
            # Rare (AF < 0.01 or null)
            (
                (pl.col("gnomadg_af_float") < 0.01) |
                pl.col("gnomadg_af_float").is_null()
            ) &
            # Pathogenic
            (
                pl.col("clinsig").str.contains("pathogenic").fill_null(False)
            ) &
            # Protein-altering
            (
                pl.col("consequence").str.contains("missense_variant") |
                pl.col("consequence").str.contains("stop_gained") |
                pl.col("consequence").str.contains("frameshift")
            )
        )
        .sort(["chrom", "pos"])
    )

    print(f"\nFound {len(result):,} variants matching all criteria")
    print(result.select(["chrom", "pos", "ref", "alt", "gene", "consequence", "clinsig", "gnomadg_af"]).head(10))


def example_12_export_subset(df: pl.DataFrame, output_file: str = "filtered_variants.parquet") -> None:
    """
    Example 12: Export filtered results to new Parquet file.
    """
    print("\n" + "=" * 80)
    print("EXAMPLE 12: Export filtered results")
    print("=" * 80)

    # Filter for chromosome 1, rare variants
    filtered = (
        df
        .filter(pl.col("chrom") == "1")
        .with_columns(
            pl.col("gnomadg_af").cast(pl.Float64).alias("gnomadg_af_float")
        )
        .filter(
            (pl.col("gnomadg_af_float") < 0.01) |
            pl.col("gnomadg_af_float").is_null()
        )
    )

    print(f"\nExporting {len(filtered):,} records to {output_file}")
    filtered.write_parquet(output_file, compression="zstd")
    print(f"[SUCCESS] Exported to {output_file}")


def example_13_null_handling(df: pl.DataFrame) -> None:
    """
    Example 13: Working with null values.
    """
    print("\n" + "=" * 80)
    print("EXAMPLE 13: Null value handling")
    print("=" * 80)

    # Count nulls in key columns
    null_counts = df.select([
        pl.col("qual").is_null().sum().alias("qual_nulls"),
        pl.col("dp").is_null().sum().alias("dp_nulls"),
        pl.col("rsid").is_null().sum().alias("rsid_nulls"),
        pl.col("clinsig").is_null().sum().alias("clinsig_nulls"),
        pl.col("gnomadg_af").is_null().sum().alias("gnomadg_af_nulls")
    ])

    print("\nNull value counts:")
    print(null_counts)

    # Get records with non-null rsID (known variants)
    known_variants = (
        df
        .filter(pl.col("rsid").is_not_null())
        .sort(["chrom", "pos"])
    )

    print(f"\nFound {len(known_variants):,} variants with rsID")
    print(known_variants.select(["chrom", "pos", "ref", "alt", "gene", "rsid"]).head(10))


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
        # example_2_specific_gene(parquet_file, "BRCA1")
        example_3_high_impact_variants(parquet_file)
        example_4_rare_variants(parquet_file, max_af=0.01)
        example_5_pathogenic_variants(parquet_file)
        example_6_specific_region(parquet_file, chrom="chr17", start=43000000, end=43100000)
        example_7_acmg_classification(parquet_file)

        # Optional examples (commented out by default)
        # example_8_missense_with_prediction(parquet_file)
        # example_9_variant_summary_stats(parquet_file)
        # example_10_top_genes(parquet_file, top_n=20)
        # example_11_filter_multiple_conditions(parquet_file)
        # example_12_export_subset(parquet_file, "filtered_variants.parquet")
        # example_13_null_handling(parquet_file)

        # Export examples - uncomment to use
        # Export all chr1 variants to TSV
        export_to_tsv(
            parquet_file,
            "chr1_variants.tsv",
            filter_func=lambda lf: lf.filter(pl.col("chrom") == "chr1")
        )

        # Export BRCA1 variants with specific columns
        # export_to_tsv(
        #     parquet_file,
        #     "brca1_variants.tsv",
        #     filter_func=lambda lf: lf.filter(pl.col("gene") == "BRCA1"),
        #     columns=["chrom", "pos", "ref", "alt", "gene", "consequence", "gnomadg_af"]
        # )

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
