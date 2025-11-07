#!/usr/bin/env python3
"""
VEP Parquet API Server

FastAPI server for querying VEP-annotated variant Parquet files.
Uses lazy loading with Polars for optimal performance on large datasets.

Usage:
    python vep_api.py --parquet results.parquet --host 0.0.0.0 --port 8000

API Documentation:
    Once running, visit: http://localhost:8000/docs
"""

import argparse
import sys
from pathlib import Path
from typing import Optional, List

import polars as pl
import uvicorn
from fastapi import FastAPI, HTTPException, Query
from fastapi.responses import JSONResponse
from pydantic import BaseModel, Field


# Global variable to store parquet file path
PARQUET_FILE: Optional[Path] = None


# ============================================================================
# Pydantic Models for API Responses
# ============================================================================

class VariantResponse(BaseModel):
    """Single variant record."""
    chrom: str
    pos: int
    ref: str
    alt: str
    gene: Optional[str] = None
    variant_id: str
    pair_id: str
    qual: Optional[float] = None
    consequence: Optional[str] = None
    clinsig: Optional[str] = None
    gnomadg_af: Optional[str] = None
    hgvsc: Optional[str] = None
    hgvsp: Optional[str] = None
    sift: Optional[str] = None
    polyphen: Optional[str] = None


class VariantsListResponse(BaseModel):
    """List of variants with metadata."""
    total_count: int
    returned_count: int
    variants: List[dict]


class StatsResponse(BaseModel):
    """Statistics response."""
    stats: dict


# ============================================================================
# FastAPI App
# ============================================================================

app = FastAPI(
    title="VEP Parquet API",
    description="REST API for querying VEP-annotated variants from Parquet files",
    version="1.0.0",
    docs_url="/docs",
    redoc_url="/redoc"
)


# ============================================================================
# Helper Functions
# ============================================================================

def get_parquet_path() -> Path:
    """Get the Parquet file path."""
    if PARQUET_FILE is None:
        raise HTTPException(status_code=500, detail="Parquet file not configured")
    return PARQUET_FILE


def execute_query(lazy_frame: pl.LazyFrame, limit: int = 100) -> dict:
    """
    Execute a lazy query and return results as dict.

    Args:
        lazy_frame: Polars LazyFrame with filters applied
        limit: Maximum number of records to return

    Returns:
        Dictionary with total_count, returned_count, and variants
    """
    # Execute query
    df = lazy_frame.head(limit).collect()

    # Convert to list of dicts
    variants = df.to_dicts()

    return {
        "total_count": len(df),
        "returned_count": len(variants),
        "variants": variants
    }


# ============================================================================
# API Endpoints
# ============================================================================

@app.get("/", tags=["General"])
async def root():
    """API root endpoint."""
    return {
        "message": "VEP Parquet API",
        "version": "1.0.0",
        "docs": "/docs",
        "parquet_file": str(PARQUET_FILE) if PARQUET_FILE else "Not configured"
    }


@app.get("/health", tags=["General"])
async def health_check():
    """Health check endpoint."""
    parquet_path = get_parquet_path()

    if not parquet_path.exists():
        raise HTTPException(status_code=503, detail="Parquet file not found")

    return {
        "status": "healthy",
        "parquet_file": str(parquet_path),
        "file_size_mb": parquet_path.stat().st_size / (1024 * 1024)
    }


@app.get("/variants/chromosome/{chrom}", tags=["Variants"], response_model=VariantsListResponse)
async def get_variants_by_chromosome(
    chrom: str = Path(..., description="Chromosome name (e.g., chr1, 1, chrX)"),
    limit: int = Query(100, ge=1, le=10000, description="Maximum number of records to return"),
    min_pos: Optional[int] = Query(None, description="Minimum position"),
    max_pos: Optional[int] = Query(None, description="Maximum position")
):
    """
    Get variants from a specific chromosome.

    **Optimal Pattern:** Lazy scan → filter chromosome → (optional position filter) → limit → collect
    """
    parquet_path = get_parquet_path()

    # Build lazy query
    lf = pl.scan_parquet(parquet_path).filter(pl.col("chrom") == chrom)

    # Apply position filters if provided
    if min_pos is not None:
        lf = lf.filter(pl.col("pos") >= min_pos)
    if max_pos is not None:
        lf = lf.filter(pl.col("pos") <= max_pos)

    # Sort by position
    lf = lf.sort("pos")

    return execute_query(lf, limit)


@app.get("/variants/gene/{gene}", tags=["Variants"], response_model=VariantsListResponse)
async def get_variants_by_gene(
    gene: str = Path(..., description="Gene symbol (e.g., BRCA1)"),
    limit: int = Query(100, ge=1, le=10000, description="Maximum number of records to return")
):
    """
    Get all variants for a specific gene.

    **Optimal Pattern:** Lazy scan → filter gene → sort → limit → collect
    """
    parquet_path = get_parquet_path()

    # Build lazy query
    lf = (
        pl.scan_parquet(parquet_path)
        .filter(pl.col("gene") == gene)
        .sort(["chrom", "pos"])
    )

    return execute_query(lf, limit)


@app.get("/variants/region", tags=["Variants"], response_model=VariantsListResponse)
async def get_variants_by_region(
    chrom: str = Query(..., description="Chromosome name"),
    start: int = Query(..., ge=0, description="Start position"),
    end: int = Query(..., ge=0, description="End position"),
    limit: int = Query(100, ge=1, le=10000, description="Maximum number of records to return")
):
    """
    Get variants in a specific genomic region.

    **Optimal Pattern:** Lazy scan → filter region → sort → limit → collect
    """
    parquet_path = get_parquet_path()

    if end <= start:
        raise HTTPException(status_code=400, detail="End position must be greater than start position")

    # Build lazy query
    lf = (
        pl.scan_parquet(parquet_path)
        .filter(
            (pl.col("chrom") == chrom) &
            (pl.col("pos") >= start) &
            (pl.col("pos") <= end)
        )
        .sort("pos")
    )

    return execute_query(lf, limit)


@app.get("/variants/pathogenic", tags=["Variants"], response_model=VariantsListResponse)
async def get_pathogenic_variants(
    limit: int = Query(100, ge=1, le=10000, description="Maximum number of records to return"),
    chrom: Optional[str] = Query(None, description="Filter by chromosome")
):
    """
    Get pathogenic or likely pathogenic variants (from ClinVar).

    **Optimal Pattern:** Lazy scan → filter clinsig → (optional chrom filter) → sort → limit → collect
    """
    parquet_path = get_parquet_path()

    # Build lazy query
    lf = pl.scan_parquet(parquet_path).filter(
        pl.col("clinsig").str.contains("pathogenic") |
        pl.col("clinsig").str.contains("likely_pathogenic")
    )

    # Optional chromosome filter
    if chrom:
        lf = lf.filter(pl.col("chrom") == chrom)

    # Select relevant columns and sort
    lf = (
        lf
        .select(["chrom", "pos", "ref", "alt", "gene", "consequence", "clinsig", "hgvsc", "variant_id", "pair_id"])
        .sort(["chrom", "pos"])
    )

    return execute_query(lf, limit)


@app.get("/variants/rare", tags=["Variants"], response_model=VariantsListResponse)
async def get_rare_variants(
    max_af: float = Query(0.01, ge=0.0, le=1.0, description="Maximum allele frequency (default: 0.01)"),
    limit: int = Query(100, ge=1, le=10000, description="Maximum number of records to return"),
    chrom: Optional[str] = Query(None, description="Filter by chromosome")
):
    """
    Get rare variants (gnomAD AF below threshold or null).

    **Optimal Pattern:** Lazy scan → filter AF → (optional chrom filter) → sort → limit → collect
    """
    parquet_path = get_parquet_path()

    # Build lazy query
    lf = pl.scan_parquet(parquet_path).filter(
        (pl.col("gnomadg_af").cast(pl.Float64) < max_af) |
        pl.col("gnomadg_af").is_null()
    )

    # Optional chromosome filter
    if chrom:
        lf = lf.filter(pl.col("chrom") == chrom)

    # Select relevant columns and sort
    lf = (
        lf
        .select(["chrom", "pos", "ref", "alt", "gene", "gnomadg_af", "consequence", "variant_id", "pair_id"])
        .sort(["chrom", "pos"])
    )

    return execute_query(lf, limit)


@app.get("/variants/high-impact", tags=["Variants"], response_model=VariantsListResponse)
async def get_high_impact_variants(
    limit: int = Query(100, ge=1, le=10000, description="Maximum number of records to return"),
    chrom: Optional[str] = Query(None, description="Filter by chromosome")
):
    """
    Get high-impact variants (stop_gained, frameshift).

    **Optimal Pattern:** Lazy scan → filter consequence → (optional chrom filter) → sort → limit → collect
    """
    parquet_path = get_parquet_path()

    # Build lazy query
    lf = pl.scan_parquet(parquet_path).filter(
        pl.col("consequence").str.contains("stop_gained") |
        pl.col("consequence").str.contains("frameshift")
    )

    # Optional chromosome filter
    if chrom:
        lf = lf.filter(pl.col("chrom") == chrom)

    # Sort
    lf = lf.sort(["chrom", "pos"])

    return execute_query(lf, limit)


@app.get("/variants/search", tags=["Variants"], response_model=VariantsListResponse)
async def search_variants(
    chrom: Optional[str] = Query(None, description="Chromosome"),
    gene: Optional[str] = Query(None, description="Gene symbol"),
    consequence: Optional[str] = Query(None, description="Consequence type (contains)"),
    max_af: Optional[float] = Query(None, ge=0.0, le=1.0, description="Maximum allele frequency"),
    pathogenic: bool = Query(False, description="Only pathogenic variants"),
    limit: int = Query(100, ge=1, le=10000, description="Maximum number of records to return")
):
    """
    Search variants with multiple optional filters.

    **Optimal Pattern:** Lazy scan → apply all filters → sort → limit → collect
    """
    parquet_path = get_parquet_path()

    # Build lazy query
    lf = pl.scan_parquet(parquet_path)

    # Apply filters
    if chrom:
        lf = lf.filter(pl.col("chrom") == chrom)

    if gene:
        lf = lf.filter(pl.col("gene") == gene)

    if consequence:
        lf = lf.filter(pl.col("consequence").str.contains(consequence))

    if max_af is not None:
        lf = lf.filter(
            (pl.col("gnomadg_af").cast(pl.Float64) < max_af) |
            pl.col("gnomadg_af").is_null()
        )

    if pathogenic:
        lf = lf.filter(
            pl.col("clinsig").str.contains("pathogenic") |
            pl.col("clinsig").str.contains("likely_pathogenic")
        )

    # Sort
    lf = lf.sort(["chrom", "pos"])

    return execute_query(lf, limit)


@app.get("/stats/chromosomes", tags=["Statistics"], response_model=StatsResponse)
async def get_chromosome_stats():
    """
    Get variant statistics per chromosome.

    **Note:** This loads filtered data for aggregation.
    """
    parquet_path = get_parquet_path()

    # Build query
    df = (
        pl.scan_parquet(parquet_path)
        .group_by("chrom")
        .agg([
            pl.col("variant_id").n_unique().alias("unique_variants"),
            pl.col("pair_id").count().alias("total_pairs"),
            pl.col("gene").n_unique().alias("unique_genes")
        ])
        .sort("chrom")
        .collect()
    )

    return {"stats": df.to_dict(as_series=False)}


@app.get("/stats/genes", tags=["Statistics"], response_model=StatsResponse)
async def get_top_genes(
    top_n: int = Query(20, ge=1, le=100, description="Number of top genes to return")
):
    """
    Get top genes with most variants.

    **Note:** This loads data for aggregation.
    """
    parquet_path = get_parquet_path()

    # Build query
    df = (
        pl.scan_parquet(parquet_path)
        .group_by("gene")
        .agg([
            pl.col("variant_id").n_unique().alias("unique_variants"),
            pl.col("pair_id").count().alias("total_pairs")
        ])
        .sort("unique_variants", descending=True)
        .head(top_n)
        .collect()
    )

    return {"stats": df.to_dict(as_series=False)}


# ============================================================================
# Main Entry Point
# ============================================================================

def main():
    """Main function to start the API server."""
    parser = argparse.ArgumentParser(description="VEP Parquet API Server")
    parser.add_argument(
        "--parquet",
        type=str,
        required=True,
        help="Path to VEP-annotated Parquet file"
    )
    parser.add_argument(
        "--host",
        type=str,
        default="0.0.0.0",
        help="Host to bind to (default: 0.0.0.0)"
    )
    parser.add_argument(
        "--port",
        type=int,
        default=8000,
        help="Port to bind to (default: 8000)"
    )
    parser.add_argument(
        "--reload",
        action="store_true",
        help="Enable auto-reload for development"
    )

    args = parser.parse_args()

    # Set global parquet file path
    global PARQUET_FILE
    PARQUET_FILE = Path(args.parquet)

    # Validate parquet file exists
    if not PARQUET_FILE.exists():
        print(f"[ERROR] Parquet file not found: {PARQUET_FILE}")
        sys.exit(1)

    file_size_mb = PARQUET_FILE.stat().st_size / (1024 * 1024)
    print(f"\n{'='*80}")
    print(f"VEP Parquet API Server")
    print(f"{'='*80}")
    print(f"[INFO] Parquet file: {PARQUET_FILE}")
    print(f"[INFO] File size: {file_size_mb:.1f} MB")
    print(f"[INFO] Host: {args.host}")
    print(f"[INFO] Port: {args.port}")
    print(f"\n[INFO] Starting server...")
    print(f"[INFO] API Documentation: http://{args.host}:{args.port}/docs")
    print(f"[INFO] ReDoc: http://{args.host}:{args.port}/redoc")
    print(f"{'='*80}\n")

    # Start server
    uvicorn.run(
        "vep_api:app",
        host=args.host,
        port=args.port,
        reload=args.reload,
        log_level="info"
    )


if __name__ == "__main__":
    main()
