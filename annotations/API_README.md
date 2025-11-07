# VEP Parquet REST API

Fast REST API for querying VEP-annotated variants from Parquet files using lazy evaluation for optimal performance.

## Features

- ✅ **Lazy Loading**: Uses `scan_parquet` for 5-10x faster queries on large files
- ✅ **RESTful API**: Standard HTTP endpoints with JSON responses
- ✅ **Auto Documentation**: Interactive Swagger UI and ReDoc
- ✅ **High Performance**: Polars-powered queries with predicate/projection pushdown
- ✅ **Flexible Filtering**: Multiple filter combinations for variants
- ✅ **Statistics**: Aggregate data by chromosome and gene

## Installation

```bash
pip install -r requirements.txt
```

Required dependencies:
- `fastapi>=0.104.0`
- `uvicorn[standard]>=0.24.0`
- `pydantic>=2.0.0`
- `polars>=0.19.0`

## Quick Start

### 1. Start the API Server

```bash
python vep_api.py --parquet results.parquet --host 0.0.0.0 --port 8000
```

**Arguments:**
- `--parquet`: Path to VEP-annotated Parquet file (required)
- `--host`: Host to bind to (default: `0.0.0.0`)
- `--port`: Port to bind to (default: `8000`)
- `--reload`: Enable auto-reload for development (optional)

### 2. Access the API

Once running:
- **API Documentation**: http://localhost:8000/docs (Swagger UI)
- **Alternative Docs**: http://localhost:8000/redoc (ReDoc)
- **Health Check**: http://localhost:8000/health

## API Endpoints

### General

#### `GET /` - Root
Returns API information and status.

```bash
curl http://localhost:8000/
```

#### `GET /health` - Health Check
Check API health and Parquet file status.

```bash
curl http://localhost:8000/health
```

### Query Variants

#### `GET /variants/chromosome/{chrom}` - Get Variants by Chromosome

Get all variants from a specific chromosome.

**Parameters:**
- `chrom` (path): Chromosome name (e.g., `chr1`, `1`, `chrX`)
- `limit` (query): Max records to return (default: 100, max: 10000)
- `min_pos` (query, optional): Minimum position
- `max_pos` (query, optional): Maximum position

**Examples:**

```bash
# Get first 100 variants from chr1
curl "http://localhost:8000/variants/chromosome/chr1"

# Get variants from chr1 with position range
curl "http://localhost:8000/variants/chromosome/chr1?min_pos=1000000&max_pos=2000000&limit=50"

# Get 500 variants from chrX
curl "http://localhost:8000/variants/chromosome/chrX?limit=500"
```

**Response:**
```json
{
  "total_count": 100,
  "returned_count": 100,
  "variants": [
    {
      "chrom": "chr1",
      "pos": 12345,
      "ref": "A",
      "alt": "G",
      "gene": "GENE1",
      "variant_id": "abc123...",
      "pair_id": "def456...",
      "consequence": "missense_variant",
      "clinsig": "benign",
      "gnomadg_af": "0.001"
    }
  ]
}
```

#### `GET /variants/gene/{gene}` - Get Variants by Gene

Get all variants for a specific gene.

**Parameters:**
- `gene` (path): Gene symbol (e.g., `BRCA1`)
- `limit` (query): Max records to return (default: 100)

**Examples:**

```bash
# Get BRCA1 variants
curl "http://localhost:8000/variants/gene/BRCA1"

# Get first 200 TP53 variants
curl "http://localhost:8000/variants/gene/TP53?limit=200"
```

#### `GET /variants/region` - Get Variants by Genomic Region

Get variants in a specific genomic region.

**Parameters:**
- `chrom` (query): Chromosome name
- `start` (query): Start position (inclusive)
- `end` (query): End position (inclusive)
- `limit` (query): Max records to return (default: 100)

**Examples:**

```bash
# Get variants in chr17:43000000-43100000 (BRCA1 region)
curl "http://localhost:8000/variants/region?chrom=chr17&start=43000000&end=43100000&limit=200"
```

#### `GET /variants/pathogenic` - Get Pathogenic Variants

Get pathogenic or likely pathogenic variants from ClinVar.

**Parameters:**
- `limit` (query): Max records to return (default: 100)
- `chrom` (query, optional): Filter by chromosome

**Examples:**

```bash
# Get pathogenic variants
curl "http://localhost:8000/variants/pathogenic?limit=50"

# Get pathogenic variants on chr13
curl "http://localhost:8000/variants/pathogenic?chrom=chr13&limit=100"
```

#### `GET /variants/rare` - Get Rare Variants

Get rare variants based on gnomAD allele frequency.

**Parameters:**
- `max_af` (query): Maximum allele frequency (default: 0.01)
- `limit` (query): Max records to return (default: 100)
- `chrom` (query, optional): Filter by chromosome

**Examples:**

```bash
# Get variants with AF < 1%
curl "http://localhost:8000/variants/rare"

# Get ultra-rare variants (AF < 0.001) on chr7
curl "http://localhost:8000/variants/rare?max_af=0.001&chrom=chr7&limit=200"
```

#### `GET /variants/high-impact` - Get High-Impact Variants

Get high-impact variants (stop_gained, frameshift).

**Parameters:**
- `limit` (query): Max records to return (default: 100)
- `chrom` (query, optional): Filter by chromosome

**Examples:**

```bash
# Get high-impact variants
curl "http://localhost:8000/variants/high-impact?limit=100"

# Get high-impact variants on chr22
curl "http://localhost:8000/variants/high-impact?chrom=chr22"
```

#### `GET /variants/search` - Search Variants with Multiple Filters

Advanced search with multiple optional filters.

**Parameters:**
- `chrom` (query, optional): Chromosome filter
- `gene` (query, optional): Gene symbol filter
- `consequence` (query, optional): Consequence type (contains match)
- `max_af` (query, optional): Maximum allele frequency
- `pathogenic` (query, optional): Only pathogenic variants (default: false)
- `limit` (query): Max records to return (default: 100)

**Examples:**

```bash
# Search for rare missense variants in BRCA1
curl "http://localhost:8000/variants/search?gene=BRCA1&consequence=missense&max_af=0.01&limit=50"

# Search for pathogenic variants on chr13 with AF < 0.05
curl "http://localhost:8000/variants/search?chrom=chr13&pathogenic=true&max_af=0.05"

# Search for frameshift variants
curl "http://localhost:8000/variants/search?consequence=frameshift&limit=100"
```

### Statistics

#### `GET /stats/chromosomes` - Get Chromosome Statistics

Get variant counts and statistics per chromosome.

**Example:**

```bash
curl "http://localhost:8000/stats/chromosomes"
```

**Response:**
```json
{
  "stats": {
    "chrom": ["chr1", "chr2", "chr3", ...],
    "unique_variants": [15234, 14567, 13890, ...],
    "total_pairs": [23456, 21234, 19876, ...],
    "unique_genes": [1234, 1156, 1089, ...]
  }
}
```

#### `GET /stats/genes` - Get Top Genes by Variant Count

Get top genes with most variants.

**Parameters:**
- `top_n` (query): Number of top genes to return (default: 20, max: 100)

**Example:**

```bash
# Get top 10 genes
curl "http://localhost:8000/stats/genes?top_n=10"
```

**Response:**
```json
{
  "stats": {
    "gene": ["TTN", "MUC16", "OBSCN", ...],
    "unique_variants": [5234, 4567, 3890, ...],
    "total_pairs": [7234, 6567, 5890, ...]
  }
}
```

## Performance Characteristics

All endpoints use **lazy evaluation** for optimal performance:

| Endpoint | Pattern | Performance |
|----------|---------|-------------|
| `/variants/chromosome/{chrom}` | scan → filter chrom → limit → collect | **Fast**: Only reads 1 chromosome |
| `/variants/gene/{gene}` | scan → filter gene → sort → limit → collect | **Fast**: Selective read with predicate pushdown |
| `/variants/region` | scan → filter region → sort → limit → collect | **Very Fast**: Minimal data read |
| `/variants/pathogenic` | scan → filter clinsig → select cols → limit → collect | **Fast**: Column + row pruning |
| `/variants/rare` | scan → filter AF → select cols → limit → collect | **Fast**: Efficient numeric filtering |
| `/variants/search` | scan → multiple filters → sort → limit → collect | **Fast**: All filters pushed down |
| `/stats/chromosomes` | scan → group_by → agg → collect | **Medium**: Aggregates full dataset |
| `/stats/genes` | scan → group_by → agg → sort → limit → collect | **Medium**: Aggregates full dataset |

**Benefits for 5M+ rows:**
- ✅ No upfront data loading
- ✅ Only reads necessary rows (predicate pushdown)
- ✅ Only reads necessary columns (projection pushdown)
- ✅ Query plan optimization
- ✅ 5-10x faster than eager loading

## Usage Examples

### Python with `requests`

```python
import requests

# Base URL
BASE_URL = "http://localhost:8000"

# Get pathogenic BRCA1 variants
response = requests.get(f"{BASE_URL}/variants/search", params={
    "gene": "BRCA1",
    "pathogenic": True,
    "limit": 50
})

variants = response.json()["variants"]
for variant in variants:
    print(f"{variant['chrom']}:{variant['pos']} {variant['gene']} - {variant['clinsig']}")
```

### JavaScript with `fetch`

```javascript
// Get rare variants on chr17
const response = await fetch('http://localhost:8000/variants/rare?chrom=chr17&limit=100');
const data = await response.json();

console.log(`Found ${data.total_count} rare variants`);
data.variants.forEach(v => {
    console.log(`${v.chrom}:${v.pos} ${v.gene} (AF: ${v.gnomadg_af})`);
});
```

### cURL Examples

```bash
# Get chromosome statistics
curl http://localhost:8000/stats/chromosomes | jq

# Search for high-impact rare variants
curl "http://localhost:8000/variants/search?consequence=stop_gained&max_af=0.01&limit=50" | jq

# Get variants in BRCA1 region
curl "http://localhost:8000/variants/region?chrom=chr17&start=43000000&end=43100000" | jq '.variants[] | {pos, ref, alt, gene}'
```

## Interactive API Documentation

FastAPI automatically generates interactive API documentation:

### Swagger UI (http://localhost:8000/docs)
- Try out API endpoints directly in the browser
- See request/response schemas
- Test with different parameters
- View response data

### ReDoc (http://localhost:8000/redoc)
- Alternative documentation interface
- Clean, readable format
- Detailed endpoint descriptions
- Response models

## Development

### Enable Auto-Reload

For development with automatic code reloading:

```bash
python vep_api.py --parquet results.parquet --reload
```

### Running in Production

For production deployment with multiple workers:

```bash
uvicorn vep_api:app --host 0.0.0.0 --port 8000 --workers 4
```

### CORS Configuration

To enable CORS for web applications, add this to `vep_api.py`:

```python
from fastapi.middleware.cors import CORSMiddleware

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # Configure appropriately for production
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)
```

## Error Handling

The API returns standard HTTP status codes:

- `200 OK`: Successful request
- `400 Bad Request`: Invalid parameters
- `404 Not Found`: Resource not found
- `500 Internal Server Error`: Server error
- `503 Service Unavailable`: Parquet file not available

Example error response:
```json
{
  "detail": "End position must be greater than start position"
}
```

## Limitations

- Maximum `limit` per request: 10,000 records
- Pagination: Currently returns first N results (no offset support yet)
- Aggregation endpoints load more data (may be slower for very large files)

## Tips for Best Performance

1. **Use specific filters**: More filters = less data read
2. **Limit results**: Use `limit` parameter to control response size
3. **Use column selection**: Statistics endpoints read only necessary columns
4. **Cache results**: Consider caching frequent queries at application level
5. **Regional queries**: Most efficient - only scan specific genomic region

## Example Use Cases

### Clinical Variant Analysis

```bash
# Find rare pathogenic variants in cancer genes
curl "http://localhost:8000/variants/search?gene=BRCA1&pathogenic=true&max_af=0.001"
```

### Population Genetics

```bash
# Get chromosome-level variant distribution
curl "http://localhost:8000/stats/chromosomes"
```

### Gene Prioritization

```bash
# Find genes with most variants
curl "http://localhost:8000/stats/genes?top_n=50"
```

### Region Analysis

```bash
# Analyze specific genomic region
curl "http://localhost:8000/variants/region?chrom=chr13&start=32000000&end=33000000&limit=1000"
```

## Troubleshooting

### Server won't start

```bash
# Check if port is already in use
lsof -i :8000

# Try a different port
python vep_api.py --parquet results.parquet --port 8080
```

### Slow queries

- Check Parquet file has statistics enabled (generated by `process_vep_parallel.py`)
- Use more specific filters to reduce data scanned
- Consider file size - very large files (>10GB) may need more memory

### Connection refused

- Ensure server is running
- Check firewall settings
- Verify host/port configuration

## License

Same as parent project.
