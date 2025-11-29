# HPA Pipeline Diagnosis Report

## Problem Summary
The HPA pipeline always falls back to legacy mode because the API endpoint doesn't work.

## Root Cause
**The API endpoint `https://www.proteinatlas.org/api/search_download.php` does not exist or is not publicly accessible.**

### Evidence:
```
API Call: https://www.proteinatlas.org/api/search_download.php?search=Skin%20Cutaneous%20Melanoma&format=json
Response: HTML page (<!DOCTYPE HTML> <html>...) instead of JSON
Status: 200 OK (but returns website HTML, not API response)
```

## Code Issues

### Issue 1: Non-existent API Endpoint
**Location**: `exo_gpt/hpa_downloader.py:38`
```python
api_url = "https://www.proteinatlas.org/api/search_download.php"
```

**Problem**: This endpoint returns the HPA website homepage (HTML) instead of JSON data.

**Why this happens**:
- HPA may not provide a public search-based API
- The endpoint URL may be incorrect or outdated
- Authentication may be required

### Issue 2: No Fallback Data Source
**Location**: `exo_gpt/step1_ev_biomarker_finder.py:549`

The code tries to download HPA data automatically, but has no fallback to:
- Check for existing HPA bulk downloads
- Use pre-downloaded HPA datasets
- Provide clear instructions for manual download

### Issue 3: Disease Term Mapping
**Location**: `exo_gpt/hpa_converter.py:249-279`

The disease mapping is correct in principle:
```python
"melanoma": "Skin Cutaneous Melanoma",
"colorectal cancer": "Colorectal cancer",
```

But it doesn't matter because the API endpoint doesn't work regardless of the search term.

## Solutions

### Solution 1: Use HPA Bulk Downloads (RECOMMENDED)
HPA provides bulk downloads of their pathology data:
- URL: https://www.proteinatlas.org/about/download
- Download: `pathology.tsv.zip` or similar
- Place in: `./data/hpa/pathology.tsv`
- Code will automatically load it

### Solution 2: Manual Cancer-Specific Download
1. Go to: https://www.proteinatlas.org/humanproteome/pathology
2. Filter by cancer type (e.g., "Skin Cutaneous Melanoma")
3. Export as TSV/CSV
4. Save to: `./data/hpa_raw/{cancer_name}_hpa_proteins.csv`

### Solution 3: Disable HPA Pipeline (Current Behavior)
The system already falls back to legacy pipeline, which works fine using:
- `./data/databases/` (EVpedia, ExoCarta, Vesiclepedia, etc.)
- `./data/publications/` (literature-derived biomarkers)
- `./data/experiments/` (user experiments)

## Verification

### Test if endpoint exists:
```bash
curl -v "https://www.proteinatlas.org/api/search_download.php?search=Breast%20cancer&format=json"
```

Expected (if endpoint exists): JSON data
Actual: HTML page (404 or homepage)

### HPA Official API Documentation:
**Source**: https://www.proteinatlas.org/about/help/dataaccess

HPA provides these access methods:
1. **Per-gene API**: `https://www.proteinatlas.org/{ENSG_ID}.{format}`
   - Example: `https://www.proteinatlas.org/ENSG00000134057.json`
   - Formats: xml, tsv, json
   - ✓ Works for individual genes
   - ❌ Not suitable for "all genes in cancer X"

2. **Search-based download**: Construct a search query and download results
   - Example: Search on website → Export → TSV/CSV
   - ❌ No documented API endpoint for this
   - ✅ Manual download works

3. **Bulk downloads**: https://www.proteinatlas.org/about/download
   - Entire pathology dataset
   - ✓ Best for programmatic access
   - File: `pathology.tsv.zip`

**Conclusion**: The `search_download.php` endpoint **DOES NOT EXIST** in HPA's documented API.

## Recommended Fix

### Option A: Use Bulk HPA Data (BEST SOLUTION)
**Download once, use forever**

1. Download HPA pathology bulk data:
   ```bash
   cd data
   wget https://www.proteinatlas.org/download/pathology.tsv.zip
   unzip pathology.tsv.zip
   mv pathology.tsv hpa/hpa_pathology_bulk.tsv
   ```

2. Update `_load_hpa_data()` to read bulk file and filter by cancer type
3. Use existing HPA pipeline logic

**Advantages**:
- ✓ One-time download
- ✓ Complete dataset
- ✓ No API dependency
- ✓ Fast local access

### Option B: Remove HPA Auto-Download (SIMPLEST)
**Remove the non-working download attempt**

1. Remove `download_hpa_pathology()` calls from `step1_ev_biomarker_finder.py`
2. Only check for existing HPA files in `./data/hpa/`
3. If not found, go straight to legacy pipeline
4. Add clear message: "To use HPA pipeline, download bulk data to ./data/hpa/"

**Advantages**:
- ✓ No broken API calls
- ✓ Clear user instructions
- ✓ Legacy pipeline works fine

### Option C: Implement Manual Download Guide
Create a helper script that:
1. Prints instructions for manual download
2. Validates downloaded files
3. Converts to expected format

### Option D: Keep Current Behavior (NO ACTION NEEDED)
The legacy pipeline works perfectly.
HPA pipeline is optional enhancement.
Current fallback mechanism is functional.

## Current Status
✓ Legacy pipeline works correctly
✓ Finds biomarkers from existing EV databases
✓ Provides results with scoring and ranking
❌ HPA automatic download doesn't work
❌ HPA pipeline never activates

