# 🔧 How to Fix HPA Pipeline

## The Problem (Root Cause)

```
❌ WRONG: Using non-existent API endpoint
   https://www.proteinatlas.org/api/search_download.php

✓ RIGHT: Use HPA's actual data access methods
   https://www.proteinatlas.org/about/download (bulk download)
   OR manual search → export from website
```

## Why It Doesn't Work

The code assumes this API exists:
```python
# This endpoint DOES NOT EXIST
GET https://www.proteinatlas.org/api/search_download.php?search=Skin%20Cutaneous%20Melanoma&format=json
```

What actually happens:
- API returns HTML homepage (not JSON)
- Code detects HTML → conversion fails
- Falls back to legacy pipeline

## Quick Fix (Choose One)

### Option 1: Download Bulk HPA Data (Recommended)

**Step 1**: Download HPA pathology data
```bash
cd /home/david/.cursor/worktrees/ExoGPT/xQgUR/data
mkdir -p hpa
cd hpa

# Download bulk pathology data
wget https://www.proteinatlas.org/download/pathology.tsv.zip
unzip pathology.tsv.zip

# Should now have: ./data/hpa/pathology.tsv
```

**Step 2**: Update `_load_hpa_data()` to read bulk file
```python
# In exo_gpt/step1_ev_biomarker_finder.py, modify _load_hpa_data():

def _load_hpa_data(hpa_dir: str, disease_tissue: str, normal_tissue: str) -> Optional[pd.DataFrame]:
    """Load HPA pathology bulk data and filter by cancer type."""
    
    # Try bulk file first
    bulk_path = os.path.join(hpa_dir, "pathology.tsv")
    if os.path.exists(bulk_path):
        print(f"Loading HPA bulk pathology data from: {bulk_path}")
        df = pd.read_csv(bulk_path, sep='\t')
        
        # Filter by cancer type if needed
        # HPA pathology.tsv has columns like:
        # Gene, Cancer, High, Medium, Low, Not detected
        # You'll need to transform this to the expected format
        
        return df
    
    # Fall back to existing logic for pre-converted files
    # ... (existing code) ...
```

**Step 3**: Test
```bash
# Run Step 1 with melanoma
# Should now use HPA pipeline (not legacy)
```

---

### Option 2: Disable HPA Auto-Download (Simplest)

Just remove the broken download attempt and rely on legacy pipeline:

```python
# In exo_gpt/step1_ev_biomarker_finder.py:

def aggregate_biomarkers(...):
    # ...
    
    if use_hpa_pipeline:
        hpa_dir = os.path.join(data_dir, "hpa")
        hpa_df = _load_hpa_data(hpa_dir, disease_tissue, normal_tissue)
        
        # REMOVE THIS ENTIRE BLOCK (lines 549-615):
        # if hpa_df is None or hpa_df.empty:
        #     from exo_gpt.hpa_downloader import download_hpa_pathology
        #     ...
        
        # ADD THIS INSTEAD:
        if hpa_df is None or hpa_df.empty:
            print("⚠️  HPA data not found.")
            print("   To use HPA pipeline:")
            print("   1. Download from: https://www.proteinatlas.org/about/download")
            print("   2. Save to: ./data/hpa/pathology.tsv")
            print("   Using legacy pipeline instead...")
            
            result = _aggregate_biomarkers_legacy(disease, biofluid, data_dir)
            result["pipeline_mode"] = "legacy_fallback"
            result["fallback_reason"] = "HPA data not found. Download bulk data to enable HPA pipeline."
            return result
```

---

### Option 3: Manual Download Per Cancer Type

For each cancer type you need:

1. Go to: https://www.proteinatlas.org/humanproteome/pathology
2. Search for your cancer (e.g., "Melanoma")
3. Click "Download" → Export as TSV
4. Save to: `./data/hpa_raw/melanoma_hpa_proteins.csv`
5. Run Step 1 (converter will process it automatically)

## What to Do Next

**My recommendation: Option 1 (Bulk Download)**

Reasons:
- ✓ Download once, use for all cancer types
- ✓ ~500MB file, but contains everything
- ✓ No API dependency
- ✓ Fast local access

**Alternative: Option 2 (Disable Auto-Download)**

If you're okay with just using the legacy pipeline (which works fine), 
remove the broken download code and use existing EV databases.

## Testing the Fix

After implementing Option 1:

```bash
# 1. Verify HPA file exists
ls -lh data/hpa/pathology.tsv

# 2. Run Step 1 with melanoma
# Should see: "Pipeline Mode: 🚀 HPA Pipeline"

# 3. Check results
# Should have log2FC, p-values, TM, EVpedia columns
```

## Summary

| What | Status | Why |
|------|--------|-----|
| Current HPA download | ❌ Broken | API endpoint doesn't exist |
| Legacy pipeline | ✓ Working | Uses existing EV databases |
| Bulk HPA download | ✓ Best solution | Official HPA data access method |
| Manual download per cancer | ✓ Works | Tedious but functional |

**Root cause**: Code uses non-existent API endpoint `search_download.php`  
**Fix**: Use HPA bulk download or disable auto-download  
**Timeline**: 15-30 minutes to download and integrate bulk data

