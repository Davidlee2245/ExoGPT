# Data Directory Structure Update

## New Directory Structure

```
ExoGPT/data/
├── protein/
│   ├── databases/
│   ├── experiments/
│   ├── hpa_raw/
│   ├── publications/
│   ├── disease_mapping.csv
│   ├── structure_mapping.csv
│   └── uniprot/
│       └── uniprot.tsv
└── mRNA/
    ├── databases/
    ├── experiments/
    ├── hpa_raw/
    ├── publications/
    ├── disease_mapping.csv
    └── structure_mapping.csv
```

## Code Updates Summary

### ✅ Updated Files

1. **Protein Modules** (use `data/protein/`):
   - `exo_gpt/step1_ev_biomarker_finder.py` - Updated DATA_DIR
   - `exo_gpt/step1_ev_biomarker_finder_sub.py` - Updated DATA_DIR
   - `exo_gpt/step2_target_epitope_curator.py` - Updated DATA_DIR
   - `exo_gpt/publication_extractor.py` - Updated uniprot path to `data/protein/uniprot/`

2. **mRNA Modules** (use `data/mRNA/`):
   - `exo_gpt/step1_mrna_biomarker_finder.py` - Updated DATA_DIR
   - `exo_gpt/step2_primer_design.py` - No data directory references (uses APIs)

3. **Data Downloader** (supports both):
   - `exo_gpt/data_downloader.py` - Added `data_type` parameter ("protein" or "mRNA")
   - Updated to use `self.data_dir` based on data_type

4. **Backend API**:
   - `backend/app.py` - Updated to use `data/protein/` for Step 1 endpoints
   - Updated structure_mapping default path to `data/protein/structure_mapping.csv`
   - Updated publication endpoints to use `data/protein/publications/`

### ✅ Additional Updates Completed

1. **Publication Upload Panel** (now context-aware):
   - Updated `PublicationUploadPanel.tsx` to accept `dataType` prop
   - Updated `App.tsx` to pass `dataType="protein"` or `dataType="mrna"` based on active step
   - Panel now uses correct endpoints based on context

2. **mRNA Publication Endpoints** (added to backend):
   - `/api/mrna/step1/publication/list` - List mRNA publications
   - `/api/mrna/step1/publication/upload` - Upload mRNA publication
   - `/api/mrna/step1/publication/details/<filename>` - Get mRNA publication details
   - `/api/mrna/step1/publication/extract` - Extract biomarkers from mRNA publication
   - `/api/mrna/step1/publication/search` - RAG search across mRNA publications

### ⚠️ Optional Future Updates

1. **HPA Modules** (if used):
   - `exo_gpt/hpa_downloader.py` - Uses hardcoded `data/hpa_raw` (could be `data/protein/hpa_raw`)
   - `exo_gpt/hpa_converter.py` - Uses hardcoded `data/hpa` (could be `data/protein/hpa`)
   - Note: HPA is primarily for protein data, so current structure may be acceptable

## Usage

### Protein Modules
- Automatically use `data/protein/` directory
- No changes needed in function calls

### mRNA Modules  
- Automatically use `data/mRNA/` directory
- No changes needed in function calls

### Data Downloader
```python
# For protein data
downloader = DataDownloader(data_type="protein")
downloader.download_all_sources(disease, biofluid)

# For mRNA data
downloader = DataDownloader(data_type="mRNA")
downloader.download_all_sources(disease, biofluid)
```

## Migration Notes

- Existing code will continue to work (defaults to protein)
- New mRNA modules automatically use correct directory
- Data files should be organized into `protein/` and `mRNA/` subdirectories

