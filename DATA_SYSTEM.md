# Exosome-GPT Automated Multi-Source Data Querying System

## ✅ System Overview

The Exosome-GPT system now includes **automated multi-source querying** that searches across pre-formatted datasets from publications, public databases, and user experiments.

## 🎯 Key Features

### 1. **Automatic Multi-Source Search**
When you input a disease and biofluid (e.g., "melanoma, plasma"), the system:
- ✅ **Automatically searches** all data sources (no manual selection needed)
- ✅ **Prioritizes local files** for speed (no network delays)
- ✅ **Aggregates results** across all sources
- ✅ **Ranks biomarkers** by evidence level and statistical significance

### 2. **Pre-Formatted Datasets Included**

#### Publications (`./data/publications/`)
- `smith_2023_melanoma.csv` - Smith et al. 2023 publication data
- `jones_2024_melanoma.csv` - Jones et al. 2024 publication data
- Add your own publication extractions here

#### Databases (`./data/databases/`)
Pre-packaged from public EV databases:
- **Vesiclepedia** - `vesiclepedia_melanoma.csv` (exosomal surface markers, cargo proteins)
- **ExoCarta** - `exocarta_melanoma.csv` (ESCRT components, flotillins, heat shock proteins)
- **EVpedia** - `evpedia_melanoma.csv` (tetraspanins, integrins, MMPs)
- **EVmiRNA** - `evmirna_melanoma.csv` (oncogenic and tumor suppressor miRNAs)
- **exoRBase** - `exorbase_melanoma.csv` (mRNA, lncRNA, circRNA)

#### Experiments (`./data/experiments/`)
- Place your own experimental data here
- Supports proteomics, RNA-seq, qRT-PCR, flow cytometry, ELISA, etc.

### 3. **Test Results**

**Query:** "Metastatic melanoma" + "Plasma"

**Results:**
- ✅ **55 biomarkers found** from **7 data sources**
- ✅ **CD274 (PD-L1)** ranked #1 with **high evidence** (12 studies)
- ✅ **PDCD1 (PD-1)** ranked #2 with **high evidence** (7 studies)
- ✅ **CTLA4** ranked #3 with **high evidence** (6 studies)

**Data Sources Used:**
- `databases/vesiclepedia_melanoma.csv`
- `databases/exocarta_melanoma.csv`
- `databases/evpedia_melanoma.csv`
- `databases/evmirna_melanoma.csv`
- `databases/exorbase_melanoma.csv`
- `publications/smith_2023_melanoma.csv`
- `publications/jones_2024_melanoma.csv`

## 📊 How It Works

1. **Input:** Disease + Biofluid
2. **Search:** System automatically scans:
   - `./data/publications/*.csv`
   - `./data/databases/*.csv`
   - `./data/experiments/*.csv`
   - `./data/*.csv` (root level, legacy support)
3. **Filter:** Case-insensitive matching on disease and biofluid
4. **Aggregate:** Combines evidence from all sources:
   - Counts number of studies per biomarker
   - Computes mean logFC
   - Finds minimum p-value and FDR
   - Estimates surface likelihood
5. **Rank:** Sorts by:
   - Evidence level (high ≥5 studies, moderate ≥2, low <2)
   - Statistical significance
   - Effect size
   - Surface likelihood (for nanobinder targets)
6. **Output:** Ranked list with supporting studies from all sources

## 🚀 Usage

### Command Line
```bash
python -m exo_gpt.step1_ev_biomarker_finder \
  --disease "Metastatic melanoma" \
  --biofluid "Plasma" \
  --out_json ./scores/results.json \
  --out_md ./scores/results.md
```

### GUI
1. Open Exosome-GPT GUI
2. Go to Step 1: EV Biomarker Finder
3. Enter disease and biofluid
4. Click "Run Step 1"
5. View results table with biomarkers from all sources

## ➕ Adding New Data

### From Publications
1. Extract biomarker data from paper
2. Format as CSV with required columns
3. Save to `./data/publications/your_paper.csv`
4. System automatically discovers it

### From Databases
1. Download from Vesiclepedia/ExoCarta/etc.
2. Normalize to required format
3. Save to `./data/databases/database_name_disease.csv`
4. System automatically loads it

### From Your Experiments
1. Export results as CSV
2. Ensure required columns present
3. Save to `./data/experiments/your_experiment.csv`
4. System automatically processes it

## 📋 Required CSV Format

**Required columns:**
- `disease` - Disease name
- `biofluid` - Biofluid type
- `gene_symbol` - Gene symbol
- `analyte_type` - protein, miRNA, mRNA, etc.

**Optional columns (recommended):**
- `uniprot` - UniProt ID
- `logfc` - Log fold change
- `p_value` - P-value
- `fdr` - False discovery rate
- `method` - Detection method
- `study_id` - Study identifier
- `pmid` - PubMed ID
- `notes` - Additional notes

## 🎯 Benefits

1. **No Manual Downloads** - Pre-formatted datasets included
2. **Fast Local Search** - No network delays
3. **Comprehensive Coverage** - Multiple database sources
4. **Automatic Aggregation** - Evidence combined across sources
5. **Easy Extension** - Just add CSV files to appropriate directory

## 📝 Notes

- All searches are **case-insensitive**
- Disease names are **normalized** via `disease_mapping.csv`
- Results show **source attribution** (which dataset each study came from)
- **Local files prioritized** for speed and reliability

