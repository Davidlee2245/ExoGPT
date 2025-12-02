# Exosome-GPT Data Directory

This directory contains **pre-formatted EV biomarker datasets** from multiple sources, automatically queried by Step 1 (EV Biomarker Finder).

## 🚀 Automated Multi-Source Querying

When you input a disease and biofluid (e.g., "melanoma, plasma"), the system **automatically searches across all data sources**:

1. **Local files first** (fast, prioritized) - Searches all subdirectories
2. **Aggregates results** - Combines evidence from all sources
3. **Ranks biomarkers** - By evidence level, statistical significance, and surface likelihood

## 📁 Directory Structure

```
data/
├── publications/      # Data extracted from published articles
│   ├── smith_2023_melanoma.csv
│   └── jones_2024_melanoma.csv
│
├── databases/         # Pre-packaged datasets from public EV databases
│   ├── vesiclepedia_melanoma.csv
│   ├── exocarta_melanoma.csv
│   ├── evpedia_melanoma.csv
│   ├── evmirna_melanoma.csv
│   └── exorbase_melanoma.csv
│
├── experiments/       # User-generated experimental data
│   └── (place your CSV files here)
│
├── disease_mapping.csv  # Disease name normalization (DOID/MeSH)
└── README.md           # This file
```

## 📊 Data Sources

### Publications (`./publications/`)
- **Extracted from published research articles**
- Each file typically represents one publication
- Includes PMID references for traceability
- Example files: `smith_2023_melanoma.csv`, `jones_2024_melanoma.csv`

### Databases (`./databases/`)
**Pre-packaged CSV files** generated from public EV databases:

- **Vesiclepedia** (`vesiclepedia_melanoma.csv`)
  - Source: http://microvesicles.org/
  - Comprehensive EV protein and RNA database
  - Contains exosomal surface markers, cargo proteins

- **ExoCarta** (`exocarta_melanoma.csv`)
  - Source: http://www.exocarta.org/
  - Exosome protein and RNA cargo database
  - ESCRT components, flotillins, heat shock proteins

- **EVpedia** (`evpedia_melanoma.csv`)
  - Source: http://evpedia.info/
  - EV protein and RNA database
  - Tetraspanins, integrins, matrix metalloproteinases

- **EVmiRNA** (`evmirna_melanoma.csv`)
  - EV miRNA database
  - Oncogenic miRNAs (miR-21, miR-155, miR-221, etc.)
  - Tumor suppressor miRNAs (miR-125b, miR-200c, let-7a)

- **exoRBase** (`exorbase_melanoma.csv`)
  - Exosome RNA database
  - mRNA, lncRNA, circRNA in exosomes
  - Stem cell markers, immune checkpoint mRNAs

### Experiments (`./experiments/`)
- **User-generated data** from your own EV experiments
- Supports: proteomics, RNA-seq, qRT-PCR, flow cytometry, ELISA, Western blot
- See format requirements below

## 📋 Data Format

All CSV/TSV files must include these **required columns**:

| Column | Description | Example |
|--------|-------------|---------|
| `disease` | Disease name | "Metastatic melanoma" |
| `biofluid` | Biofluid type | "Plasma", "Serum", "Urine" |
| `gene_symbol` | Gene symbol | "CD274", "PDCD1" |
| `analyte_type` | Type of analyte | "protein", "miRNA", "mRNA" |

**Optional columns** (recommended for better ranking):

| Column | Description | Example |
|--------|-------------|---------|
| `uniprot` | UniProt ID | "Q9NZQ7" |
| `logfc` | Log fold change | 1.25 |
| `p_value` | P-value | 0.0001 |
| `fdr` | False discovery rate | 0.001 |
| `method` | Detection method | "ELISA", "Flow cytometry", "Mass spectrometry" |
| `study_id` | Study/experiment identifier | "study_001" |
| `pmid` | PubMed ID | "PMID:31234567" |
| `notes` | Additional notes | "PD-L1 exosomal surface marker" |

## 🔍 How It Works

Step 1 automatically:

1. **Scans all subdirectories** (`publications/`, `databases/`, `experiments/`)
2. **Loads all CSV/TSV files** with required columns
3. **Filters by disease and biofluid** (case-insensitive, supports synonyms)
4. **Aggregates evidence** across all sources (counts studies, computes mean logFC, min p-value/FDR)
5. **Ranks biomarkers** by:
   - Evidence level (high ≥5 studies, moderate ≥2, low <2)
   - Statistical significance (p-value, FDR)
   - Effect size (logFC)
   - Surface likelihood (for nanobinder targets)
6. **Returns ranked list** with supporting studies from all sources

## ➕ Adding Your Own Data

### From Publications
1. Extract biomarker tables from papers
2. Format as CSV with required columns
3. Save to `./publications/your_paper_name.csv`
4. System automatically discovers it on next run

### From Databases
1. Download data from Vesiclepedia/ExoCarta/etc.
2. Normalize columns to match format
3. Save to `./databases/database_name_disease.csv`
4. System automatically loads it

### From Experiments
1. Export your results as CSV
2. Ensure required columns are present
3. Save to `./experiments/your_experiment_name.csv`
4. System automatically processes it

## 🎯 Example Query

**Input:**
- Disease: "Metastatic melanoma"
- Biofluid: "Plasma"

**System automatically:**
- Searches `publications/smith_2023_melanoma.csv`
- Searches `publications/jones_2024_melanoma.csv`
- Searches `databases/vesiclepedia_melanoma.csv`
- Searches `databases/exocarta_melanoma.csv`
- Searches `databases/evpedia_melanoma.csv`
- Searches `databases/evmirna_melanoma.csv`
- Searches `databases/exorbase_melanoma.csv`
- Searches `experiments/*.csv` (if any)

**Output:**
- Aggregated biomarkers ranked by evidence
- CD274 (PD-L1) typically ranked #1 with high evidence (multiple studies)
- All supporting studies listed with source attribution

## 📝 Notes

- **Local files are prioritized** for speed (no network delays)
- **All sources are searched** automatically (no manual selection needed)
- **Results are aggregated** across sources (same biomarker from multiple sources = stronger evidence)
- **Disease normalization** via `disease_mapping.csv` handles synonyms (e.g., "melanoma" = "Metastatic melanoma")

