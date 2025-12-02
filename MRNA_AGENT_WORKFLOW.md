# mRNA Agent Modules - Workflow Explanation

## Overview

The mRNA Agent modules provide an automated pipeline for EV-derived mRNA biomarker discovery and primer design. This document explains the workflow and specifies which components are **REAL WORKING CODE** vs **DUMMY/PLACEHOLDER** implementations.

---

## Workflow Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                    mRNA AGENT PIPELINE                       │
└─────────────────────────────────────────────────────────────┘

STEP 1: EV mRNA Biomarker Discovery
├── Data Source Queries
│   ├── ExoRBase/ExoRBase2      [✅ REAL - Local CSV files]
│   ├── ExoCarta/Vesiclepedia   [✅ REAL - Local CSV files]
│   ├── TCGA                     [❌ DUMMY - Returns None]
│   ├── GEO                      [❌ DUMMY - Returns None]
│   └── PubMed                   [✅ REAL - Live API calls]
├── Ensembl ID Lookup            [✅ REAL - Live API calls]
├── Scoring Algorithm            [✅ REAL - Fully implemented]
└── Output JSON                  [✅ REAL - Complete]

STEP 2: Primer Design
├── Ensembl Transcript Query     [✅ REAL - Live API calls]
├── Primer3 Integration          [✅ REAL - Requires installation]
├── Primer-BLAST Off-target      [❌ DUMMY - Returns 0]
└── Output JSON                  [✅ REAL - Complete]
```

---

## STEP 1: EV mRNA Biomarker Discovery

### Workflow Steps

1. **Input**: Disease name + Biofluid type
2. **Candidate Gene Discovery**: 
   - Scans local CSV files in `data/databases/` for ExoRBase data
   - Extracts all unique gene symbols with mRNA/RNA analyte type
3. **For Each Candidate Gene**:
   - Query ExoRBase for EV logFC and p-values
   - Query ExoCarta/Vesiclepedia for EV evidence
   - Query TCGA for tumor vs normal logFC
   - Query GEO for disease vs control logFC
   - Query PubMed for literature hits
   - Query Ensembl for gene ID
4. **Scoring**: Compute composite score using weighted formula
5. **Ranking**: Sort by score (descending)
6. **Output**: Top N biomarkers in JSON format

### Implementation Status

#### ✅ **REAL WORKING CODE**

1. **ExoRBase Query** (`query_exorbase()`)
   - **Status**: ✅ REAL
   - **Implementation**: Reads local CSV files from `data/databases/exorbase_*.csv`
   - **Functionality**: 
     - Filters by gene symbol and analyte_type (mRNA/RNA)
     - Extracts logFC, p-value, evidence count
     - Aggregates multiple entries
   - **Limitation**: Requires pre-downloaded CSV files (no live API)

2. **ExoCarta/Vesiclepedia Query** (`query_ev_databases()`)
   - **Status**: ✅ REAL
   - **Implementation**: Searches local CSV files using glob patterns
   - **Functionality**:
     - Counts evidence across multiple database files
     - Calculates presence score (0/1/2) based on evidence count
   - **Limitation**: Requires pre-downloaded CSV files

3. **PubMed Query** (`query_pubmed()`)
   - **Status**: ✅ REAL - **LIVE API CALLS**
   - **Implementation**: Uses NCBI E-utilities API
   - **Functionality**:
     - Constructs query: `"GENE" AND exosome AND "DISEASE"`
     - Returns publication count
   - **Note**: This is a **real working API integration**

4. **Ensembl ID Lookup** (`query_ensembl_id()`)
   - **Status**: ✅ REAL - **LIVE API CALLS**
   - **Implementation**: Uses Ensembl REST API
   - **Functionality**: Converts gene symbol to Ensembl ID
   - **Note**: This is a **real working API integration**

5. **Scoring Algorithm** (`compute_mrna_score()`)
   - **Status**: ✅ REAL - **FULLY IMPLEMENTED**
   - **Formula**: 
     ```
     score = w1*EV_presence + w2*EV_logFC + w3*Disease_logFC + w4*Literature_hits
     ```
   - **Weights**:
     - w1 = 2.0 (EV presence, 0-2 scale)
     - w2 = 1.5 (EV logFC)
     - w3 = 1.5 (Disease logFC)
     - w4 = 0.5 (Literature hits)
   - **Normalization**: LogFC values clamped to [-5, 5], literature capped at 50

6. **Main Discovery Function** (`discover_mrna_biomarkers()`)
   - **Status**: ✅ REAL - **FULLY FUNCTIONAL**
   - **Output**: Complete JSON with ranked biomarkers

#### ❌ **DUMMY/PLACEHOLDER**

1. **TCGA Query** (`query_tcga()`)
   - **Status**: ❌ DUMMY
   - **Current**: Returns `None` always
   - **TODO**: Implement TCGA API query or local database lookup
   - **Impact**: Disease logFC from TCGA will always be missing

2. **GEO Query** (`query_geo()`)
   - **Status**: ❌ DUMMY
   - **Current**: Returns `None` always
   - **TODO**: Implement GEO API query for RNA-seq datasets
   - **Impact**: Disease logFC from GEO will always be missing

---

## STEP 2: Primer Design

### Workflow Steps

1. **Input**: Biomarkers from Step 1 (JSON)
2. **For Each Biomarker**:
   - Query Ensembl for transcript information
   - Select target transcript (dominant isoform or isoform-specific)
   - Extract sequence and exon structure
   - Design primers using Primer3
   - Check off-target hits using Primer-BLAST
3. **Output**: Primer designs in JSON format

### Implementation Status

#### ✅ **REAL WORKING CODE**

1. **Ensembl Transcript Query** (`query_ensembl_transcript()`)
   - **Status**: ✅ REAL - **LIVE API CALLS**
   - **Implementation**: Uses Ensembl REST API
   - **Functionality**:
     - Retrieves transcript IDs for a gene
     - Gets cDNA sequence
     - Gets exon structure with coordinates
     - Identifies canonical/dominant isoform
   - **Note**: This is a **real working API integration**

2. **Primer3 Integration** (`design_primers_primer3()`)
   - **Status**: ✅ REAL - **REQUIRES INSTALLATION**
   - **Implementation**: Calls `primer3_core` command-line tool via subprocess
   - **Functionality**:
     - Creates Primer3 input file with parameters
     - Runs Primer3 executable
     - Parses output to extract primer sequences, Tm, GC%, amplicon size
   - **Requirements**: 
     - Primer3 must be installed: `conda install -c bioconda primer3`
     - If not installed, function returns `None` with warning message
   - **Parameters Used**:
     - Primer size: 18-24 bp (opt: 20)
     - Tm: 58-62°C (opt: 60)
     - GC: 40-60%
     - Amplicon: 70-150 bp
     - GC clamp: 1 base

3. **Primer Design Orchestration** (`design_primers_for_biomarker()`)
   - **Status**: ✅ REAL - **FULLY FUNCTIONAL**
   - **Functionality**:
     - Selects target transcript based on strategy
     - Attempts exon-exon junction targeting
     - Falls back to full sequence if no exon structure
   - **Limitation**: Exon junction targeting is simplified (doesn't map sequence coordinates to exon boundaries)

4. **Batch Processing** (`design_primers_batch()`)
   - **Status**: ✅ REAL - **FULLY FUNCTIONAL**
   - **Output**: Complete JSON with primer designs

#### ❌ **DUMMY/PLACEHOLDER**

1. **Primer-BLAST Off-target Check** (`check_off_target_primer_blast()`)
   - **Status**: ❌ DUMMY
   - **Current**: Always returns `0` (no off-target hits)
   - **TODO**: Implement NCBI Primer-BLAST API call or automated browser interaction
   - **Impact**: Off-target screening is not performed (always shows 0 hits)

---

## Data Flow Summary

### Step 1 Flow (mRNA Biomarker Discovery)

```
User Input (Disease, Biofluid)
    ↓
Load Local CSV Files (ExoRBase, ExoCarta, Vesiclepedia)
    ↓
Extract Candidate Genes
    ↓
For Each Gene:
    ├─→ Query ExoRBase (✅ REAL - CSV files)
    ├─→ Query ExoCarta/Vesiclepedia (✅ REAL - CSV files)
    ├─→ Query TCGA (❌ DUMMY - returns None)
    ├─→ Query GEO (❌ DUMMY - returns None)
    ├─→ Query PubMed (✅ REAL - Live API)
    └─→ Query Ensembl ID (✅ REAL - Live API)
    ↓
Compute Score (✅ REAL - Full algorithm)
    ↓
Rank & Output JSON (✅ REAL - Complete)
```

### Step 2 Flow (Primer Design)

```
Step 1 JSON Input
    ↓
For Each Biomarker:
    ├─→ Query Ensembl Transcripts (✅ REAL - Live API)
    ├─→ Select Target Transcript (✅ REAL - Strategy-based)
    ├─→ Design Primers with Primer3 (✅ REAL - Requires installation)
    └─→ Check Off-targets (❌ DUMMY - returns 0)
    ↓
Output JSON (✅ REAL - Complete)
```

---

## What Works Right Now

### ✅ Fully Functional Components

1. **mRNA Biomarker Discovery** (Step 1)
   - ✅ Reads from local CSV files (ExoRBase, ExoCarta, Vesiclepedia)
   - ✅ Queries PubMed API for literature evidence
   - ✅ Queries Ensembl API for gene IDs
   - ✅ Computes scores using full algorithm
   - ✅ Ranks and outputs biomarkers

2. **Primer Design** (Step 2)
   - ✅ Queries Ensembl API for transcript sequences
   - ✅ Designs primers using Primer3 (if installed)
   - ✅ Outputs complete primer designs

### ⚠️ Partial/Placeholder Components

1. **TCGA Integration**: Not implemented (returns None)
2. **GEO Integration**: Not implemented (returns None)
3. **Primer-BLAST Off-target**: Not implemented (returns 0)

---

## Requirements for Full Functionality

### Already Satisfied
- ✅ Python dependencies (requests, pandas, numpy)
- ✅ Ensembl API access (public, no key needed)
- ✅ PubMed API access (public, no key needed)

### Needs Installation
- ⚠️ **Primer3**: `conda install -c bioconda primer3`
  - Without this, Step 2 will fail to design primers

### Needs Implementation
- ❌ TCGA API integration or local database
- ❌ GEO API integration
- ❌ Primer-BLAST API integration

---

## Testing the Modules

### Test Step 1 (mRNA Biomarker Discovery)

```bash
# Requires CSV files in data/databases/
python -m exo_gpt.step1_mrna_biomarker_finder \
  --disease "Metastatic melanoma" \
  --biofluid "Plasma" \
  --out_json test_mrna_biomarkers.json \
  --top_n 20
```

**What will work:**
- ✅ Reads ExoRBase CSV files (if present)
- ✅ Reads ExoCarta/Vesiclepedia CSV files (if present)
- ✅ Queries PubMed API (live)
- ✅ Queries Ensembl API (live)
- ✅ Computes scores
- ✅ Outputs JSON

**What won't work:**
- ❌ TCGA data (returns None)
- ❌ GEO data (returns None)

### Test Step 2 (Primer Design)

```bash
# Requires Step 1 output JSON
python -m exo_gpt.step2_primer_design \
  --biomarkers_json test_mrna_biomarkers.json \
  --out_json test_primers.json \
  --strategy dominant_isoform
```

**What will work:**
- ✅ Queries Ensembl API for transcripts (live)
- ✅ Designs primers with Primer3 (if installed)
- ✅ Outputs JSON

**What won't work:**
- ❌ Off-target screening (always returns 0)
- ⚠️ Primer3 must be installed first

---

## Summary Table

| Component | Status | Type | Notes |
|-----------|--------|------|-------|
| **Step 1: Biomarker Discovery** |
| ExoRBase Query | ✅ REAL | Local CSV | Requires pre-downloaded files |
| ExoCarta/Vesiclepedia | ✅ REAL | Local CSV | Requires pre-downloaded files |
| PubMed Query | ✅ REAL | Live API | Fully functional |
| Ensembl ID Lookup | ✅ REAL | Live API | Fully functional |
| TCGA Query | ❌ DUMMY | Placeholder | Returns None |
| GEO Query | ❌ DUMMY | Placeholder | Returns None |
| Scoring Algorithm | ✅ REAL | Full Implementation | Complete |
| **Step 2: Primer Design** |
| Ensembl Transcript Query | ✅ REAL | Live API | Fully functional |
| Primer3 Integration | ✅ REAL | External Tool | Requires installation |
| Primer-BLAST Off-target | ❌ DUMMY | Placeholder | Returns 0 |
| Exon Junction Targeting | ⚠️ PARTIAL | Simplified | Works but not fully mapped |

---

## Conclusion

**Overall Status**: The mRNA agent modules are **~70% functional** with real working code for:
- Local database queries (CSV files)
- Live API integrations (PubMed, Ensembl)
- Complete scoring and ranking algorithms
- Primer design with Primer3

**Missing/Placeholder**: 
- TCGA and GEO integrations (affects disease logFC)
- Primer-BLAST off-target screening (affects primer validation)

The modules are **production-ready for the implemented features** but will need the placeholder components implemented for full functionality as specified in the requirements.

