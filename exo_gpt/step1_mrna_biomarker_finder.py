"""
Step 1 — EV mRNA Biomarker Finder
----------------------------------

Automated pipeline for EV-derived mRNA biomarker discovery.

Data Sources:
- ExoRBase/ExoRBase2: EV RNA-seq metrics (logFC, p-value)
- ExoCarta/Vesiclepedia: EV mRNA evidence
- TCGA: Tumor vs normal mRNA expression
- GEO: Disease vs control RNA-seq datasets
- PubMed: Literature evidence

Scoring Formula:
EV_mRNA_score = w1 * EV_presence + w2 * EV_logFC + w3 * Disease_logFC + w4 * Literature_hits
"""

from __future__ import annotations
import argparse
import json
import os
from dataclasses import dataclass, asdict
from typing import List, Dict, Any, Optional
import math

import numpy as np
import pandas as pd
import requests
from urllib.parse import quote

# mRNA data directory
DATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "data", "mRNA"))

# Cache directory for intermediate results
CACHE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "data", "mRNA", ".cache"))
os.makedirs(CACHE_DIR, exist_ok=True)

# Scoring weights
# Increased EV presence weight to emphasize mRNAs actually detected in EVs
W1_EV_PRESENCE = 5.0  # Weight for EV presence score (0-2 scale) - HIGHER weight for EV-detected mRNAs
W2_EV_LOGFC = 1.5     # Weight for EV logFC
W3_DISEASE_LOGFC = 1.5  # Weight for disease logFC
W4_LITERATURE = 0.5    # Weight for literature hits


@dataclass
class MRnaBiomarker:
    """mRNA biomarker evidence record"""
    gene_symbol: str
    ensembl_id: Optional[str] = None
    ev_logfc: Optional[float] = None
    tcga_logfc: Optional[float] = None
    geo_logfc: Optional[float] = None
    ev_presence_score: int = 0  # 0 = not detected, 1 = one dataset, 2 = multiple datasets
    ev_evidence_count: int = 0  # Number of supporting EV datasets
    literature_hits: int = 0
    score: float = 0.0
    data_sources: List[str] = None
    
    def __post_init__(self):
        if self.data_sources is None:
            self.data_sources = []


def query_exorbase(gene_symbol: str, disease: str) -> Dict[str, Any]:
    """
    Query ExoRBase/ExoRBase2 for EV RNA-seq data.
    
    Reads from ExoRBase2 .txt files in ./data/mRNA/databases/exorbase2/
    Calculates logFC by comparing disease samples to healthy controls.
    
    Returns dict with:
    - ev_logfc: log2 fold change in EVs (disease vs control)
    - p_value: p-value (calculated if possible)
    - evidence_count: number of datasets
    """
    result = {
        "ev_logfc": None,
        "p_value": None,
        "evidence_count": 0
    }
    
    # Map disease names to ExoRBase2 file names
    disease_mapping = {
        "colorectal cancer": "CRC",
        "crc": "CRC",
        "breast cancer": "BRCA",
        "brca": "BRCA",
        "pancreatic cancer": "PAAD",
        "paad": "PAAD",
        "liver cancer": "HCC",
        "hcc": "HCC",
        "hepatocellular carcinoma": "HCC",
        "brain cancer": "GBM",
        "gbm": "GBM",
        "glioblastoma": "GBM",
    }
    
    # Find matching ExoRBase2 file
    disease_lower = disease.lower().strip()
    project_code = None
    
    for key, code in disease_mapping.items():
        if key in disease_lower:
            project_code = code
            break
    
    if not project_code:
        return result
    
    exorbase_dir = os.path.join(DATA_DIR, "databases", "exorbase2")
    disease_file = os.path.join(exorbase_dir, f"{project_code}_longRNAs.txt")
    healthy_file = os.path.join(exorbase_dir, "Healthy_longRNAs.txt")
    
    if not os.path.exists(disease_file):
        return result
    
    try:
        # Read disease samples
        disease_df = pd.read_csv(disease_file, sep="\t", low_memory=False)
        
        # Read healthy controls
        healthy_df = None
        if os.path.exists(healthy_file):
            healthy_df = pd.read_csv(healthy_file, sep="\t", low_memory=False)
        
        # Find gene row (first column is Gene.symbol)
        disease_gene_col = disease_df.columns[0]
        gene_row = disease_df[disease_df[disease_gene_col].str.upper() == gene_symbol.upper()]
        
        if gene_row.empty:
            return result
        
        # Get expression values from disease samples (all columns except first are samples)
        disease_sample_cols = [col for col in disease_df.columns if col != disease_gene_col]
        disease_expr_values = gene_row[disease_sample_cols].values[0]
        # Convert to pandas Series to use dropna()
        disease_expr = pd.Series(disease_expr_values)
        disease_expr = pd.to_numeric(disease_expr, errors='coerce').dropna()
        
        if len(disease_expr) == 0:
            return result
        
        # Calculate logFC: disease vs healthy
        if healthy_df is not None:
            # Get gene column from healthy dataframe (should be same name)
            healthy_gene_col = healthy_df.columns[0]
            healthy_gene_row = healthy_df[healthy_df[healthy_gene_col].str.upper() == gene_symbol.upper()]
            if not healthy_gene_row.empty:
                # Get expression values from healthy samples (all columns except first are samples)
                healthy_sample_cols = [col for col in healthy_df.columns if col != healthy_gene_col]
                healthy_expr_values = healthy_gene_row[healthy_sample_cols].values[0]
                # Convert to pandas Series to use dropna()
                healthy_expr = pd.Series(healthy_expr_values)
                healthy_expr = pd.to_numeric(healthy_expr, errors='coerce').dropna()
                
                if len(healthy_expr) > 0:
                    # Calculate mean expression
                    disease_mean = disease_expr.mean()
                    healthy_mean = healthy_expr.mean()
                    
                    # Calculate log2 fold change
                    # Add pseudocount to avoid log(0)
                    if healthy_mean > 0 and disease_mean > 0:
                        logfc = np.log2((disease_mean + 0.1) / (healthy_mean + 0.1))
                        result["ev_logfc"] = float(logfc)
                        result["evidence_count"] = len(disease_expr) + len(healthy_expr)
                        
                        # Simple t-test for p-value (if scipy available)
                        try:
                            from scipy import stats
                            if len(disease_expr) > 1 and len(healthy_expr) > 1:
                                _, p_val = stats.ttest_ind(disease_expr, healthy_expr)
                                result["p_value"] = float(p_val)
                        except ImportError:
                            pass
        else:
            # If no healthy controls, use disease mean as expression value
            # (can't calculate logFC without controls)
            disease_mean = disease_expr.mean()
            if disease_mean > 0:
                result["evidence_count"] = len(disease_expr)
                # Store as expression value (not logFC)
                result["ev_logfc"] = float(np.log2(disease_mean + 0.1))
        
    except Exception as e:
        print(f"⚠️  Error reading ExoRBase2 data: {e}")
        import traceback
        traceback.print_exc()
    
    return result


def _query_ev_databases_from_preloaded(gene_symbol: str, ev_databases_data: Dict[str, pd.DataFrame]) -> Dict[str, Any]:
    """
    Query EV databases from pre-loaded DataFrames (optimized version).
    
    Uses the full Vesiclepedia database DataFrame (pre-loaded).
    
    Args:
        gene_symbol: Gene symbol to query
        ev_databases_data: Dictionary containing pre-loaded Vesiclepedia DataFrame
    
    Returns:
        Dict with ev_presence_score and ev_evidence_count
    """
    result = {
        "ev_presence_score": 0,
        "ev_evidence_count": 0
    }
    
    evidence_count = 0
    
    # Search through pre-loaded Vesiclepedia database
    for df in ev_databases_data.values():
        try:
            # Normalize column names
            df.columns = [c.strip().upper() for c in df.columns]
            
            # Find gene symbol and content type columns
            gene_symbol_col = None
            content_type_col = None
            
            for col in df.columns:
                if "GENE SYMBOL" in col or "GENE_SYMBOL" in col:
                    gene_symbol_col = col
                elif "CONTENT TYPE" in col or "CONTENT_TYPE" in col:
                    content_type_col = col
            
            if gene_symbol_col and content_type_col:
                # Filter for mRNA entries only and matching gene symbol
                gene_matches = df[df[gene_symbol_col].astype(str).str.upper() == gene_symbol.upper()]
                
                if len(gene_matches) > 0:
                    # Filter for mRNA content type only (exclude protein entries)
                    mrna_matches = gene_matches[
                        gene_matches[content_type_col].astype(str).str.lower().str.strip() == "mrna"
                    ]
                    evidence_count += len(mrna_matches)
        except Exception:
            pass
    
    result["ev_evidence_count"] = evidence_count
    
    if evidence_count == 0:
        result["ev_presence_score"] = 0
    elif evidence_count == 1:
        result["ev_presence_score"] = 1
    else:
        result["ev_presence_score"] = 2
    
    return result


def query_ev_databases(gene_symbol: str) -> Dict[str, Any]:
    """
    Query Vesiclepedia for EV mRNA evidence.
    
    Uses the full VESICLEPEDIA_PROTEIN_MRNA_DETAILS_5.1.txt file to search for mRNA entries.
    
    Returns dict with:
    - ev_presence_score: 0 (not detected), 1 (one dataset), 2 (multiple datasets)
    - ev_evidence_count: number of supporting datasets
    """
    result = {
        "ev_presence_score": 0,
        "ev_evidence_count": 0
    }
    
    # Use full Vesiclepedia database file
    vesiclepedia_file = os.path.join(DATA_DIR, "databases", "VESICLEPEDIA_PROTEIN_MRNA_DETAILS_5.1.txt")
    
    if not os.path.exists(vesiclepedia_file):
        return result
    
    evidence_count = 0
    
    try:
        # Read Vesiclepedia file (tab-separated)
        df = pd.read_csv(vesiclepedia_file, sep="\t", low_memory=False)
        
        # Normalize column names
        df.columns = [c.strip().upper() for c in df.columns]
        
        # Find gene symbol and content type columns
        gene_symbol_col = None
        content_type_col = None
        
        for col in df.columns:
            if "GENE SYMBOL" in col or "GENE_SYMBOL" in col:
                gene_symbol_col = col
            elif "CONTENT TYPE" in col or "CONTENT_TYPE" in col:
                content_type_col = col
        
        if gene_symbol_col and content_type_col:
            # Filter for mRNA entries only and matching gene symbol
            gene_matches = df[df[gene_symbol_col].astype(str).str.upper() == gene_symbol.upper()]
            
            # Filter for mRNA content type only (exclude protein entries)
            if len(gene_matches) > 0:
                mrna_matches = gene_matches[
                    gene_matches[content_type_col].astype(str).str.lower().str.strip() == "mrna"
                ]
                evidence_count = len(mrna_matches)
        else:
            # Fallback: if column names don't match expected format, try to find gene symbol
            # and assume all matches are mRNA (less accurate but better than nothing)
            for col in df.columns:
                if "SYMBOL" in col or "GENE" in col:
                    gene_symbol_col = col
                    break
            
            if gene_symbol_col:
                gene_matches = df[df[gene_symbol_col].astype(str).str.upper() == gene_symbol.upper()]
                evidence_count = len(gene_matches)
    
    except Exception as e:
        print(f"⚠️  Error reading Vesiclepedia file: {e}")
        import traceback
        traceback.print_exc()
    
    result["ev_evidence_count"] = evidence_count
    
    if evidence_count == 0:
        result["ev_presence_score"] = 0
    elif evidence_count == 1:
        result["ev_presence_score"] = 1
    else:
        result["ev_presence_score"] = 2
    
    return result


def _query_tcga_from_dataframe(gene_symbol: str, df: pd.DataFrame) -> Optional[float]:
    """
    Query TCGA data from a pre-loaded DataFrame.
    
    Args:
        gene_symbol: Gene symbol to query
        df: Pre-loaded TCGA DataFrame
    
    Returns:
        log2 fold change (tumor vs normal) or None if not found
    """
    if df is None or len(df) == 0:
        return None
    
    try:
        # Normalize column names (lowercase, strip whitespace)
        df.columns = [c.lower().strip() for c in df.columns]
        
        # Find gene symbol column (could be gene_symbol, gene, symbol, etc.)
        gene_col = None
        for col in ["gene_symbol", "gene", "symbol", "gene_name"]:
            if col in df.columns:
                gene_col = col
                break
        
        if not gene_col:
            # If no gene column found, try first column
            gene_col = df.columns[0]
        
        # Filter for gene (case-insensitive)
        gene_df = df[df[gene_col].astype(str).str.upper() == gene_symbol.upper()]
        
        if gene_df.empty:
            return None
        
        # Check if logFC column already exists
        for logfc_col in ["logfc", "log2_fc", "log_fold_change", "fold_change"]:
            if logfc_col in gene_df.columns:
                logfc_value = gene_df[logfc_col].dropna()
                if len(logfc_value) > 0:
                    # Return mean if multiple values, or first value
                    mean_logfc = logfc_value.mean()
                    if not pd.isna(mean_logfc):
                        return float(mean_logfc)
        
        # Calculate logFC from sample_type if available
        if "sample_type" in gene_df.columns:
            tumor_df = gene_df[gene_df["sample_type"].astype(str).str.contains("Tumor", case=False, na=False)]
            normal_df = gene_df[gene_df["sample_type"].astype(str).str.contains("Normal", case=False, na=False)]
            
            # Find expression column
            expr_col = None
            for col in ["log2_fpkm", "log2_expression", "expression_value", "fpkm", "expression", "log2fpkm"]:
                if col in gene_df.columns:
                    expr_col = col
                    break
            
            if expr_col:
                tumor_expr = tumor_df[expr_col].dropna()
                normal_expr = normal_df[expr_col].dropna()
                
                if len(tumor_expr) > 0 and len(normal_expr) > 0:
                    tumor_mean = tumor_expr.mean()
                    normal_mean = normal_expr.mean()
                    
                    # Calculate logFC: tumor - normal (already in log2 space)
                    if normal_mean > 0:
                        logfc = tumor_mean - normal_mean
                        return float(logfc)
                    elif tumor_mean > 0:
                        # If normal is 0, use tumor value as logFC
                        return float(tumor_mean)
                elif len(tumor_expr) > 0:
                    # Only tumor samples available
                    tumor_mean = tumor_expr.mean()
                    if not pd.isna(tumor_mean):
                        return float(tumor_mean)
        
        # Fallback: if expression column exists but no sample_type, return mean
        for expr_col in ["log2_fpkm", "log2_expression", "expression_value", "fpkm", "expression"]:
            if expr_col in gene_df.columns:
                mean_expr = gene_df[expr_col].dropna().mean()
                if not pd.isna(mean_expr):
                    return float(mean_expr)
    
    except Exception as e:
        # Silently fail for individual genes
        pass
    
    return None


def query_tcga(gene_symbol: str, disease: str) -> Optional[float]:
    """
    Query TCGA for tumor vs normal mRNA expression (logFC).
    
    Reads from TCGA CSV files in ./data/mRNA/databases/TCGA/
    
    Returns log2 fold change or None if not found.
    """
    # Check for local TCGA data files
    tcga_dir = os.path.join(DATA_DIR, "databases", "TCGA")
    
    if not os.path.exists(tcga_dir):
        return None
    
    # Look for TCGA CSV files
    import glob
    disease_safe = disease.lower().replace(" ", "_").replace("/", "_")
    
    # Try to find matching TCGA file
    tcga_files = glob.glob(os.path.join(tcga_dir, f"tcga_*_{disease_safe}.csv"))
    if not tcga_files:
        # Try any TCGA file
        tcga_files = glob.glob(os.path.join(tcga_dir, "tcga_*.csv"))
    
    if not tcga_files:
        return None
    
    # Try each file until we find the gene
    for tcga_file in tcga_files:
        try:
            # Read TCGA CSV file
            df = pd.read_csv(tcga_file, low_memory=False)
            
            # Skip if file is empty (only headers)
            if len(df) == 0:
                continue
            
            # Use optimized function
            result = _query_tcga_from_dataframe(gene_symbol, df)
            if result is not None:
                return result
        
        except Exception as e:
            print(f"⚠️  Error reading TCGA file {tcga_file}: {e}")
            continue
    
    return None


def query_geo(gene_symbol: str, disease: str) -> Optional[float]:
    """
    Query GEO for disease vs control RNA-seq comparisons.
    
    Returns log2 fold change or None if not found.
    """
    # TODO: Implement GEO API query
    # For now, return None (placeholder)
    # In production, this would:
    # - Query GEO API for relevant datasets
    # - Extract logFC from disease vs control comparisons
    return None


def query_pubmed(gene_symbol: str, disease: str) -> int:
    """
    Query PubMed for literature evidence.
    
    Searches for: "GENE + exosome + [disease]"
    
    Returns number of matching publications.
    """
    try:
        # PubMed E-utilities API
        base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        
        # Construct query: gene + exosome + disease
        query = f'"{gene_symbol}" AND exosome AND "{disease}"'
        params = {
            "db": "pubmed",
            "term": query,
            "retmode": "json",
            "retmax": 10000
        }
        
        response = requests.get(base_url, params=params, timeout=10)
        if response.status_code == 200:
            data = response.json()
            if "esearchresult" in data and "count" in data["esearchresult"]:
                count = int(data["esearchresult"]["count"])
                return count
    except Exception as e:
        print(f"⚠️  Error querying PubMed: {e}")
    
    return 0


# Gene symbol aliases/alternative names mapping
GENE_SYMBOL_ALIASES = {
    "IL8": "CXCL8",  # Interleukin 8 → C-X-C motif chemokine ligand 8
    "ALIX": "PDCD6IP",  # ALIX → Programmed cell death 6 interacting protein
    "HSA-MIR-21": None,  # MicroRNA, not a protein-coding gene
    # Add more aliases as needed
}

# Global cache for BioMart data (loaded once)
_BIOMART_LOOKUP: Optional[Dict[str, Dict[str, str]]] = None


def load_biomart_lookup() -> Dict[str, Dict[str, str]]:
    """
    Load BioMart export file and create lookup dictionaries.
    
    Returns:
        Dictionary mapping gene_symbol (upper) -> {
            'ensembl_id': 'ENSG00000...',
            'biotype': 'protein_coding' | 'lncRNA' | etc.
        }
    """
    global _BIOMART_LOOKUP
    
    if _BIOMART_LOOKUP is not None:
        return _BIOMART_LOOKUP
    
    biomart_path = os.path.join(DATA_DIR, "databases", "mart_export.txt")
    
    if not os.path.exists(biomart_path):
        print(f"⚠️  BioMart file not found at {biomart_path}")
        print("   Falling back to API queries (slower)")
        _BIOMART_LOOKUP = {}
        return _BIOMART_LOOKUP
    
    print(f"📥 Loading BioMart data from {biomart_path}...")
    lookup = {}
    
    try:
        # Read BioMart file (tab-separated)
        df = pd.read_csv(biomart_path, sep="\t", low_memory=False)
        
        # Expected columns: Gene stable ID, HGNC symbol, Gene type
        # Map to standard column names
        df.columns = [c.strip() for c in df.columns]
        
        # Find the correct column names (handle variations)
        gene_id_col = None
        gene_symbol_col = None
        biotype_col = None
        
        for col in df.columns:
            col_lower = col.lower()
            if "gene stable id" in col_lower and "version" not in col_lower:
                gene_id_col = col
            elif "hgnc symbol" in col_lower or "gene name" in col_lower or "symbol" in col_lower:
                gene_symbol_col = col
            elif "gene type" in col_lower or "biotype" in col_lower:
                biotype_col = col
        
        if not gene_id_col or not gene_symbol_col or not biotype_col:
            print(f"⚠️  Could not find required columns in BioMart file")
            print(f"   Found columns: {list(df.columns)}")
            print("   Falling back to API queries (slower)")
            _BIOMART_LOOKUP = {}
            return _BIOMART_LOOKUP
        
        # Process each row
        for _, row in df.iterrows():
            gene_symbol = str(row[gene_symbol_col]).strip()
            ensembl_id = str(row[gene_id_col]).strip()
            biotype = str(row[biotype_col]).strip()
            
            # Skip if missing data
            if pd.isna(gene_symbol) or gene_symbol == "nan" or not gene_symbol:
                continue
            if pd.isna(ensembl_id) or ensembl_id == "nan" or not ensembl_id:
                continue
            if pd.isna(biotype) or biotype == "nan" or not biotype:
                continue
            
            gene_symbol_upper = gene_symbol.upper()
            
            # If we already have this gene, prefer protein_coding entries
            if gene_symbol_upper in lookup:
                existing_biotype = lookup[gene_symbol_upper].get("biotype", "")
                # If existing is not protein_coding but new one is, update
                if existing_biotype != "protein_coding" and biotype == "protein_coding":
                    lookup[gene_symbol_upper] = {
                        "ensembl_id": ensembl_id,
                        "biotype": biotype
                    }
                # If both are protein_coding or both are non-coding, keep first (they should be same)
            else:
                lookup[gene_symbol_upper] = {
                    "ensembl_id": ensembl_id,
                    "biotype": biotype
                }
        
        print(f"✓ Loaded {len(lookup)} genes from BioMart")
        _BIOMART_LOOKUP = lookup
        return lookup
        
    except Exception as e:
        print(f"⚠️  Error loading BioMart file: {e}")
        import traceback
        traceback.print_exc()
        print("   Falling back to API queries (slower)")
        _BIOMART_LOOKUP = {}
        return _BIOMART_LOOKUP

def query_ensembl_id(gene_symbol: str, biomart_lookup: Optional[Dict[str, Dict[str, str]]] = None) -> Optional[str]:
    """
    Get Ensembl ID for a gene symbol, using BioMart lookup first, then API fallback.
    
    Args:
        gene_symbol: Gene symbol
        biomart_lookup: Optional pre-loaded BioMart lookup dictionary
    
    Returns:
        Ensembl ID or None
    """
    # Try BioMart lookup first (fast)
    if biomart_lookup is None:
        biomart_lookup = load_biomart_lookup()
    
    gene_symbol_upper = gene_symbol.upper()
    if gene_symbol_upper in biomart_lookup:
        return biomart_lookup[gene_symbol_upper].get("ensembl_id")
    
    # Try alias in BioMart
    alias = GENE_SYMBOL_ALIASES.get(gene_symbol_upper)
    if alias and alias.upper() in biomart_lookup:
        return biomart_lookup[alias.upper()].get("ensembl_id")
    
    # Fallback to API (slower)
    try:
        base_url = "https://rest.ensembl.org"
        lookup_url = f"{base_url}/lookup/symbol/homo_sapiens/{gene_symbol}"
        
        headers = {"Content-Type": "application/json"}
        response = requests.get(lookup_url, headers=headers, timeout=10)
        
        if response.status_code == 200:
            data = response.json()
            return data.get("id")
        else:
            # Try alias/alternative name if available
            if alias:
                lookup_url = f"{base_url}/lookup/symbol/homo_sapiens/{alias}"
                response = requests.get(lookup_url, headers=headers, timeout=10)
                if response.status_code == 200:
                    data = response.json()
                    return data.get("id")
    except Exception as e:
        pass  # Silently fail
    
    return None


def query_ensembl_biotype(
    gene_symbol: str, 
    ensembl_id: Optional[str] = None,
    biomart_lookup: Optional[Dict[str, Dict[str, str]]] = None
) -> Optional[str]:
    """
    Get gene biotype for a gene symbol, using BioMart lookup first, then API fallback.
    
    Args:
        gene_symbol: Gene symbol
        ensembl_id: Optional Ensembl ID (if already known)
        biomart_lookup: Optional pre-loaded BioMart lookup dictionary
    
    Returns:
        Gene biotype (e.g., "protein_coding", "lncRNA", "pseudogene") or None
    """
    # Try BioMart lookup first (fast)
    if biomart_lookup is None:
        biomart_lookup = load_biomart_lookup()
    
    gene_symbol_upper = gene_symbol.upper()
    if gene_symbol_upper in biomart_lookup:
        return biomart_lookup[gene_symbol_upper].get("biotype")
    
    # Try alias in BioMart
    alias = GENE_SYMBOL_ALIASES.get(gene_symbol_upper)
    if alias and alias.upper() in biomart_lookup:
        return biomart_lookup[alias.upper()].get("biotype")
    
    # Fallback to API (slower)
    try:
        base_url = "https://rest.ensembl.org"
        
        # If we have Ensembl ID, use it directly
        if ensembl_id:
            lookup_url = f"{base_url}/lookup/id/{ensembl_id}"
        else:
            # First get Ensembl ID from gene symbol
            lookup_url = f"{base_url}/lookup/symbol/homo_sapiens/{gene_symbol}"
        
        headers = {"Content-Type": "application/json"}
        response = requests.get(lookup_url, headers=headers, timeout=10)
        
        if response.status_code == 200:
            data = response.json()
            # Get biotype from the lookup result
            biotype = data.get("biotype")
            if biotype:
                return biotype
            
            # If biotype not in lookup, try to get it from gene info
            gene_id = data.get("id")
            if gene_id:
                # Get gene info with biotype
                info_url = f"{base_url}/overlap/id/{gene_id}"
                params = {"feature": "gene", "content-type": "application/json"}
                info_response = requests.get(info_url, params=params, headers=headers, timeout=10)
                if info_response.status_code == 200:
                    info_data = info_response.json()
                    if info_data and len(info_data) > 0:
                        return info_data[0].get("biotype")
        else:
            # Try alias if direct lookup failed
            if alias and not ensembl_id:
                lookup_url = f"{base_url}/lookup/symbol/homo_sapiens/{alias}"
                response = requests.get(lookup_url, headers=headers, timeout=10)
                if response.status_code == 200:
                    data = response.json()
                    return data.get("biotype")
    except Exception as e:
        pass  # Silently fail for individual genes
    
    return None


def filter_protein_coding_genes(
    candidate_genes: set,
    biomart_lookup: Optional[Dict[str, Dict[str, str]]] = None
) -> set:
    """
    Filter candidate genes to keep only protein-coding genes (mRNAs).
    
    Uses BioMart lookup (fast) or API fallback to get gene biotype and filters out 
    non-coding RNAs (lncRNA, antisense RNA, pseudogenes, etc.).
    
    Args:
        candidate_genes: Set of gene symbols to filter
        biomart_lookup: Optional pre-loaded BioMart lookup dictionary
    
    Returns:
        Set of protein-coding gene symbols
    """
    print(f"🔬 Filtering for protein-coding genes (mRNAs) only...")
    print(f"   Starting with {len(candidate_genes)} candidate genes")
    
    # Load BioMart lookup if not provided
    if biomart_lookup is None:
        biomart_lookup = load_biomart_lookup()
    
    if not biomart_lookup:
        print("   ⚠️  BioMart lookup not available, using API (slower)...")
    
    protein_coding_genes = set()
    processed = 0
    skipped = 0
    
    gene_list = list(candidate_genes)
    
    # Process all genes (fast lookup from BioMart)
    for gene_symbol in gene_list:
        try:
            gene_symbol_upper = gene_symbol.upper()
            
            # Try BioMart lookup first (fast)
            biotype = None
            if biomart_lookup and gene_symbol_upper in biomart_lookup:
                biotype = biomart_lookup[gene_symbol_upper].get("biotype")
            else:
                # Try alias
                alias = GENE_SYMBOL_ALIASES.get(gene_symbol_upper)
                if alias and biomart_lookup and alias.upper() in biomart_lookup:
                    biotype = biomart_lookup[alias.upper()].get("biotype")
                else:
                    # Fallback to API (slower, but only for genes not in BioMart)
                    biotype = query_ensembl_biotype(gene_symbol, biomart_lookup=biomart_lookup)
            
            if biotype == "protein_coding":
                protein_coding_genes.add(gene_symbol)
            
            processed += 1
            if processed % 5000 == 0:
                print(f"   Processed {processed}/{len(gene_list)} genes, found {len(protein_coding_genes)} protein-coding...")
        
        except Exception as e:
            # If query fails, skip this gene (conservative approach)
            skipped += 1
    
    print(f"✓ Filtered to {len(protein_coding_genes)} protein-coding genes (mRNAs)")
    print(f"   Removed {len(candidate_genes) - len(protein_coding_genes)} non-coding RNAs (lncRNA, pseudogenes, etc.)")
    if skipped > 0:
        print(f"   Skipped {skipped} genes due to errors")
    
    return protein_coding_genes


def _get_cache_key(disease: str, biofluid: str, use_exorbase2: bool, use_tcga: bool) -> str:
    """Generate cache key from parameters"""
    import hashlib
    key_str = f"{disease}_{biofluid}_{use_exorbase2}_{use_tcga}"
    return hashlib.md5(key_str.encode()).hexdigest()


def _get_cache_path(cache_key: str, cache_type: str) -> str:
    """Get cache file path"""
    return os.path.join(CACHE_DIR, f"{cache_key}_{cache_type}.json")


def _load_cache(cache_key: str, cache_type: str) -> Optional[Any]:
    """Load cached data (can be dict, list, or other JSON-serializable types)"""
    cache_path = _get_cache_path(cache_key, cache_type)
    if os.path.exists(cache_path):
        try:
            with open(cache_path, "r") as f:
                return json.load(f)
        except Exception as e:
            print(f"⚠️  Error loading cache {cache_type}: {e}")
    return None


def _save_cache(cache_key: str, cache_type: str, data: Any) -> None:
    """Save data to cache"""
    cache_path = _get_cache_path(cache_key, cache_type)
    try:
        with open(cache_path, "w") as f:
            json.dump(data, f, indent=2)
    except Exception as e:
        print(f"⚠️  Error saving cache {cache_type}: {e}")


def compute_mrna_score(
    ev_presence_score: int,
    ev_logfc: Optional[float],
    disease_logfc: Optional[float],
    literature_hits: int
) -> float:
    """
    Compute EV mRNA biomarker score.
    
    Formula:
    score = w1 * EV_presence + w2 * EV_logFC + w3 * Disease_logFC + w4 * Literature_hits
    """
    score = 0.0
    
    # EV presence (0-2 scale)
    score += W1_EV_PRESENCE * ev_presence_score
    
    # EV logFC (normalize to reasonable range, e.g., -5 to 5)
    if ev_logfc is not None:
        normalized_ev_logfc = max(-5.0, min(5.0, ev_logfc))
        score += W2_EV_LOGFC * normalized_ev_logfc
    
    # Disease logFC
    if disease_logfc is not None:
        normalized_disease_logfc = max(-5.0, min(5.0, disease_logfc))
        score += W3_DISEASE_LOGFC * normalized_disease_logfc
    
    # Literature hits (normalize to avoid excessive weight)
    normalized_literature = min(50, literature_hits)  # Cap at 50
    score += W4_LITERATURE * normalized_literature
    
    return score


def discover_mrna_biomarkers(
    disease: str,
    biofluid: str = "Plasma",
    top_n: int = 50,
    auto_download_tcga: bool = True,
    use_exorbase2: bool = True,
    use_tcga: bool = True
) -> Dict[str, Any]:
    """
    Main function to discover EV mRNA biomarkers.
    
    Args:
        disease: Disease name
        biofluid: Biofluid type (default: Plasma)
        top_n: Number of top biomarkers to return
        auto_download_tcga: If True, automatically download TCGA data if not found
    
    Returns:
        JSON-serializable dict with ranked biomarkers
    """
    print(f"🔍 Discovering EV mRNA biomarkers for {disease} in {biofluid}...")
    
    # Generate cache key
    cache_key = _get_cache_key(disease, biofluid, use_exorbase2, use_tcga)
    print(f"💾 Cache key: {cache_key[:8]}...")
    
    # Auto-download TCGA data if enabled
    if auto_download_tcga and use_tcga:
        tcga_dir = os.path.join(DATA_DIR, "databases", "TCGA")
        disease_safe = disease.lower().replace(" ", "_").replace("/", "_")
        
        import glob
        # Check if TCGA data already exists for this disease
        existing_files = glob.glob(os.path.join(tcga_dir, f"tcga_*_{disease_safe}.csv"))
        
        # Check if existing file has data (not just headers)
        has_data = False
        if existing_files:
            try:
                test_df = pd.read_csv(existing_files[0], low_memory=False)
                if len(test_df) > 0:
                    has_data = True
            except Exception:
                pass
        
        if not has_data:
            print(f"📥 Auto-downloading TCGA data for {disease}...")
            try:
                from exo_gpt.tcga_downloader import download_tcga_mrna_data
                downloaded_file = download_tcga_mrna_data(disease, tcga_dir)
                if downloaded_file:
                    print(f"✓ TCGA data downloaded successfully")
                else:
                    print(f"⚠️  TCGA download returned no file, continuing without TCGA data...")
            except Exception as e:
                print(f"⚠️  Failed to auto-download TCGA data: {e}")
                print("   Continuing without TCGA data...")
                import traceback
                traceback.print_exc()
    
    # Get candidate genes from available data sources
    candidate_genes = set()
    
    # Load from ExoRBase2 data files
    exorbase2_dir = os.path.join(DATA_DIR, "databases", "exorbase2")
    
    # Map disease to ExoRBase2 file
    disease_mapping = {
        "colorectal cancer": "CRC",
        "crc": "CRC",
        "breast cancer": "BRCA",
        "brca": "BRCA",
        "pancreatic cancer": "PAAD",
        "paad": "PAAD",
        "liver cancer": "HCC",
        "hcc": "HCC",
        "hepatocellular carcinoma": "HCC",
        "brain cancer": "GBM",
        "gbm": "GBM",
        "glioblastoma": "GBM",
    }
    
    disease_lower = disease.lower().strip()
    project_code = None
    
    for key, code in disease_mapping.items():
        if key in disease_lower:
            project_code = code
            break
    
    if use_exorbase2 and project_code and os.path.exists(exorbase2_dir):
        exorbase2_file = os.path.join(exorbase2_dir, f"{project_code}_longRNAs.txt")
        if os.path.exists(exorbase2_file):
            try:
                df = pd.read_csv(exorbase2_file, sep="\t", low_memory=False)
                # First column is Gene.symbol
                if len(df.columns) > 0:
                    gene_col = df.columns[0]
                    candidate_genes.update(df[gene_col].dropna().astype(str).str.upper().unique())
                    print(f"✓ Loaded {len(candidate_genes)} genes from ExoRBase2 ({project_code})")
            except Exception as e:
                print(f"⚠️  Error loading ExoRBase2 data: {e}")
    
    # Also check legacy ExoRBase CSV files (if they exist and ExoRBase2 is enabled)
    if use_exorbase2:
        exorbase_path = os.path.join(DATA_DIR, "databases", f"exorbase_{disease.lower().replace(' ', '_')}.csv")
        if os.path.exists(exorbase_path):
            try:
                df = pd.read_csv(exorbase_path)
                df.columns = [c.lower() for c in df.columns]
                if "gene_symbol" in df.columns:
                    # Filter for mRNA/RNA
                    if "analyte_type" in df.columns:
                        rna_df = df[df["analyte_type"].str.lower().isin(["mrna", "rna", "transcript"])]
                        candidate_genes.update(rna_df["gene_symbol"].dropna().str.upper().unique())
                    else:
                        candidate_genes.update(df["gene_symbol"].dropna().str.upper().unique())
            except Exception as e:
                print(f"⚠️  Error loading ExoRBase data: {e}")
    
    # Also check ExoCarta/Vesiclepedia for additional candidates (always include these)
    import glob
    for path in glob.glob(os.path.join(DATA_DIR, "databases", "exocarta_*.csv")):
        try:
            df = pd.read_csv(path)
            df.columns = [c.lower() for c in df.columns]
            if "gene_symbol" in df.columns:
                candidate_genes.update(df["gene_symbol"].dropna().str.upper().unique())
        except Exception:
            pass
    
    # If no candidate genes found and TCGA is enabled, we'll need to get genes from TCGA
    # For now, if no candidates found, return empty result
    
    if not candidate_genes:
        print("⚠️  No candidate genes found in data sources")
        return {
            "disease": disease,
            "biofluid": biofluid,
            "ev_mrna_biomarkers": [],
            "data_sources": []
        }
    
    print(f"✓ Found {len(candidate_genes)} candidate genes")
    
    # OPTIMIZATION 1: Pre-load all data into memory to avoid repeated file reads
    print("📦 Pre-loading data into memory...")
    exorbase2_data = None
    healthy_data = None
    tcga_data = None
    ev_databases_data = {}
    
    if use_exorbase2:
        exorbase2_dir = os.path.join(DATA_DIR, "databases", "exorbase2")
        disease_lower = disease.lower().strip()
        project_code = None
        for key, code in {
            "colorectal cancer": "CRC", "crc": "CRC",
            "breast cancer": "BRCA", "brca": "BRCA",
            "pancreatic cancer": "PAAD", "paad": "PAAD",
            "liver cancer": "HCC", "hcc": "HCC",
            "hepatocellular carcinoma": "HCC",
            "brain cancer": "GBM", "gbm": "GBM",
            "glioblastoma": "GBM",
        }.items():
            if key in disease_lower:
                project_code = code
                break
        
        if project_code:
            disease_file = os.path.join(exorbase2_dir, f"{project_code}_longRNAs.txt")
            healthy_file = os.path.join(exorbase2_dir, "Healthy_longRNAs.txt")
            if os.path.exists(disease_file):
                try:
                    exorbase2_data = pd.read_csv(disease_file, sep="\t", low_memory=False)
                    print(f"  ✓ Loaded ExoRBase2 disease data: {len(exorbase2_data)} genes")
                except Exception as e:
                    print(f"  ⚠️  Error loading ExoRBase2 disease data: {e}")
            if os.path.exists(healthy_file):
                try:
                    healthy_data = pd.read_csv(healthy_file, sep="\t", low_memory=False)
                    print(f"  ✓ Loaded ExoRBase2 healthy data: {len(healthy_data)} genes")
                except Exception as e:
                    print(f"  ⚠️  Error loading ExoRBase2 healthy data: {e}")
    
    if use_tcga:
        tcga_dir = os.path.join(DATA_DIR, "databases", "TCGA")
        import glob
        disease_safe = disease.lower().replace(" ", "_").replace("/", "_")
        tcga_files = glob.glob(os.path.join(tcga_dir, f"tcga_*_{disease_safe}.csv"))
        if not tcga_files:
            tcga_files = glob.glob(os.path.join(tcga_dir, "tcga_*.csv"))
        if tcga_files:
            try:
                tcga_data = pd.read_csv(tcga_files[0], low_memory=False)
                print(f"  ✓ Loaded TCGA data: {len(tcga_data)} records")
            except Exception as e:
                print(f"  ⚠️  Error loading TCGA data: {e}")
    
    # Pre-load Vesiclepedia database (full file, not disease-specific CSVs)
    vesiclepedia_file = os.path.join(DATA_DIR, "databases", "VESICLEPEDIA_PROTEIN_MRNA_DETAILS_5.1.txt")
    if os.path.exists(vesiclepedia_file):
        try:
            df = pd.read_csv(vesiclepedia_file, sep="\t", low_memory=False)
            ev_databases_data["vesiclepedia_full"] = df
            print(f"  ✓ Loaded Vesiclepedia database: {len(df)} entries")
        except Exception as e:
            print(f"  ⚠️  Error loading Vesiclepedia database: {e}")
    
    print("✓ Data pre-loading complete")
    
    # STEP 1: Load BioMart data and filter for protein-coding genes only
    print("\n" + "="*60)
    print("STEP 1: Filtering for protein-coding genes (mRNAs) only")
    print("="*60)
    print("   Excluding: lncRNA, antisense RNA, pseudogenes, and other non-coding RNAs")
    
    # Load BioMart lookup (fast, local file)
    biomart_lookup = load_biomart_lookup()
    
    # Step 1a: Build Ensembl ID cache from BioMart (fast, no API calls)
    print("🔍 Step 1a: Building Ensembl ID cache from BioMart...")
    ensembl_cache = {}
    if biomart_lookup:
        gene_list = list(candidate_genes)
        total_genes = len(gene_list)
        not_found_count = 0
        
        # Fast lookup: only use BioMart, skip API calls during cache building
        for i, gene_symbol in enumerate(gene_list):
            gene_symbol_upper = gene_symbol.upper()
            
            # Direct lookup
            if gene_symbol_upper in biomart_lookup:
                ensembl_id = biomart_lookup[gene_symbol_upper].get("ensembl_id")
                if ensembl_id:
                    ensembl_cache[gene_symbol] = ensembl_id
                    continue
            
            # Try alias
            alias = GENE_SYMBOL_ALIASES.get(gene_symbol_upper)
            if alias and alias.upper() in biomart_lookup:
                ensembl_id = biomart_lookup[alias.upper()].get("ensembl_id")
                if ensembl_id:
                    ensembl_cache[gene_symbol] = ensembl_id
                    continue
            
            # Not found in BioMart (will query API later if needed)
            not_found_count += 1
            
            # Progress update every 10000 genes (more frequent for large datasets)
            if (i + 1) % 10000 == 0 or (i + 1) == total_genes:
                print(f"  Processed {i + 1}/{total_genes} genes, found {len(ensembl_cache)} Ensembl IDs...")
        
        print(f"  ✓ Retrieved {len(ensembl_cache)} Ensembl IDs from BioMart")
        if not_found_count > 0:
            print(f"    ({not_found_count} genes not found in BioMart - will query API later if needed)")
    else:
        print("  ⚠️  BioMart not available, using API (slower)...")
        # Fallback to API batch querying
        gene_list = list(candidate_genes)
        batch_size = 50
        for i in range(0, min(len(gene_list), 2000), batch_size):
            batch = gene_list[i:i + batch_size]
            for gene_symbol in batch:
                ensembl_id = query_ensembl_id(gene_symbol, biomart_lookup=biomart_lookup)
                if ensembl_id:
                    ensembl_cache[gene_symbol] = ensembl_id
            if (i + batch_size) % 500 == 0:
                print(f"  Queried {min(i + batch_size, len(gene_list))}/{min(len(gene_list), 2000)} genes...")
            import time
            time.sleep(0.1)
        print(f"  ✓ Retrieved {len(ensembl_cache)} Ensembl IDs from API")
    
    # Step 1b: Filter for protein-coding only (using BioMart, fast)
    # Check cache first - use a cache key based on candidate genes to ensure consistency
    import hashlib
    candidate_genes_sorted = sorted(candidate_genes)
    candidate_genes_hash = hashlib.md5(str(candidate_genes_sorted).encode()).hexdigest()[:8]
    protein_coding_cache_file = os.path.join(CACHE_DIR, f"protein_coding_{candidate_genes_hash}.json")
    protein_coding_genes = None
    
    if os.path.exists(protein_coding_cache_file):
        try:
            with open(protein_coding_cache_file, "r") as f:
                cached_data = json.load(f)
                # Check if cache is still valid (same candidate genes)
                if isinstance(cached_data, dict) and cached_data.get("candidate_genes_hash") == candidate_genes_hash:
                    protein_coding_genes = set(cached_data.get("protein_coding_genes", []))
                    print(f"  ✓ Loaded {len(protein_coding_genes)} protein-coding genes from cache")
                    print(f"    (Skipping time-consuming filtering step)")
                elif isinstance(cached_data, list):
                    # Backward compatibility: old format was just a list
                    protein_coding_genes = set(cached_data)
                    print(f"  ✓ Loaded {len(protein_coding_genes)} protein-coding genes from cache (old format)")
        except Exception as e:
            print(f"  ⚠️  Error loading protein-coding cache: {e}")
    
    if protein_coding_genes is None:
        print("🔍 Step 1b: Filtering for protein-coding genes using BioMart...")
        print(f"    This may take a while for {len(candidate_genes)} genes...")
        protein_coding_genes = filter_protein_coding_genes(
            candidate_genes, 
            biomart_lookup=biomart_lookup
        )
        
        # Save to cache with metadata
        try:
            cache_data = {
                "candidate_genes_hash": candidate_genes_hash,
                "candidate_genes_count": len(candidate_genes),
                "protein_coding_genes": list(protein_coding_genes),
                "protein_coding_count": len(protein_coding_genes),
                "disease": disease,
                "biofluid": biofluid
            }
            with open(protein_coding_cache_file, "w") as f:
                json.dump(cache_data, f, indent=2)
            print(f"  💾 Saved {len(protein_coding_genes)} protein-coding genes to cache")
            print(f"    Cache file: {os.path.basename(protein_coding_cache_file)}")
        except Exception as e:
            print(f"  ⚠️  Error saving protein-coding cache: {e}")
    
    if not protein_coding_genes:
        print("⚠️  No protein-coding genes found after filtering")
        return {
            "disease": disease,
            "biofluid": biofluid,
            "ev_mrna_biomarkers": [],
            "data_sources": [],
            "filter_info": {
                "filtered_by_biotype": True,
                "protein_coding_only": True,
                "original_candidate_count": len(candidate_genes),
                "filtered_count": 0
            }
        }
    
    print(f"✓ Proceeding with {len(protein_coding_genes)} protein-coding genes (mRNAs)\n")
    
    # Update candidate_genes to only include protein-coding genes
    candidate_genes = protein_coding_genes
    
    # STEP 2: Calculate ExoRBase2 logFC for protein-coding genes only
    print("\n" + "="*60)
    print("STEP 2: Calculating EV logFC (ExoRBase2 disease vs healthy)")
    print("="*60)
    
    # Try to load cached ExoRBase2 results
    cached_exorbase2 = _load_cache(cache_key, "exorbase2_results")
    exorbase2_results = {}
    if cached_exorbase2:
        print("  ✓ Loaded cached ExoRBase2 results")
        exorbase2_results = cached_exorbase2
    elif use_exorbase2 and exorbase2_data is not None and healthy_data is not None:
        print("📊 Calculating ExoRBase2 logFC for protein-coding genes (vectorized)...")
        try:
            disease_gene_col = exorbase2_data.columns[0]
            healthy_gene_col = healthy_data.columns[0]
            
            # Set gene symbol as index for faster lookup
            exorbase2_data_indexed = exorbase2_data.set_index(exorbase2_data[disease_gene_col].str.upper())
            healthy_data_indexed = healthy_data.set_index(healthy_data[healthy_gene_col].str.upper())
            
            disease_sample_cols = [col for col in exorbase2_data.columns if col != disease_gene_col]
            healthy_sample_cols = [col for col in healthy_data.columns if col != healthy_gene_col]
            
            # Process only protein-coding genes
            batch_size = 1000
            for batch_start in range(0, len(candidate_genes), batch_size):
                batch_genes = list(candidate_genes)[batch_start:batch_start + batch_size]
                
                for gene_symbol in batch_genes:
                    gene_upper = gene_symbol.upper()
                    if gene_upper in exorbase2_data_indexed.index and gene_upper in healthy_data_indexed.index:
                        disease_row = exorbase2_data_indexed.loc[gene_upper]
                        healthy_row = healthy_data_indexed.loc[gene_upper]
                        
                        disease_expr = pd.to_numeric(disease_row[disease_sample_cols], errors='coerce').dropna()
                        healthy_expr = pd.to_numeric(healthy_row[healthy_sample_cols], errors='coerce').dropna()
                        
                        if len(disease_expr) > 0 and len(healthy_expr) > 0:
                            disease_mean = disease_expr.mean()
                            healthy_mean = healthy_expr.mean()
                            
                            if healthy_mean > 0 and disease_mean > 0:
                                logfc = np.log2((disease_mean + 0.1) / (healthy_mean + 0.1))
                                exorbase2_results[gene_symbol] = {
                                    "ev_logfc": float(logfc),
                                    "evidence_count": len(disease_expr) + len(healthy_expr),
                                    "p_value": None
                                }
                                
                                # Calculate p-value if scipy available
                                try:
                                    from scipy import stats
                                    if len(disease_expr) > 1 and len(healthy_expr) > 1:
                                        _, p_val = stats.ttest_ind(disease_expr, healthy_expr)
                                        exorbase2_results[gene_symbol]["p_value"] = float(p_val)
                                except ImportError:
                                    pass
                
                if (batch_start + batch_size) % 5000 == 0:
                    print(f"  Processed {min(batch_start + batch_size, len(candidate_genes))}/{len(candidate_genes)} protein-coding genes...")
            
            print(f"  ✓ Calculated ExoRBase2 logFC for {len(exorbase2_results)} protein-coding genes")
            # Save to cache
            _save_cache(cache_key, "exorbase2_results", exorbase2_results)
            print("  💾 Saved ExoRBase2 results to cache")
        except Exception as e:
            print(f"  ⚠️  Error in vectorized ExoRBase2 calculation: {e}")
            import traceback
            traceback.print_exc()
    
    # STEP 3: Batch query PubMed for protein-coding genes only
    print("\n" + "="*60)
    print("STEP 3: Querying PubMed for literature evidence")
    print("="*60)
    
    # Try to load cached PubMed results
    cached_pubmed = _load_cache(cache_key, "pubmed_cache")
    pubmed_cache = {}
    if cached_pubmed:
        print("  ✓ Loaded cached PubMed results")
        pubmed_cache = cached_pubmed
    else:
        try:
            # Query PubMed in batches (only for protein-coding genes)
            batch_size = 100
            gene_list = list(candidate_genes)  # Already filtered to protein-coding
            for batch_start in range(0, min(len(gene_list), 1000), batch_size):  # Limit to 1000 to avoid rate limits
                batch_genes = gene_list[batch_start:batch_start + batch_size]
                for gene_symbol in batch_genes:
                    literature_hits = query_pubmed(gene_symbol, disease)
                    pubmed_cache[gene_symbol] = literature_hits
                if (batch_start + batch_size) % 200 == 0:
                    print(f"  Queried {min(batch_start + batch_size, len(gene_list))} protein-coding genes...")
            print(f"  ✓ Retrieved PubMed hits for {len(pubmed_cache)} protein-coding genes")
            # Save to cache
            _save_cache(cache_key, "pubmed_cache", pubmed_cache)
            print("  💾 Saved PubMed results to cache")
        except Exception as e:
            print(f"  ⚠️  Error batch querying PubMed: {e}")
    
    # STEP 4: Compute biomarker scores for protein-coding genes
    print("\n" + "="*60)
    print("STEP 4: Computing biomarker scores and ranking")
    print("="*60)
    
    # Try to load cached biomarkers
    cached_biomarkers = _load_cache(cache_key, "biomarkers")
    biomarkers: List[MRnaBiomarker] = []
    data_sources = set()
    
    if cached_biomarkers:
        print(f"  ✓ Loaded {len(cached_biomarkers)} cached biomarkers")
        # Convert cached dicts back to MRnaBiomarker objects
        for b_dict in cached_biomarkers:
            biomarker = MRnaBiomarker(**b_dict)
            biomarkers.append(biomarker)
            data_sources.update(biomarker.data_sources)
    else:
        # Process each protein-coding gene (using pre-loaded data)
        print(f"⚡ Processing {len(candidate_genes)} protein-coding genes (using pre-loaded data)...")
        for i, gene_symbol in enumerate(candidate_genes):
            if (i + 1) % 1000 == 0:
                print(f"  Processing gene {i + 1}/{len(candidate_genes)}: {gene_symbol}")
            
            # Use pre-loaded data instead of querying each time
            exorbase_data = exorbase2_results.get(gene_symbol, {"ev_logfc": None, "p_value": None, "evidence_count": 0}) if use_exorbase2 else {"ev_logfc": None, "p_value": None, "evidence_count": 0}
            
            # Query EV databases (using pre-loaded data)
            ev_data = _query_ev_databases_from_preloaded(gene_symbol, ev_databases_data)
            
            # Query TCGA (using pre-loaded data)
            tcga_logfc = None
            if use_tcga and tcga_data is not None:
                tcga_logfc = _query_tcga_from_dataframe(gene_symbol, tcga_data)
            
            geo_logfc = None  # GEO placeholder
            
            # Use cached results
            literature_hits = pubmed_cache.get(gene_symbol, 0)
            # Get Ensembl ID from cache or BioMart lookup
            ensembl_id = ensembl_cache.get(gene_symbol)
            if not ensembl_id and biomart_lookup:
                gene_symbol_upper = gene_symbol.upper()
                if gene_symbol_upper in biomart_lookup:
                    ensembl_id = biomart_lookup[gene_symbol_upper].get("ensembl_id")
            
            # Use TCGA logFC if available, otherwise GEO
            disease_logfc = tcga_logfc if tcga_logfc is not None else geo_logfc
            
            # Compute score
            score = compute_mrna_score(
                ev_presence_score=ev_data["ev_presence_score"],
                ev_logfc=exorbase_data["ev_logfc"],
                disease_logfc=disease_logfc,
                literature_hits=literature_hits
            )
            
            # Only include if there's some evidence
            if score > 0 or exorbase_data["ev_logfc"] is not None or disease_logfc is not None:
                biomarker = MRnaBiomarker(
                    gene_symbol=gene_symbol,
                    ensembl_id=ensembl_id,
                    ev_logfc=exorbase_data["ev_logfc"],
                    tcga_logfc=tcga_logfc,
                    geo_logfc=geo_logfc,
                    ev_presence_score=ev_data["ev_presence_score"],
                    ev_evidence_count=ev_data["ev_evidence_count"],
                    literature_hits=literature_hits,
                    score=score,
                    data_sources=["ExoRBase"] if exorbase_data["ev_logfc"] is not None else []
                )
                
                if ev_data["ev_evidence_count"] > 0:
                    biomarker.data_sources.append("ExoCarta/Vesiclepedia")
                if tcga_logfc is not None:
                    biomarker.data_sources.append("TCGA")
                if geo_logfc is not None:
                    biomarker.data_sources.append("GEO")
                if literature_hits > 0:
                    biomarker.data_sources.append("PubMed")
                
                biomarkers.append(biomarker)
                data_sources.update(biomarker.data_sources)
                
                # Save progress every 1000 genes
                if (i + 1) % 1000 == 0:
                    # Save intermediate results
                    biomarkers_dict = [asdict(b) for b in biomarkers]
                    _save_cache(cache_key, "biomarkers", biomarkers_dict)
                    print(f"  💾 Saved progress: {len(biomarkers)} biomarkers processed")
    
    # Sort by score (descending)
    biomarkers.sort(key=lambda x: x.score, reverse=True)
    
    # Take top N
    top_biomarkers = biomarkers[:top_n]
    
    # Save final results to cache (always save, even if loaded from cache, to ensure cache is up to date)
    biomarkers_dict = [asdict(b) for b in biomarkers]
    _save_cache(cache_key, "biomarkers", biomarkers_dict)
    if not cached_biomarkers:
        print("  💾 Saved final biomarkers to cache")
    else:
        print("  💾 Updated biomarkers cache")
    
    print(f"\n✓ Found {len(top_biomarkers)} ranked protein-coding mRNA biomarkers")
    print(f"   (Filtered from {len(candidate_genes)} protein-coding genes)")
    
    # Convert to JSON-serializable format
    result = {
        "disease": disease,
        "biofluid": biofluid,
        "ev_mrna_biomarkers": [asdict(b) for b in top_biomarkers],
        "data_sources": list(data_sources),
        "filter_info": {
            "filtered_by_biotype": True,
            "protein_coding_only": True,
            "total_protein_coding_genes": len(candidate_genes),
            "biomarkers_found": len(biomarkers),
            "top_n_returned": len(top_biomarkers)
        }
    }
    
    return result


def biomarkers_to_markdown(biomarkers: List[Dict[str, Any]], top_n: int = 50) -> str:
    """Convert biomarker list to markdown table"""
    if not biomarkers:
        return "No biomarkers found."
    
    top = biomarkers[:top_n]
    
    lines = [
        "## EV mRNA Biomarkers",
        "",
        "| Rank | Gene | Ensembl ID | EV logFC | Disease logFC | EV Presence | Literature | Score |",
        "|------|------|------------|----------|---------------|-------------|------------|-------|"
    ]
    
    for i, b in enumerate(top, 1):
        gene = b.get("gene_symbol", "N/A")
        ensembl = b.get("ensembl_id", "N/A") or "N/A"
        ev_logfc = f"{b.get('ev_logfc', 0):.2f}" if b.get("ev_logfc") is not None else "N/A"
        disease_logfc = f"{b.get('tcga_logfc') or b.get('geo_logfc', 0):.2f}" if (b.get("tcga_logfc") or b.get("geo_logfc")) is not None else "N/A"
        ev_presence = b.get("ev_presence_score", 0)
        literature = b.get("literature_hits", 0)
        score = f"{b.get('score', 0):.2f}"
        
        lines.append(f"| {i} | {gene} | {ensembl} | {ev_logfc} | {disease_logfc} | {ev_presence} | {literature} | {score} |")
    
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="EV mRNA Biomarker Finder")
    parser.add_argument("--disease", required=True, help="Disease name")
    parser.add_argument("--biofluid", default="Plasma", help="Biofluid type")
    parser.add_argument("--out_json", help="Output JSON file")
    parser.add_argument("--out_md", help="Output markdown file")
    parser.add_argument("--top_n", type=int, default=50, help="Number of top biomarkers")
    
    args = parser.parse_args()
    
    result = discover_mrna_biomarkers(args.disease, args.biofluid, args.top_n)
    
    if args.out_json:
        with open(args.out_json, "w") as f:
            json.dump(result, f, indent=2)
        print(f"✓ Saved JSON to {args.out_json}")
    
    if args.out_md:
        md = biomarkers_to_markdown(result["ev_mrna_biomarkers"], args.top_n)
        with open(args.out_md, "w") as f:
            f.write(f"# EV mRNA Biomarkers: {args.disease} in {args.biofluid}\n\n")
            f.write(md)
        print(f"✓ Saved markdown to {args.out_md}")
    
    print(f"\n✓ Found {len(result['ev_mrna_biomarkers'])} biomarkers")


if __name__ == "__main__":
    main()

