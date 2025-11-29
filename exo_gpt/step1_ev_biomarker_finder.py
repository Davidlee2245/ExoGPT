"""
Step 1 — EV Biomarker Finder (HPA → UniProt → EVpedia Pipeline)
----------------------------------------------------------------

Automated pipeline that reproduces the workflow:

Human Protein Atlas (Disease vs Normal)
        ↓ Marker expression analysis
UniProt filtering for transmembrane proteins
        ↓
EVpedia filtering for proteins found in EVs

All processing is fully automatic without manual downloads.
Assumes bulk datasets from HPA, UniProt, and EVpedia are stored in `./data/`:
- `./data/hpa/` — HPA RNA/protein expression TSV files
- `./data/uniprot/` — UniProt SwissProt TSV export
- `./data/evpedia/` — EVpedia proteins TSV export
"""

from __future__ import annotations
import argparse
import json
import os
from dataclasses import dataclass, asdict
from typing import List, Dict, Any, Optional, Tuple
import math
import re

import numpy as np
import pandas as pd
from scipy import stats

DATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "data"))

# Minimal heuristic threshold for considering an extracellular domain
# "long enough" to be an EV surface marker candidate.
# PROM1-style long loops are typically >150 aa, but we use a conservative
# lower bound so shorter but still usable loops are not discarded.
MIN_EXTRACELLULAR_DOMAIN_LENGTH = 50


def _load_ev_tables(data_dir: str, search_subdirs: bool = True) -> List[Tuple[str, pd.DataFrame]]:
    """
    Load EV-like tables from data_dir and all subdirectories.
    
    This function is kept for backward compatibility with data_downloader.py.
    It loads the old format EV biomarker tables from publications/databases/experiments.
    
    Automatically searches:
    - ./data/publications/ - Publication-extracted data
    - ./data/databases/ - Pre-packaged database extracts
    - ./data/experiments/ - User experimental data
    - ./data/ - Root level files (legacy support)
    
    Returns list of (dataset_name, DataFrame) tuples.
    """
    tables: List[Tuple[str, pd.DataFrame]] = []
    if not os.path.exists(data_dir):
        return tables

    # Define search directories (prioritize subdirectories for organization)
    search_dirs = []
    if search_subdirs:
        subdirs = ["publications", "databases", "experiments"]
        for subdir in subdirs:
            subdir_path = os.path.join(data_dir, subdir)
            if os.path.exists(subdir_path):
                search_dirs.append(subdir_path)
    # Also search root data_dir
    search_dirs.append(data_dir)

    for search_dir in search_dirs:
        if not os.path.exists(search_dir):
            continue
            
        for fname in os.listdir(search_dir):
            fpath = os.path.join(search_dir, fname)
            if not os.path.isfile(fpath):
                continue
            if not (fname.lower().endswith(".csv") or fname.lower().endswith(".tsv")):
                continue

            try:
                df = pd.read_csv(fpath, sep=None, engine="python")
            except Exception:
                continue

            # Require at least these columns
            required = {"disease", "biofluid", "gene_symbol", "analyte_type"}
            if required.issubset({c.lower() for c in df.columns}):
                df.columns = [c.lower() for c in df.columns]
                # Include relative path for traceability
                rel_path = os.path.relpath(fpath, data_dir)
                tables.append((rel_path, df))
    
    return tables


@dataclass
class BiomarkerEvidence:
    """Biomarker evidence from HPA → UniProt → EVpedia pipeline."""
    gene_symbol: str
    uniprot: Optional[str]
    analyte_type: str
    log2fc: Optional[float]
    p_value: Optional[float]
    fdr: Optional[float]
    mean_expression_disease: Optional[float]
    mean_expression_normal: Optional[float]
    is_transmembrane: bool
    in_evpedia: bool
    evidence_level: str
    score: float
    supporting_studies: List[Dict[str, Any]]
    # Compatibility fields for Step 2
    protein_symbol: Optional[str] = None  # e.g., CD133 for PROM1
    num_studies: int = 1
    mean_logfc: Optional[float] = None  # Alias for log2fc
    min_p_value: Optional[float] = None  # Alias for p_value
    min_fdr: Optional[float] = None  # Alias for fdr
    surface_likelihood: Optional[float] = None  # Derived from is_transmembrane
    z_score_pathology: Optional[float] = None  # Z-score for pathology/disease tissue
    z_score_normal: Optional[float] = None  # Z-score for normal tissue
    # Publication evidence (from data/publications/publications_data.xlsx, if available)
    in_publications_excel: bool = False
    # UniProt / topology details for EV surface marker assessment
    in_uniprot: Optional[bool] = None
    tm_helix_count: Optional[int] = None
    tm_helix_regions: Optional[List[Dict[str, int]]] = None  # [{"start": int, "end": int}, ...]
    topo_domains: Optional[List[Dict[str, Any]]] = None  # [{"start": int, "end": int, "region_type": str}, ...]
    longest_extracellular_domain: Optional[int] = None
    has_long_extracellular_domain: Optional[bool] = None
    is_good_ev_surface_candidate: Optional[bool] = None


def _load_hpa_bulk_pathology(hpa_dir: str, disease_tissue: str, normal_tissue: str, data_dir: str = None) -> Optional[pd.DataFrame]:
    """
    Load HPA bulk pathology.tsv file and convert to expression format.
    
    HPA pathology.tsv format has:
    - Gene, Gene name, Cancer columns
    - Each row is a gene-cancer combination
    - High, Medium, Low, Not detected columns (counts of samples)
    
    Converts to: gene_symbol, tissue, expression_value format.
    
    Checks multiple locations:
    - ./data/hpa/pathology.tsv (preferred)
    - ./data/hpa_raw/pathology.tsv (raw data directory)
    - /home/david/.cursor-tutor/ExoGPT/data/hpa_raw/pathology.tsv (absolute path fallback)
    """
    # Try multiple locations
    possible_paths = []
    
    # Preferred location
    bulk_path = os.path.join(hpa_dir, "pathology.tsv")
    possible_paths.append(bulk_path)
    
    # hpa_raw directory
    if data_dir:
        hpa_raw_dir = os.path.join(data_dir, "hpa_raw")
        bulk_path_raw = os.path.join(hpa_raw_dir, "pathology.tsv")
        possible_paths.append(bulk_path_raw)
    
    # Absolute path fallback (for different workspace locations)
    fallback_paths = [
        "/home/david/.cursor-tutor/ExoGPT/data/hpa_raw/pathology.tsv",
        os.path.expanduser("~/.cursor-tutor/ExoGPT/data/hpa_raw/pathology.tsv"),
    ]
    possible_paths.extend(fallback_paths)
    
    # Find the first existing path
    bulk_path = None
    for path in possible_paths:
        if os.path.exists(path):
            bulk_path = path
            break
    
    if bulk_path is None:
        return None
    
    try:
        print(f"Loading HPA bulk pathology data from: {bulk_path}")
        df = pd.read_csv(bulk_path, sep="\t", low_memory=False)
        
        # Normalize column names
        df.columns = [c.strip() for c in df.columns]
        
        # Find gene column - prefer "Gene name" over "Gene" (Ensembl ID)
        gene_col = None
        for col in df.columns:
            col_lower = col.lower()
            if col_lower == "gene name":
                gene_col = col
                break
            elif col_lower in ["gene", "gene_name", "ensembl", "ensembl id"]:
                if gene_col is None:  # Use as fallback
                    gene_col = col
        
        if gene_col is None:
            print(f"⚠️  Could not find gene column in pathology.tsv")
            print(f"   Available columns: {list(df.columns)}")
            return None
        
        print(f"✓ Using gene column: '{gene_col}'")
        
        # Check if "Cancer" column exists (this is the actual format)
        if "Cancer" not in df.columns:
            print(f"⚠️  Could not find 'Cancer' column in pathology.tsv")
            print(f"   Available columns: {list(df.columns)}")
            return None
        
        # Extract disease name from disease_tissue (e.g., "colorectal cancer tissue" -> "colorectal cancer")
        disease_name = disease_tissue.lower().replace(" tissue", "").strip()
        
        # Map disease names to HPA cancer type names (case-insensitive)
        # Filter rows where Cancer column matches
        cancer_filter = df["Cancer"].astype(str).str.lower() == disease_name.lower()
        
        # Also try common variations
        if not cancer_filter.any():
            # Try mapping common disease names
            disease_mappings = {
                "colorectal": "colorectal cancer",
                "crc": "colorectal cancer",
                "melanoma": "skin cutaneous melanoma",
                "metastatic melanoma": "skin cutaneous melanoma",
                "breast": "breast cancer",
                "lung": "lung cancer",
                "prostate": "prostate cancer",
            }
            
            mapped_name = disease_mappings.get(disease_name.lower(), disease_name.lower() + " cancer")
            cancer_filter = df["Cancer"].astype(str).str.lower() == mapped_name.lower()
        
        if not cancer_filter.any():
            # Show available cancer types
            available_cancers = df["Cancer"].unique()[:10]
            print(f"⚠️  Could not find cancer type '{disease_name}' in pathology.tsv")
            print(f"   Available cancer types (sample): {list(available_cancers)}")
            return None
        
        cancer_df = df[cancer_filter].copy()
        print(f"✓ Found {len(cancer_df)} genes for cancer type: {cancer_df['Cancer'].iloc[0] if not cancer_df.empty else 'N/A'}")
        
        # --- Optional: derive gene-wise normal expression from normal_tissue.tsv ---
        normal_baseline: Dict[str, float] = {}
        if data_dir:
            try:
                hpa_raw_dir = os.path.join(data_dir, "hpa_raw")
                normal_path = os.path.join(hpa_raw_dir, "normal_tissue.tsv")
                if os.path.exists(normal_path):
                    normal_df = pd.read_csv(normal_path, sep="\t", low_memory=False)
                    # Prefer "Gene name" column (HGNC symbol)
                    norm_gene_col = None
                    for col in normal_df.columns:
                        cl = col.lower()
                        if cl == "gene name":
                            norm_gene_col = col
                            break
                        if cl in ["gene", "gene_name", "ensembl", "ensembl id"] and norm_gene_col is None:
                            norm_gene_col = col
                    level_col = None
                    for col in normal_df.columns:
                        if col.lower() == "level":
                            level_col = col
                            break
                    if norm_gene_col and level_col:
                        # Map qualitative levels to numeric scores
                        level_map = {
                            "not detected": 0.0,
                            "low": 1.0,
                            "medium": 2.0,
                            "high": 3.0,
                        }
                        scores: Dict[str, List[float]] = {}
                        for _, nrow in normal_df.iterrows():
                            g = str(nrow.get(norm_gene_col, "") or "").strip()
                            if not g:
                                continue
                            lvl = str(nrow.get(level_col, "") or "").strip().lower()
                            if lvl not in level_map:
                                continue
                            val = level_map[lvl]
                            if g not in scores:
                                scores[g] = []
                            scores[g].append(val)
                        for g, vals in scores.items():
                            if vals:
                                # Mean qualitative expression across all normal tissues / cell types
                                normal_baseline[g] = float(np.mean(vals))
            except Exception as e:
                print(f"⚠️  Failed to load normal_tissue.tsv for baseline expression: {e}")
        
        # Convert pathology counts to expression values
        # Use High count as primary signal, Medium as secondary
        # Expression value = weighted sum of High/Medium/Low counts
        converted_data = []
        
        for _, row in cancer_df.iterrows():
            gene = str(row[gene_col]).strip()
            if pd.isna(gene) or gene == "" or gene == "nan":
                continue
            
            # Get pathology counts
            high_count = int(row.get("High", 0)) if pd.notna(row.get("High")) else 0
            medium_count = int(row.get("Medium", 0)) if pd.notna(row.get("Medium")) else 0
            low_count = int(row.get("Low", 0)) if pd.notna(row.get("Low")) else 0
            not_detected = int(row.get("Not detected", 0)) if pd.notna(row.get("Not detected")) else 0
            
            total_samples = high_count + medium_count + low_count + not_detected
            if total_samples == 0:
                continue
            
            # Use ratio-based expression with weights so that:
            #   High > Medium > Low > Not detected
            #   e.g. 3H + 2M > 2H + 3M
            # Weighted signal = 3*High + 2*Medium
            # Weighted background = 1*Low + 0.25*Not detected
            # expression_value = (signal + 1) / (background + 1)
            signal = 3.0 * high_count + 2.0 * medium_count
            background = 1.0 * low_count + 0.25 * not_detected
            expression_value = (signal + 1.0) / (background + 1.0)
            
            # Only include if there is at least some signal (ratio > 1 ~ enrichment)
            if expression_value > 1.0:
                # Add disease tissue row with ratio-based expression
                converted_data.append({
                    "gene_symbol": gene,
                    "tissue": disease_tissue,
                    "expression_value": expression_value
                })
                
                # Add normal tissue row using baseline derived from normal_tissue.tsv if available.
                # This penalizes biomarkers that are also highly expressed in many normal tissues.
                baseline = normal_baseline.get(gene)
                if baseline is None:
                    # Default neutral baseline if we have no normal-tissue info for this gene
                    baseline = 1.0
                converted_data.append({
                    "gene_symbol": gene,
                    "tissue": normal_tissue,
                    "expression_value": baseline
                })
        
        if not converted_data:
            print(f"⚠️  No valid data extracted from pathology.tsv for '{disease_name}'")
            return None
        
        result_df = pd.DataFrame(converted_data)
        print(f"✓ Converted {len(result_df) // 2} genes from bulk pathology data")
        return result_df
        
    except Exception as e:
        print(f"❌ Error loading bulk pathology data: {e}")
        import traceback
        traceback.print_exc()
        return None


def _load_hpa_data(hpa_dir: str, disease_tissue: str, normal_tissue: str, data_dir: str = None) -> Optional[pd.DataFrame]:
    """
    Load HPA expression data from local CSV/TSV files.
    
    First tries to load bulk pathology.tsv file (if available).
    Then falls back to pre-converted expression format files.
    
    Expected formats:
    1. Bulk pathology.tsv: Gene, Cancer columns → converted to expression format
    2. Pre-converted: gene_symbol, tissue, expression_value columns
    
    Returns DataFrame with normalized gene names and expression values.
    """
    # Try bulk pathology.tsv first (even if hpa_dir doesn't exist, check hpa_raw)
    bulk_df = _load_hpa_bulk_pathology(hpa_dir, disease_tissue, normal_tissue, data_dir)
    if bulk_df is not None and not bulk_df.empty:
        print(f"✓ Successfully loaded bulk HPA data: {len(bulk_df)} rows")
        return bulk_df
    else:
        print(f"⚠️  Bulk HPA data loader returned: {type(bulk_df)} (empty: {bulk_df.empty if bulk_df is not None else 'N/A'})")
    
    # If hpa_dir doesn't exist, we can't load pre-converted files
    if not os.path.exists(hpa_dir):
        return None
    
    # Fall back to pre-converted expression format files
    hpa_files = []
    for fname in os.listdir(hpa_dir):
        fpath = os.path.join(hpa_dir, fname)
        if os.path.isfile(fpath) and (fname.lower().endswith(".csv") or fname.lower().endswith(".tsv")):
            # Skip bulk pathology file (already tried)
            if "pathology.tsv" in fname.lower():
                continue
            hpa_files.append(fpath)
    
    if not hpa_files:
        return None
    
    # Load all HPA files and concatenate
    dfs = []
    for fpath in hpa_files:
        try:
            # Try CSV first, then TSV
            try:
                df = pd.read_csv(fpath, sep=",", engine="python")
            except:
                df = pd.read_csv(fpath, sep="\t", engine="python")
            
            # Normalize column names
            df.columns = [c.lower().strip() for c in df.columns]
            
            # Check for required columns
            required_cols = ["gene_symbol", "tissue"]
            if not all(col in df.columns for col in required_cols):
                continue
            
            # Find expression column (could be "expression", "rna", "protein", "tpm", etc.)
            expr_cols = [c for c in df.columns if any(x in c for x in ["expression", "rna", "protein", "tpm", "fpkm", "count"])]
            if not expr_cols:
                continue
            
            expr_col = expr_cols[0]
            df["expression_value"] = pd.to_numeric(df[expr_col], errors="coerce")
            df = df.dropna(subset=["expression_value"])
            
            dfs.append(df[["gene_symbol", "tissue", "expression_value"]])
        except Exception:
            continue
    
    if not dfs:
        return None
    
    all_hpa = pd.concat(dfs, ignore_index=True)
    return all_hpa


def _perform_differential_expression(
    hpa_df: pd.DataFrame,
    disease_tissue: str,
    normal_tissue: str
) -> pd.DataFrame:
    """
    Perform differential expression analysis between disease and normal tissues.
    
    Computes:
    - mean expression in disease
    - mean expression in normal
    - log2FC
    - t-test p-value
    - FDR correction (Benjamini-Hochberg)
    
    Returns DataFrame with one row per gene.
    """
    disease_lower = disease_tissue.lower().strip()
    normal_lower = normal_tissue.lower().strip()
    
    # Filter by tissue (case-insensitive)
    disease_mask = hpa_df["tissue"].astype(str).str.lower().str.contains(disease_lower, na=False, case=False)
    normal_mask = hpa_df["tissue"].astype(str).str.lower().str.contains(normal_lower, na=False, case=False)
    
    disease_df = hpa_df[disease_mask].copy()
    normal_df = hpa_df[normal_mask].copy()
    
    if disease_df.empty or normal_df.empty:
        return pd.DataFrame()
    
    
    # Calculate global statistics for z-score computation
    all_disease_values = disease_df["expression_value"].values
    all_normal_values = normal_df["expression_value"].values
    
    disease_mean = float(np.mean(all_disease_values))
    disease_std = float(np.std(all_disease_values)) if len(all_disease_values) > 1 else 1.0
    if disease_std == 0:
        disease_std = 1.0
    
    normal_mean = float(np.mean(all_normal_values))
    normal_std = float(np.std(all_normal_values)) if len(all_normal_values) > 1 else 1.0
    if normal_std == 0:
        normal_std = 1.0
    
# Group by gene and compute statistics
    results = []
    
    for gene in hpa_df["gene_symbol"].unique():
        if pd.isna(gene):
            continue
        
        gene_disease = disease_df[disease_df["gene_symbol"] == gene]["expression_value"]
        gene_normal = normal_df[normal_df["gene_symbol"] == gene]["expression_value"]
        
        if len(gene_disease) == 0 or len(gene_normal) == 0:
            continue
        
        mean_disease = float(gene_disease.mean())
        mean_normal = float(gene_normal.mean())
        
        # Avoid log(0) or division by zero
        if mean_disease <= 0 or mean_normal <= 0:
            continue
        
        log2fc = math.log2(mean_disease / mean_normal)
        
        # Calculate z-scores: z = (value - mean) / std
        z_score_pathology = (mean_disease - disease_mean) / disease_std
        z_score_normal = (mean_normal - normal_mean) / normal_std
        
        # T-test
        try:
            t_stat, p_value = stats.ttest_ind(gene_disease, gene_normal)
            p_value = float(p_value) if not (np.isnan(p_value) or np.isinf(p_value)) else None
        except:
            p_value = None
        
        results.append({
            "gene_symbol": str(gene),
            "mean_expression_disease": mean_disease,
            "mean_expression_normal": mean_normal,
            "log2fc": log2fc,
            "z_score_pathology": z_score_pathology,
            "z_score_normal": z_score_normal,
            "p_value": p_value,
        })
    
    if not results:
        return pd.DataFrame()
    
    markers_df = pd.DataFrame(results)
    
    # FDR correction (Benjamini-Hochberg)
    if markers_df["p_value"].notna().any():
        valid_p = markers_df["p_value"].notna()
        p_vals = markers_df.loc[valid_p, "p_value"].values
        
        # Simple Benjamini-Hochberg FDR correction
        n = len(p_vals)
        sorted_indices = np.argsort(p_vals)
        sorted_p = p_vals[sorted_indices]
        p_adjusted = np.zeros(n)
        
        for i in range(n - 1, -1, -1):
            if i == n - 1:
                p_adjusted[sorted_indices[i]] = sorted_p[i]
            else:
                p_adjusted[sorted_indices[i]] = min(
                    sorted_p[i] * n / (i + 1),
                    p_adjusted[sorted_indices[i + 1]]
                )
        
        markers_df.loc[valid_p, "fdr"] = p_adjusted
    else:
        markers_df["fdr"] = pd.NA
    
    return markers_df


def _load_uniprot_data(uniprot_dir: str) -> Optional[pd.DataFrame]:
    """
    Load UniProt SwissProt TSV export.
    
    Expected columns:
    - Entry (UniProt ID)
    - Gene names (gene symbols)
    - Subcellular location (for membrane annotation)
    - Transmembrane (TRANSMEM annotations)
    - Topological domain (for TM helix predictions)
    """
    if not os.path.exists(uniprot_dir):
        print("⚠️  UniProt directory not found:", uniprot_dir)
        return None
    
    uniprot_files: List[str] = []
    for fname in os.listdir(uniprot_dir):
        fpath = os.path.join(uniprot_dir, fname)
        if not os.path.isfile(fpath):
            continue
        lower = fname.lower()
        if lower.endswith(".tsv") or lower.endswith(".csv"):
            uniprot_files.append(fpath)
    
    if not uniprot_files:
        print("⚠️  No UniProt .tsv/.csv files found in", uniprot_dir)
        return None
    
    # Prefer a file named uniprot.tsv if present, otherwise use the first one
    preferred = None
    for path in uniprot_files:
        if os.path.basename(path).lower() == "uniprot.tsv":
            preferred = path
            break
    if preferred is None:
        preferred = sorted(uniprot_files)[0]
    
    try:
        # Use C engine when possible; python engine complains about low_memory
        df = pd.read_csv(preferred, sep="\t", engine="c", low_memory=False)
    except Exception as e_tsv:
        print(f"⚠️  Failed to read UniProt TSV ({preferred}) with C engine: {e_tsv}")
        try:
            df = pd.read_csv(preferred, sep="\t", engine="python")
        except Exception as e_tsv_py:
            print(f"⚠️  Failed to read UniProt TSV ({preferred}) with python engine: {e_tsv_py}")
            try:
                df = pd.read_csv(preferred, sep=",", engine="c", low_memory=False)
            except Exception as e_csv:
                print(f"❌ Failed to read UniProt file ({preferred}) as CSV as well: {e_csv}")
                return None
    
    if df is None or df.empty:
        print(f"⚠️  UniProt file {preferred} loaded but is empty.")
        return None
    
    # Normalize column names (lowercase, underscores)
    df.columns = [c.lower().strip().replace(" ", "_") for c in df.columns]
    print(f"✓ Loaded UniProt data from {preferred} with shape {df.shape}")
    return df


def _build_uniprot_topology_index(
    uniprot_df: pd.DataFrame,
    min_extracellular_len: int = MIN_EXTRACELLULAR_DOMAIN_LENGTH,
) -> Dict[str, Dict[str, Any]]:
    """
    Build a per-gene UniProt topology index for downstream EV surface filtering.

    For each gene symbol (and UniProt Entry), extract:
    - Presence in UniProt
    - List of TRANSMEM helices (start, end)
    - List of TOPO_DOM regions with region_type:
        "Extracellular", "Cytoplasmic", or "Other"
    - Longest extracellular domain length
    - has_long_extracellular_domain (>= min_extracellular_len)
    - is_good_ev_surface_candidate: alias for has_long_extracellular_domain

    Returned mapping keys:
      - gene symbol (uppercased), e.g. "PROM1"
      - UniProt entry ID (uppercased), e.g. "Q9Y5Y6"
    """
    if uniprot_df is None or uniprot_df.empty:
        return {}

    # We assume columns are already normalized to lowercase with underscores
    cols = set(uniprot_df.columns)

    # Identify key columns
    entry_col = None
    for cand in ["entry", "uniprot", "uniprot_id", "accession"]:
        if cand in cols:
            entry_col = cand
            break

    gene_col = None
    for cand in ["gene_names", "gene_symbol", "gene", "gene_name"]:
        if cand in cols:
            gene_col = cand
            break

    tm_col = None
    topo_col = None
    protein_col = None
    for c in cols:
        cl = c.lower()
        if tm_col is None and "transmembrane" in cl:
            tm_col = c
        if topo_col is None and ("topological_domain" in cl or "topological domain" in cl or "topo_dom" in cl):
            topo_col = c
        if protein_col is None and ("protein_names" in cl or "protein name" in cl):
            protein_col = c

    if entry_col is None and gene_col is None:
        return {}

    tm_pattern = re.compile(r"TRANSMEM\s+(\d+)\.\.(\d+)", re.IGNORECASE)
    topo_pattern = re.compile(
        r"TOPO_DOM\s+(\d+)\.\.(\d+);\s*/note=\"([^\"]+)\"",
        re.IGNORECASE,
    )

    topology_index: Dict[str, Dict[str, Any]] = {}

    def _merge_record(existing: Optional[Dict[str, Any]], new: Dict[str, Any]) -> Dict[str, Any]:
        """
        If multiple UniProt entries map to the same gene, keep the one
        with the longest extracellular domain; otherwise keep the first.
        """
        if existing is None:
            return new
        old_len = existing.get("longest_extracellular_domain") or 0
        new_len = new.get("longest_extracellular_domain") or 0
        if new_len > old_len:
            return new
        return existing

    for _, row in uniprot_df.iterrows():
        entry_id = str(row[entry_col]).strip() if entry_col and pd.notna(row.get(entry_col)) else None
        tm_text = str(row[tm_col]) if tm_col and pd.notna(row.get(tm_col)) else ""
        topo_text = str(row[topo_col]) if topo_col and pd.notna(row.get(topo_col)) else ""
        protein_name = str(row[protein_col]) if protein_col and pd.notna(row.get(protein_col)) else ""

        # Parse TRANSMEM helices
        tm_helices: List[Dict[str, int]] = []
        for m in tm_pattern.finditer(tm_text):
            try:
                start = int(m.group(1))
                end = int(m.group(2))
            except (TypeError, ValueError):
                continue
            if start <= 0 or end < start:
                continue
            tm_helices.append({"start": start, "end": end})

        # Parse TOPO_DOM regions
        topo_domains: List[Dict[str, Any]] = []
        longest_extracellular = 0
        for m in topo_pattern.finditer(topo_text):
            try:
                start = int(m.group(1))
                end = int(m.group(2))
            except (TypeError, ValueError):
                continue
            if start <= 0 or end < start:
                continue
            note = m.group(3) or ""
            note_lower = note.lower()
            if "extracellular" in note_lower:
                region_type = "Extracellular"
            elif "cytoplasmic" in note_lower:
                region_type = "Cytoplasmic"
            else:
                region_type = note or "Other"

            length = end - start + 1
            if region_type == "Extracellular" and length > longest_extracellular:
                longest_extracellular = length

            topo_domains.append(
                {
                    "start": start,
                    "end": end,
                    "region_type": region_type,
                    "note": note,
                }
            )

        tm_helix_count = len(tm_helices) if tm_helices else 0
        has_long_extracellular = longest_extracellular >= min_extracellular_len if longest_extracellular > 0 else False

        record: Dict[str, Any] = {
            "entry": entry_id,
            "protein_name": protein_name or None,
            "tm_helix_count": tm_helix_count,
            "tm_helix_regions": tm_helices if tm_helices else None,
            "topo_domains": topo_domains if topo_domains else None,
            "longest_extracellular_domain": longest_extracellular if longest_extracellular > 0 else None,
            "has_long_extracellular_domain": has_long_extracellular if topo_domains else None,
            "is_good_ev_surface_candidate": has_long_extracellular if topo_domains else None,
        }

        # Index by UniProt Entry ID
        if entry_id:
            key = entry_id.strip().upper()
            topology_index[key] = _merge_record(topology_index.get(key), record)

        # Index by each gene symbol / alias
        if gene_col and pd.notna(row.get(gene_col)):
            genes_raw = str(row[gene_col])
            # Split on commas, semicolons, or whitespace to handle
            # multi-name cells like "PROM1 PROML1 MSTP061"
            for gene in re.split(r"[;,\s]+", genes_raw):
                gene = gene.strip().upper()
                if not gene:
                    continue
                topology_index[gene] = _merge_record(topology_index.get(gene), record)

    return topology_index


def _filter_transmembrane_proteins(
    markers_df: pd.DataFrame,
    uniprot_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Filter markers to keep only transmembrane/plasma membrane proteins.
    
    Criteria:
    - "Cell membrane" / "Plasma membrane" in subcellular location
    - TRANSMEM annotations
    - TM helix predictions (if available)
    """
    if uniprot_df is None or uniprot_df.empty:
        # If no UniProt data, return all markers (can't filter)
        markers_df["is_transmembrane"] = False
        return markers_df
    
    # Find UniProt ID column
    uniprot_id_col = None
    for col in ["entry", "uniprot", "uniprot_id", "entry_name"]:
        if col in uniprot_df.columns:
            uniprot_id_col = col
            break
    
    # Find gene symbol column
    gene_col = None
    for col in ["gene_names", "gene_symbol", "gene", "gene_name"]:
        if col in uniprot_df.columns:
            gene_col = col
            break
    
    if not uniprot_id_col:
        markers_df["is_transmembrane"] = False
        return markers_df
    
    # Find location/annotation columns
    location_cols = [c for c in uniprot_df.columns if any(x in c for x in ["location", "subcellular", "topology", "transmembrane", "tm"])]
    has_location_info = len(location_cols) > 0
    
    # Create mapping: gene_symbol -> is_transmembrane
    tm_map = {}
    
    for _, row in uniprot_df.iterrows():
        # Default assumption:
        # - If we have explicit location/annotation columns, infer TM from text.
        # - If we do NOT have those columns, assume every entry in this file is TM
        #   (user may have pre-filtered UniProt to only transmembrane/topology proteins).
        is_tm = not has_location_info
        
        if has_location_info:
            # Check subcellular location / topology annotations
            for loc_col in location_cols:
                if pd.notna(row.get(loc_col)):
                    loc_text = str(row[loc_col]).lower()
                    if any(x in loc_text for x in ["cell membrane", "plasma membrane", "membrane", "transmembrane"]):
                        is_tm = True
                        break
            
            # Check for TRANSMEM annotation in any column
            if not is_tm:
                for col in uniprot_df.columns:
                    if pd.notna(row.get(col)):
                        col_text = str(row[col]).lower()
                        if "transmem" in col_text or "tm helix" in col_text:
                            is_tm = True
                            break
        
        # Get gene symbol(s) and UniProt ID
        uniprot_id = str(row[uniprot_id_col]) if pd.notna(row.get(uniprot_id_col)) else None
        
        if gene_col and pd.notna(row.get(gene_col)):
            genes_raw = str(row[gene_col])
            # Split on commas, semicolons, or whitespace – handles "PROM1 PROML1 MSTP061"
            for gene in re.split(r"[;,\s]+", genes_raw):
                gene = gene.strip().upper()
                if gene:
                    tm_map[gene] = is_tm
                    # Also map by UniProt ID if available
                    if uniprot_id:
                        tm_map[uniprot_id] = is_tm
        
        # Also map by UniProt ID directly
        if uniprot_id:
            tm_map[uniprot_id] = is_tm
    
    # Apply filter to markers
    markers_df["is_transmembrane"] = markers_df["gene_symbol"].str.upper().map(
        lambda x: tm_map.get(x, False)
    )
    
    # Also check by UniProt ID if available
    if "uniprot" in markers_df.columns:
        markers_df["is_transmembrane"] = markers_df.apply(
            lambda row: tm_map.get(str(row["uniprot"]).upper(), row["is_transmembrane"]),
            axis=1
        )
    
    return markers_df


def _load_evpedia_data(evpedia_dir: str, data_dir: str) -> Optional[pd.DataFrame]:
    """
    Load EV/vesicle protein data for EV-enrichment filtering.
    
    Priority 1 (non-dummy, recommended):
      - data/databases/VESICLEPEDIA_PROTEIN_MRNA_DETAILS_5.1.txt
        Columns (tab-separated):
          CONTENT ID, CONTENT TYPE, ENTREZ GENE ID, GENE SYMBOL, SPECIES, EXPERIMENT ID, METHODS
        We use rows where:
          - CONTENT TYPE == 'protein' (case-insensitive)
          - GENE SYMBOL is non-empty
        Returned columns:
          - gene_symbol
    
    Legacy fallback (only if the master VESICLEPEDIA file is missing):
      - data/evpedia/*.csv, *.tsv
      - data/databases/*evpedia*.csv, *.tsv
    """
    # --- Priority 1: VESICLEPEDIA master file (non-dummy source) ---
    databases_dir = os.path.join(data_dir, "databases")
    master_candidates = []
    if os.path.exists(databases_dir):
        for fname in os.listdir(databases_dir):
            lower = fname.lower()
            if "vesiclepedia_protein_mrna_details" in lower and lower.endswith(".txt"):
                master_candidates.append(os.path.join(databases_dir, fname))
    
    if master_candidates:
        master_path = sorted(master_candidates)[0]
        try:
            df = pd.read_csv(master_path, sep="\t", engine="python", dtype=str)
        except Exception:
            df = None
        
        if df is not None and not df.empty:
            # Normalize column names
            df.columns = [c.strip() for c in df.columns]
            cols = {c.upper(): c for c in df.columns}
            content_col = cols.get("CONTENT TYPE") or cols.get("CONTENT_TYPE")
            gene_col = cols.get("GENE SYMBOL") or cols.get("GENE_SYMBOL")
            
            if content_col and gene_col:
                # Filter to protein content type
                mask_protein = df[content_col].astype(str).str.lower().str.strip() == "protein"
                dfp = df[mask_protein].copy()
                dfp["gene_symbol"] = dfp[gene_col].astype(str).str.strip()
                dfp = dfp[dfp["gene_symbol"] != ""]
                if not dfp.empty:
                    return dfp[["gene_symbol"]].drop_duplicates().reset_index(drop=True)
    
    # --- Legacy fallback: old EVpedia CSV/TSV files (for backward compatibility) ---
    evpedia_files = []
    
    # Check dedicated evpedia directory
    if os.path.exists(evpedia_dir):
        for fname in os.listdir(evpedia_dir):
            fpath = os.path.join(evpedia_dir, fname)
            if os.path.isfile(fpath) and (fname.lower().endswith(".csv") or fname.lower().endswith(".tsv")):
                evpedia_files.append(fpath)
    
    # Also check databases directory for EVpedia files
    if os.path.exists(databases_dir):
        for fname in os.listdir(databases_dir):
            if "evpedia" in fname.lower():
                fpath = os.path.join(databases_dir, fname)
                if os.path.isfile(fpath) and (fname.lower().endswith(".csv") or fname.lower().endswith(".tsv")):
                    evpedia_files.append(fpath)
    
    if not evpedia_files:
        return None
    
    dfs = []
    for fpath in evpedia_files:
        try:
            df = pd.read_csv(fpath, sep=None, engine="python")
            df.columns = [c.lower().strip() for c in df.columns]
            if "gene_symbol" in df.columns:
                cols_to_keep = ["gene_symbol"]
                if "uniprot" in df.columns:
                    cols_to_keep.append("uniprot")
                dfs.append(df[cols_to_keep].drop_duplicates())
        except Exception:
            continue
    
    if not dfs:
        return None
    
    evpedia_df = pd.concat(dfs, ignore_index=True).drop_duplicates()
    return evpedia_df


def _filter_evpedia_proteins(
    markers_df: pd.DataFrame,
    evpedia_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Filter markers to keep only proteins found in EVpedia.
    
    Inner-joins the UniProt-filtered list with EVpedia.
    """
    if evpedia_df is None or evpedia_df.empty:
        # If no EVpedia data, return all markers (can't filter)
        markers_df["in_evpedia"] = False
        return markers_df
    
    # Create set of EVpedia genes/UniProt IDs
    evpedia_genes = set(evpedia_df["gene_symbol"].str.upper().dropna().unique())
    evpedia_uniprots = set()
    if "uniprot" in evpedia_df.columns:
        evpedia_uniprots = set(evpedia_df["uniprot"].str.upper().dropna().unique())
    
    # Mark which markers are in EVpedia
    markers_df["in_evpedia"] = (
        markers_df["gene_symbol"].str.upper().isin(evpedia_genes) |
        (markers_df.get("uniprot", pd.Series()).astype(str).str.upper().isin(evpedia_uniprots) if "uniprot" in markers_df.columns else False)
    )
    
    return markers_df


def _load_publication_biomarker_index(
    data_dir: str,
) -> Optional[Dict[str, set]]:
    """
    Build an index of (disease, biofluid, gene_symbol) triples that appear
    in data/publications/publications_data.xlsx (Biomarkers sheet).
    
    Returns:
        {
          "by_disease_biofluid": { (disease_lower, biofluid_lower): {GENE1, GENE2, ...} }
        }
    or None if Excel is missing or unreadable.
    """
    publications_dir = os.path.join(data_dir, "publications")
    excel_path = os.path.join(publications_dir, "publications_data.xlsx")
    
    if not os.path.exists(excel_path):
        return None
    
    try:
        bio_df = pd.read_excel(excel_path, sheet_name="Biomarkers")
    except Exception:
        return None
    
    if bio_df is None or bio_df.empty:
        return None
    
    # Normalize columns if present
    cols = {c.lower(): c for c in bio_df.columns}
    disease_col = cols.get("disease")
    biofluid_col = cols.get("biofluid")
    gene_col = cols.get("gene_symbol") or cols.get("gene")
    reported_col = cols.get("reported_marker")  # from publication_extractor
    gene_norm_col = cols.get("gene_symbol_normalized")
    
    if not (gene_col or reported_col or gene_norm_col):
        return None
    
    index: Dict[str, set] = {}
    key_map: Dict[tuple, set] = {}
    
    for _, row in bio_df.iterrows():
        # Collect all possible identifiers for this biomarker:
        # - gene symbol (PROM1)
        # - normalized gene symbol
        # - reported marker (e.g., CD133, PD-L1) with and without spaces/hyphens
        names: set[str] = set()
        
        if gene_col:
            g = str(row.get(gene_col, "") or "").strip().upper()
            if g:
                names.add(g)
        
        if gene_norm_col:
            g_norm = str(row.get(gene_norm_col, "") or "").strip().upper()
            if g_norm:
                names.add(g_norm)
        
        if reported_col:
            rep = str(row.get(reported_col, "") or "").strip().upper()
            if rep:
                names.add(rep)
                # Also add versions without spaces or hyphens, so "CD 133" or "PD-L1"
                # will match "CD133" / "PDL1" in downstream queries.
                rep_compact = rep.replace(" ", "").replace("-", "")
                if rep_compact:
                    names.add(rep_compact)
        
        if not names:
            continue
        disease = str(row.get(disease_col, "") or "").strip().lower() if disease_col else ""
        biofluid = str(row.get(biofluid_col, "") or "").strip().lower() if biofluid_col else ""
        key = (disease, biofluid)
        if key not in key_map:
            key_map[key] = set()
        key_map[key].update(names)
    
    return {"by_disease_biofluid": key_map}


def _aggregate_biomarkers_legacy(
    disease: str,
    biofluid: str,
    data_dir: str = DATA_DIR,
) -> Dict[str, Any]:
    """
    Legacy aggregation function: searches existing EV biomarker tables.
    
    This is a fallback when HPA data is not available.
    Uses the old workflow that searches publications/databases/experiments.
    """
    # Import and call the legacy aggregate_biomarkers from the sub module
    try:
        from exo_gpt.step1_ev_biomarker_finder_sub import aggregate_biomarkers as legacy_aggregate
        result = legacy_aggregate(disease, biofluid, data_dir, search_subdirs=True)
        result["pipeline_mode"] = "legacy"
        return result
    except ImportError:
        # If submodule doesn't exist, return empty result
        return {
            "disease": disease,
            "biofluid": biofluid,
            "ev_biomarkers": [],
            "data_sources": [],
            "pipeline_mode": "legacy",
            "error": "Legacy module not available"
        }


def aggregate_biomarkers(
    disease: str,
    biofluid: str,
    data_dir: str = DATA_DIR,
    disease_tissue: Optional[str] = None,
    normal_tissue: Optional[str] = None,
    use_hpa_pipeline: bool = True,
) -> Dict[str, Any]:
    """
    Main aggregation function: HPA → UniProt → EVpedia pipeline with fallback.
    
    Args:
        disease: Disease name (e.g., "CRC")
        biofluid: Biofluid/tissue type (e.g., "CRC tissue" or "Plasma")
        data_dir: Data directory root
        disease_tissue: Optional explicit disease tissue name (defaults to disease + " tissue")
        normal_tissue: Optional explicit normal tissue name (defaults to "Normal tissue")
        use_hpa_pipeline: If True, try HPA pipeline first; if False or HPA data missing, use legacy
    
    Returns:
        JSON-serializable dict with ranked biomarkers.
    """
    # Determine tissue names
    if disease_tissue is None:
        disease_tissue = f"{disease} tissue"
    if normal_tissue is None:
        normal_tissue = "Normal tissue"
    
    # Step 1: Load HPA data and perform differential expression
    if use_hpa_pipeline:
        hpa_dir = os.path.join(data_dir, "hpa")
        hpa_df = _load_hpa_data(hpa_dir, disease_tissue, normal_tissue, data_dir)
        
        # If HPA data not found, provide instructions for manual download
        hpa_download_error = None
        hpa_conversion_error = None
        
        if hpa_df is None or hpa_df.empty:
            # Check both locations and provide helpful message
            bulk_path1 = os.path.join(hpa_dir, "pathology.tsv")
            bulk_path2 = os.path.join(data_dir, "hpa_raw", "pathology.tsv")
            
            print("⚠️  HPA data not found.")
            print(f"   Checked: {bulk_path1}")
            print(f"   Checked: {bulk_path2}")
            print("   To enable HPA pipeline:")
            print("   1. Download bulk pathology data from: https://www.proteinatlas.org/about/download")
            print("   2. Look for 'Pathology' section and download 'pathology.tsv.zip'")
            print("   3. Extract and save as: ./data/hpa/pathology.tsv OR ./data/hpa_raw/pathology.tsv")
            print("   Using legacy pipeline instead...")
            
            hpa_download_error = (
                f"HPA bulk data not found. Checked:\n"
                f"  - {bulk_path1}\n"
                f"  - {bulk_path2}\n"
                f"To use HPA pipeline:\n"
                f"1. Download pathology.tsv from https://www.proteinatlas.org/about/download\n"
                f"2. Save to: ./data/hpa/pathology.tsv OR ./data/hpa_raw/pathology.tsv\n"
                f"3. The pipeline will automatically convert and use it."
            )
            
            # If still no data, fallback to legacy
            if hpa_df is None or hpa_df.empty:
                result = _aggregate_biomarkers_legacy(disease, biofluid, data_dir)
                result["pipeline_mode"] = "legacy_fallback"
                result["fallback_reason"] = (
                    f"HPA bulk pathology data not found. Checked locations:\n"
                    f"  - {bulk_path1}\n"
                    f"  - {bulk_path2}\n"
                    f"Using existing EV biomarker data files.\n"
                    f"To enable HPA pipeline: Download pathology.tsv from https://www.proteinatlas.org/about/download "
                    f"and save to ./data/hpa/pathology.tsv OR ./data/hpa_raw/pathology.tsv"
                )
                
                if hpa_download_error:
                    result["hpa_download_error"] = hpa_download_error
                if hpa_conversion_error:
                    result["hpa_conversion_error"] = hpa_conversion_error
                
                return result
    
    markers_df = _perform_differential_expression(hpa_df, disease_tissue, normal_tissue)
    
    if markers_df.empty:
        return {
            "disease": disease,
            "biofluid": biofluid,
            "disease_tissue": disease_tissue,
            "normal_tissue": normal_tissue,
            "ev_biomarkers": [],
            "data_sources": [],
            "pipeline_stage": "differential_expression_failed",
            "error": "No differential expression results. Check tissue names match HPA data."
        }
    
    # Step 2: Filter by UniProt transmembrane annotations
    uniprot_dir = os.path.join(data_dir, "uniprot")
    uniprot_df = _load_uniprot_data(uniprot_dir)
    # Build UniProt topology index for downstream EV surface marker analysis
    uniprot_topology_index: Dict[str, Dict[str, Any]] = {}
    if uniprot_df is not None and not uniprot_df.empty:
        uniprot_topology_index = _build_uniprot_topology_index(uniprot_df)

    markers_df = _filter_transmembrane_proteins(markers_df, uniprot_df)
    
    # Step 3: Filter by EVpedia
    evpedia_dir = os.path.join(data_dir, "evpedia")
    evpedia_df = _load_evpedia_data(evpedia_dir, data_dir)
    markers_df = _filter_evpedia_proteins(markers_df, evpedia_df)
    
    # Filter to keep only transmembrane + EVpedia proteins
    final_markers = markers_df[
        (markers_df["is_transmembrane"] == True) &
        (markers_df["in_evpedia"] == True)
    ].copy()
    
    if final_markers.empty:
        # If no matches, return all markers with flags (for debugging)
        final_markers = markers_df.copy()
    
    # Convert to BiomarkerEvidence format
    biomarkers: List[BiomarkerEvidence] = []
    
    for _, row in final_markers.iterrows():
        gene_symbol = str(row["gene_symbol"])
        gene_key = gene_symbol.strip().upper()
        log2fc = float(row["log2fc"]) if pd.notna(row.get("log2fc")) else None
        p_value = float(row["p_value"]) if pd.notna(row.get("p_value")) else None
        fdr = float(row["fdr"]) if pd.notna(row.get("fdr")) else None
        mean_disease = float(row["mean_expression_disease"]) if pd.notna(row.get("mean_expression_disease")) else None
        mean_normal = float(row["mean_expression_normal"]) if pd.notna(row.get("mean_expression_normal")) else None
        
        # Extract z-scores directly from row
        # Initialize with None first to avoid "referenced before assignment" errors
        z_score_pathology = None
        z_score_normal = None
        
        if "z_score_pathology" in row.index:
            z_val = row["z_score_pathology"]
            if z_val is not None and not pd.isna(z_val):
                z_score_pathology = float(z_val)
        
        if "z_score_normal" in row.index:
            z_val = row["z_score_normal"]
            if z_val is not None and not pd.isna(z_val):
                z_score_normal = float(z_val)
        
        # Get UniProt ID if available (from EVpedia or other sources)
        # Initialize with None first to avoid "referenced before assignment" errors
        uniprot = None
        if "uniprot" in row.index and pd.notna(row.get("uniprot")):
            uniprot = str(row["uniprot"])
        elif evpedia_df is not None and not evpedia_df.empty:
            evpedia_match = evpedia_df[evpedia_df["gene_symbol"].str.upper() == gene_symbol.upper()]
            if not evpedia_match.empty and "uniprot" in evpedia_match.columns:
                uniprot = str(evpedia_match.iloc[0]["uniprot"]) if pd.notna(evpedia_match.iloc[0]["uniprot"]) else None
        
        is_transmembrane = bool(row.get("is_transmembrane", False))
        in_evpedia = bool(row.get("in_evpedia", False))

        # UniProt topology-based surface marker assessment
        topo_info = None
        in_uniprot = False
        tm_helix_count = None
        tm_helix_regions = None
        topo_domains = None
        longest_extracellular_domain = None
        has_long_extracellular_domain = None
        is_good_ev_surface_candidate = None

        # Prefer gene-based lookup; if missing, try UniProt entry ID
        if gene_key in uniprot_topology_index:
            topo_info = uniprot_topology_index[gene_key]
        elif uniprot:
            topo_info = uniprot_topology_index.get(str(uniprot).strip().upper())

        protein_symbol = None
        if topo_info is not None:
            in_uniprot = True
            tm_helix_count = topo_info.get("tm_helix_count")
            tm_helix_regions = topo_info.get("tm_helix_regions")
            topo_domains = topo_info.get("topo_domains")
            longest_extracellular_domain = topo_info.get("longest_extracellular_domain")
            has_long_extracellular_domain = topo_info.get("has_long_extracellular_domain")
            is_good_ev_surface_candidate = topo_info.get("is_good_ev_surface_candidate")
            # Derive a compact protein symbol from protein name when available,
            # e.g. extract CD antigens like "CD133" from "CD antigen CD133".
            protein_name = topo_info.get("protein_name") or ""
            if protein_name:
                # Look for CD-style antigen names, case-insensitive
                m_cd = re.search(r"\bCD\d+[A-Z]?\b", protein_name, flags=re.IGNORECASE)
                if m_cd:
                    protein_symbol = m_cd.group(0).upper()
                elif gene_symbol:
                    # Fallback: use gene symbol as protein symbol when no CD alias
                    protein_symbol = gene_symbol
        
        # Evidence level based on statistical significance
        if fdr is not None and fdr < 0.001:
            evidence_level = "high"
        elif fdr is not None and fdr < 0.05:
            evidence_level = "moderate"
        elif p_value is not None and p_value < 0.05:
            evidence_level = "moderate"
        else:
            evidence_level = "low"
        
        # Compute score
        # New scoring scheme (emphasizes differential z-scores and surface/EV evidence):
        #   score = ((z_pathology - z_normal) * 2)
        #           + 5 * TM
        #           + 5 * EVpedia
        score = 0.0
        z_path = z_score_pathology if z_score_pathology is not None else 0.0
        z_norm = z_score_normal if z_score_normal is not None else 0.0
        score += (z_path - z_norm) * 2.0
        if is_transmembrane:
            score += 5.0
        if in_evpedia:
            score += 5.0
        
        # Surface likelihood:
        # - Highest if UniProt topology suggests a long extracellular domain
        # - Otherwise, high if transmembrane
        # - Otherwise, moderate if in EVpedia
        # - Fallback low baseline
        if is_good_ev_surface_candidate:
            surface_likelihood = 0.95
        elif is_transmembrane:
            surface_likelihood = 0.8
        elif in_evpedia:
            surface_likelihood = 0.5
        else:
            surface_likelihood = 0.2
        
        biomarkers.append(
            BiomarkerEvidence(
                gene_symbol=gene_symbol,
                uniprot=uniprot,
                protein_symbol=protein_symbol,
                analyte_type="protein",
                log2fc=log2fc,
                p_value=p_value,
                fdr=fdr,
                mean_expression_disease=mean_disease,
                mean_expression_normal=mean_normal,
                is_transmembrane=is_transmembrane,
                in_evpedia=in_evpedia,
                evidence_level=evidence_level,
                score=float(score),
                supporting_studies=[{
                    "source": "HPA",
                    "disease_tissue": disease_tissue,
                    "normal_tissue": normal_tissue,
                    "log2fc": log2fc,
                    "z_score_pathology": z_score_pathology,
                    "z_score_normal": z_score_normal,
            "p_value": p_value,
                    "fdr": fdr,
                }],
                num_studies=1,
                mean_logfc=log2fc,  # Compatibility alias
                min_p_value=p_value,  # Compatibility alias
                min_fdr=fdr,  # Compatibility alias
                surface_likelihood=surface_likelihood,
                z_score_pathology=z_score_pathology,  # Z-score for pathology/disease tissue
                z_score_normal=z_score_normal,  # Z-score for normal tissue
                in_uniprot=in_uniprot,
                tm_helix_count=tm_helix_count,
                tm_helix_regions=tm_helix_regions,
                topo_domains=topo_domains,
                longest_extracellular_domain=longest_extracellular_domain,
                has_long_extracellular_domain=has_long_extracellular_domain,
                is_good_ev_surface_candidate=is_good_ev_surface_candidate,
            )
        )
    
    # Sort by score (highest first)
    biomarkers_sorted = sorted(biomarkers, key=lambda b: b.score, reverse=True)
    
    # Optional: mark which biomarkers are present in publications_data.xlsx
    pub_index = _load_publication_biomarker_index(data_dir)
    if pub_index is not None:
        by_db = pub_index.get("by_disease_biofluid", {})
        disease_key = (str(disease).strip().lower(), str(biofluid).strip().lower())
        alt_key = (str(disease).strip().lower(), "")
        genes_for_pair = by_db.get(disease_key, set()) | by_db.get(alt_key, set())
        if genes_for_pair:
            for b in biomarkers_sorted:
                g = (b.gene_symbol or "").strip().upper()
                p = (getattr(b, "protein_symbol", None) or "").strip().upper()
                # Also consider compact versions without spaces/hyphens
                g_compact = g.replace(" ", "").replace("-", "") if g else ""
                p_compact = p.replace(" ", "").replace("-", "") if p else ""
                in_pubs = (
                    (g and (g in genes_for_pair or g_compact in genes_for_pair))
                    or (p and (p in genes_for_pair or p_compact in genes_for_pair))
                )
                if in_pubs:
                    b.in_publications_excel = True
                    # Bonus for literature evidence
                    try:
                        b.score = float(b.score) + 1.0
                    except Exception:
                        pass

            # Re-sort after updating scores so rank order matches final score
            biomarkers_sorted = sorted(biomarkers_sorted, key=lambda b: b.score, reverse=True)
    
    # Convert to dicts
    biomarkers_dicts = [asdict(b) for b in biomarkers_sorted]

    
    # Clean NaN values for JSON
    def clean_value_for_json(val):
        if pd.isna(val) if hasattr(pd, 'isna') else (isinstance(val, float) and str(val) == 'nan'):
            return None
        return val
    
    biomarkers_dicts_cleaned = []
    for b_dict in biomarkers_dicts:
        cleaned = {}
        for k, v in b_dict.items():
            if isinstance(v, list):
                cleaned[k] = [
                    {k2: clean_value_for_json(v2) for k2, v2 in item.items()} if isinstance(item, dict) else clean_value_for_json(item)
                    for item in v
                ]
            else:
                cleaned[k] = clean_value_for_json(v)
        biomarkers_dicts_cleaned.append(cleaned)
    
    result = {
        "disease": disease,
        "biofluid": biofluid,
        "disease_tissue": disease_tissue,
        "normal_tissue": normal_tissue,
        "ev_biomarkers": biomarkers_dicts_cleaned,
        "data_sources": ["hpa", "uniprot", "evpedia"],
        "pipeline_mode": "hpa",
        "pipeline_stage": "complete",
    }
    
    return result


def biomarkers_to_markdown(ev_biomarkers: List[Dict[str, Any]], top_n: int = 20) -> str:
    """Generate markdown table of top biomarkers with score components exposed."""
    lines = []
    # Columns:
    # - Core identifiers and stats
    # - Individual score components (effect size, p, FDR, TM bonus, EV bonus)
    # - Surface / topology flags
    lines.append(
        "|Rank|Gene|UniProt|Log2FC|"
        "Effect (|log2FC|×1)|"
        "P-value|P term (-log10·0.5)|"
        "FDR|FDR term (-log10·0.3)|"
        "TM|TM bonus|EVpedia|EV bonus|"
        "Topo surface?|Surface likelihood|"
        "Z-score (Path)|Z-score (Norm)|Evidence|Total score|"
    )
    lines.append(
        "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"
    )
    
    for idx, bm in enumerate(ev_biomarkers[:top_n], 1):
        rank = idx
        gene = bm.get("gene_symbol", "—")
        uniprot = bm.get("uniprot") or "—"
        log2fc = f"{bm.get('log2fc', 0):.2f}" if bm.get("log2fc") is not None else "—"
        raw_log2fc = bm.get("log2fc")

        # P/FDR formatting
        raw_p = bm.get("p_value")
        raw_fdr = bm.get("fdr")
        p_value = f"{raw_p:.2e}" if raw_p is not None else "—"
        fdr = f"{raw_fdr:.2e}" if raw_fdr is not None else "—"

        # Z-scores
        z_score_path = f"{bm.get('z_score_pathology', 0):.2f}" if bm.get("z_score_pathology") is not None else "—"
        z_score_norm = f"{bm.get('z_score_normal', 0):.2f}" if bm.get("z_score_normal") is not None else "—"
        evidence = bm.get("evidence_level", "—")

        # Boolean flags
        is_tm = bool(bm.get("is_transmembrane"))
        in_evpedia = bool(bm.get("in_evpedia"))
        topo_surface = "✓" if bm.get("is_good_ev_surface_candidate") else "—"

        # Score components (recomputed to make table self-documenting)
        effect_term = abs(raw_log2fc) * 1.0 if raw_log2fc is not None else 0.0
        if raw_p is not None and raw_p > 0:
            p_term = -math.log10(raw_p) * 0.5
        else:
            p_term = 0.0
        if raw_fdr is not None and raw_fdr > 0:
            fdr_term = -math.log10(raw_fdr) * 0.3
        else:
            fdr_term = 0.0
        tm_bonus = 5.0 if is_tm else 0.0
        ev_bonus = 5.0 if in_evpedia else 0.0

        tm = "✓" if is_tm else "—"
        evpedia = "✓" if in_evpedia else "—"
        surface_likelihood = bm.get("surface_likelihood")
        surface_likelihood_str = f"{surface_likelihood:.2f}" if surface_likelihood is not None else "—"
        score = f"{bm.get('score', 0):.2f}"
        effect_term_str = f"{effect_term:.2f}"
        p_term_str = f"{p_term:.2f}"
        fdr_term_str = f"{fdr_term:.2f}"
        tm_bonus_str = f"{tm_bonus:.1f}" if tm_bonus else "0.0"
        ev_bonus_str = f"{ev_bonus:.1f}" if ev_bonus else "0.0"
        
        lines.append(
            f"|{rank}|{gene}|{uniprot}|{log2fc}|"
            f"{effect_term_str}|"
            f"{p_value}|{p_term_str}|"
            f"{fdr}|{fdr_term_str}|"
            f"{tm}|{tm_bonus_str}|{evpedia}|{ev_bonus_str}|"
            f"{topo_surface}|{surface_likelihood_str}|"
            f"{z_score_path}|{z_score_norm}|{evidence}|{score}|"
        )
    
    return "\n".join(lines)


def run_cli() -> None:
    """CLI entrypoint for Step 1."""
    parser = argparse.ArgumentParser(description="Step 1 — EV Biomarker Finder (HPA → UniProt → EVpedia)")
    parser.add_argument("--disease", required=True, help="Disease name, e.g. 'CRC'")
    parser.add_argument("--biofluid", required=True, help="Biofluid/tissue, e.g. 'CRC tissue'")
    parser.add_argument("--disease_tissue", default=None, help="Explicit disease tissue name (defaults to '{disease} tissue')")
    parser.add_argument("--normal_tissue", default="Normal tissue", help="Normal tissue name (default: 'Normal tissue')")
    parser.add_argument("--data_dir", default=DATA_DIR, help="Directory with data subdirectories")
    parser.add_argument("--out_json", required=True, help="Path to write JSON result")
    parser.add_argument("--out_md", required=True, help="Path to write markdown summary table")
    
    args = parser.parse_args()
    
    res = aggregate_biomarkers(
        args.disease,
        args.biofluid,
        data_dir=args.data_dir,
        disease_tissue=args.disease_tissue,
        normal_tissue=args.normal_tissue,
    )
    
    os.makedirs(os.path.dirname(os.path.abspath(args.out_json)), exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(args.out_md)), exist_ok=True)
    
    # Clean NaN values for JSON serialization
    def clean_for_json(obj):
        if isinstance(obj, dict):
            return {k: clean_for_json(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [clean_for_json(item) for item in obj]
        elif pd.isna(obj) if hasattr(pd, 'isna') else (isinstance(obj, float) and str(obj) == 'nan'):
            return None
        elif isinstance(obj, (float, int)) and (pd.isna(obj) if hasattr(pd, 'isna') else False):
            return None
        return obj
    
    with open(args.out_json, "w") as f:
        cleaned_res = clean_for_json(res)
        json.dump(cleaned_res, f, indent=2)
    
    md = biomarkers_to_markdown(res["ev_biomarkers"])
    with open(args.out_md, "w") as f:
        f.write(md + "\n")
    
    print(f"✓ Found {len(res['ev_biomarkers'])} biomarkers")
    print(f"✓ Wrote JSON to: {args.out_json}")
    print(f"✓ Wrote markdown to: {args.out_md}")


if __name__ == "__main__":
    run_cli()

