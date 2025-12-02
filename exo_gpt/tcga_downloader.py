"""
TCGA Data Downloader for mRNA Biomarker Discovery
Downloads TCGA mRNA expression data for specific cancer types
"""

from __future__ import annotations
import os
import json
import time
from typing import Optional, Dict, Any
import pandas as pd
import requests

# TCGA project code mapping (cancer type -> TCGA project code)
TCGA_PROJECT_CODES = {
    "colorectal cancer": "COAD",
    "breast cancer": "BRCA",
    "lung cancer": "LUAD",
    "prostate cancer": "PRAD",
    "skin cutaneous melanoma": "SKCM",
    "pancreatic cancer": "PAAD",
    "gastric cancer": "STAD",
    "liver cancer": "LIHC",
    "ovarian cancer": "OV",
    "cervical cancer": "CESC",
    "endometrial cancer": "UCEC",
    "renal cancer": "KIRC",
    "bladder cancer": "BLCA",
    "thyroid cancer": "THCA",
    "brain cancer": "GBM",
    "metastatic melanoma": "SKCM",
    "melanoma": "SKCM",
}

# GDC API endpoints
GDC_API_BASE = "https://api.gdc.cancer.gov"
GDC_DATA_BASE = "https://api.gdc.cancer.gov/data"


def map_disease_to_tcga_project(disease: str) -> Optional[str]:
    """
    Map disease name to TCGA project code.
    
    Args:
        disease: Disease name (e.g., "colorectal cancer", "melanoma")
    
    Returns:
        TCGA project code (e.g., "COAD", "SKCM") or None if not found
    """
    disease_lower = disease.lower().strip()
    
    # Direct mapping
    if disease_lower in TCGA_PROJECT_CODES:
        return TCGA_PROJECT_CODES[disease_lower]
    
    # Fuzzy matching
    for key, code in TCGA_PROJECT_CODES.items():
        if key in disease_lower or disease_lower in key:
            return code
    
    return None


def download_tcga_mrna_data(
    disease: str,
    output_dir: str,
    progress_callback: Optional[Any] = None
) -> Optional[str]:
    """
    Download TCGA mRNA expression data for a given disease using GDC API.
    
    This function uses the GDC (Genomic Data Commons) API to download
    TCGA mRNA expression data directly from the GDC data portal.
    
    Args:
        disease: Disease name (e.g., "colorectal cancer")
        output_dir: Output directory (e.g., "./data/mRNA/databases/TCGA")
        progress_callback: Optional callback for progress updates
    
    Returns:
        Path to downloaded CSV file, or None if download failed
    """
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Map disease to TCGA project code
    project_code = map_disease_to_tcga_project(disease)
    if not project_code:
        if progress_callback:
            progress_callback(f"⚠️  No TCGA project code found for: {disease}", 0.0)
        print(f"⚠️  No TCGA project code found for: {disease}")
        return None
    
    if progress_callback:
        progress_callback(f"Downloading TCGA data for {disease} (project: {project_code})...", 10.0)
    print(f"📥 Downloading TCGA data for {disease} (project: {project_code})...")
    
    # Output file path
    disease_safe = disease.lower().replace(" ", "_").replace("/", "_")
    output_file = os.path.join(output_dir, f"tcga_{project_code}_{disease_safe}.csv")
    
    # Check if file already exists and has data
    if os.path.exists(output_file):
        try:
            df_check = pd.read_csv(output_file, low_memory=False)
            if len(df_check) > 0:
                if progress_callback:
                    progress_callback(f"✓ TCGA data already exists: {output_file}", 100.0)
                print(f"✓ TCGA data already exists: {output_file}")
                return output_file
        except Exception:
            pass  # File exists but might be corrupted, re-download
    
    try:
        # Step 1: Query GDC API for files
        if progress_callback:
            progress_callback(f"Querying GDC API for TCGA-{project_code}...", 20.0)
        
        # Query for Gene Expression Quantification files (mRNA expression data)
        # Note: GDC has moved from HTSeq to STAR workflow, so we don't filter by workflow_type
        query = {
            "filters": {
                "op": "and",
                "content": [
                    {
                        "op": "=",
                        "content": {
                            "field": "cases.project.project_id",
                            "value": [f"TCGA-{project_code}"]
                        }
                    },
                    {
                        "op": "=",
                        "content": {
                            "field": "files.data_category",
                            "value": ["Transcriptome Profiling"]
                        }
                    },
                    {
                        "op": "=",
                        "content": {
                            "field": "files.data_type",
                            "value": ["Gene Expression Quantification"]
                        }
                    }
                ]
            },
            "format": "JSON",
            "size": 10000,
            "fields": "file_id,file_name,data_format,cases.samples.submitter_id,cases.samples.sample_type"
        }
        
        files_endpoint = f"{GDC_API_BASE}/files"
        response = requests.post(
            files_endpoint,
            json=query,
            headers={"Content-Type": "application/json"},
            timeout=60
        )
        
        if response.status_code != 200:
            raise Exception(f"GDC API query failed: {response.status_code}")
        
        files_data = response.json()
        files_list = files_data.get("data", {}).get("hits", [])
        
        if not files_list:
            raise Exception(f"No TCGA files found for project TCGA-{project_code}")
        
        if progress_callback:
            progress_callback(f"Found {len(files_list)} TCGA files, downloading...", 30.0)
        print(f"  Found {len(files_list)} TCGA files")
        
        # Build sample type mapping from file metadata
        sample_type_map = {}
        file_metadata = {}
        for file_info in files_list:
            file_id = file_info.get("id")
            file_name = file_info.get("file_name", "")
            cases = file_info.get("cases", [])
            
            # Extract sample information from cases
            for case in cases:
                samples = case.get("samples", [])
                for sample in samples:
                    sample_id = sample.get("submitter_id", "")
                    sample_type = sample.get("sample_type", "")
                    # Map TCGA sample types to our format
                    if "tumor" in sample_type.lower() or "primary" in sample_type.lower():
                        sample_type_map[sample_id] = "Primary Tumor"
                    elif "normal" in sample_type.lower() or "solid tissue" in sample_type.lower():
                        sample_type_map[sample_id] = "Solid Tissue Normal"
            
            # Store file metadata
            file_metadata[file_id] = {
                "file_name": file_name,
                "cases": cases
            }
        
        # Step 3: Download and process files (limit to first 50 for performance)
        if progress_callback:
            progress_callback("Downloading expression data (this may take a while)...", 50.0)
        
        all_data = []
        max_files = min(50, len(files_list))  # Limit to 50 files for performance
        
        for i, file_info in enumerate(files_list[:max_files]):
            try:
                file_id = file_info.get("id")
                file_name = file_info.get("file_name", "")
                
                # Get sample ID from file metadata
                sample_id = None
                cases = file_info.get("cases", [])
                if cases:
                    samples = cases[0].get("samples", [])
                    if samples:
                        sample_id = samples[0].get("submitter_id", "")
                
                # Determine sample type
                sample_type = sample_type_map.get(sample_id, "Primary Tumor")  # Default to tumor
                
                # Download file data
                data_url = f"{GDC_DATA_BASE}/{file_id}"
                file_response = requests.get(data_url, timeout=120)
                
                if file_response.status_code == 200:
                    # Parse TSV data
                    lines = file_response.text.strip().split('\n')
                    if not lines:
                        continue
                    
                    # Skip comment lines (starting with #)
                    data_lines = [line for line in lines if not line.strip().startswith('#')]
                    if not data_lines:
                        continue
                    
                    # Find header line (first non-comment line)
                    header_line = data_lines[0]
                    header = header_line.split('\t')
                    
                    # Find gene_id and gene_name columns
                    gene_id_idx = None
                    gene_name_idx = None
                    fpkm_col_idx = None
                    
                    for idx, col in enumerate(header):
                        col_lower = col.lower()
                        if 'gene_id' in col_lower and 'gene_name' not in col_lower:
                            gene_id_idx = idx
                        if 'gene_name' in col_lower:
                            gene_name_idx = idx
                        # Prefer FPKM over raw counts
                        if 'fpkm' in col_lower and 'unstranded' in col_lower:
                            fpkm_col_idx = idx
                        elif 'fpkm' in col_lower and fpkm_col_idx is None:
                            fpkm_col_idx = idx
                    
                    # Fallback to unstranded counts if no FPKM
                    if fpkm_col_idx is None:
                        for idx, col in enumerate(header):
                            col_lower = col.lower()
                            if 'unstranded' in col_lower and 'fpkm' not in col_lower:
                                fpkm_col_idx = idx
                                break
                    
                    if gene_id_idx is None or gene_name_idx is None:
                        print(f"  ⚠️  Could not find gene_id or gene_name columns in {file_name}")
                        continue
                    
                    # Process data lines (skip header)
                    for line in data_lines[1:]:
                        parts = line.split('\t')
                        if len(parts) < max(gene_id_idx, gene_name_idx) + 1:
                            continue
                        
                        gene_id = parts[gene_id_idx] if gene_id_idx < len(parts) else ""
                        gene_name = parts[gene_name_idx] if gene_name_idx < len(parts) else ""
                        
                        # Skip if no gene name or if it's a special row (N_unmapped, etc.)
                        if not gene_name or gene_name == "" or gene_id.startswith('N_'):
                            continue
                        
                        # Extract Ensembl ID (remove version number)
                        ensembl_id = gene_id.split('.')[0] if '.' in gene_id else gene_id
                        
                        # Get FPKM or count value
                        expr_value = None
                        if fpkm_col_idx is not None and fpkm_col_idx < len(parts):
                            try:
                                expr_value = float(parts[fpkm_col_idx])
                            except (ValueError, IndexError):
                                pass
                        
                        # If FPKM not found, try to find any numeric column after gene_name
                        if expr_value is None:
                            for val in parts[gene_name_idx + 1:]:
                                try:
                                    expr_value = float(val)
                                    break
                                except ValueError:
                                    continue
                        
                        if expr_value is not None and expr_value >= 0:
                            # If it's already FPKM, convert to log2(FPKM + 1)
                            # If it's counts, also convert to log2(count + 1)
                            import numpy as np
                            log2_fpkm = np.log2(expr_value + 1)
                            
                            all_data.append({
                                "gene_symbol": gene_name,
                                "ensembl_id": ensembl_id,
                                "sample_type": sample_type,
                                "log2_fpkm": log2_fpkm,
                                "project_code": project_code
                            })
                
                if progress_callback:
                    progress = 50.0 + (i + 1) / max_files * 40.0
                    progress_callback(f"Processed {i+1}/{max_files} files...", progress)
                
            except Exception as e:
                print(f"  ⚠️  Error processing file {file_id}: {e}")
                continue
        
        # Step 4: Save to CSV
        if progress_callback:
            progress_callback("Saving TCGA data to CSV...", 95.0)
        
        if not all_data:
            raise Exception("No data extracted from TCGA files")
        
        df = pd.DataFrame(all_data)
        df.to_csv(output_file, index=False)
        
        if progress_callback:
            progress_callback(f"✓ TCGA data downloaded: {output_file} ({len(df)} records)", 100.0)
        print(f"✓ TCGA data downloaded: {output_file}")
        print(f"  Total records: {len(df)}")
        print(f"  Unique genes: {df['gene_symbol'].nunique()}")
        print(f"  Tumor samples: {len(df[df['sample_type'].str.contains('Tumor', case=False)])}")
        print(f"  Normal samples: {len(df[df['sample_type'].str.contains('Normal', case=False)])}")
        
        return output_file
        
    except Exception as e:
        error_msg = f"Error downloading TCGA data: {str(e)}"
        if progress_callback:
            progress_callback(f"❌ {error_msg}", 0.0)
        print(f"❌ {error_msg}")
        import traceback
        traceback.print_exc()
        return None


def load_tcga_data(tcga_file: str, gene_symbol: str) -> Optional[float]:
    """
    Load TCGA data from CSV and calculate logFC for a gene.
    
    Args:
        tcga_file: Path to TCGA CSV file
        gene_symbol: Gene symbol to query
    
    Returns:
        log2 fold change (tumor vs normal) or None if not found
    """
    if not os.path.exists(tcga_file):
        return None
    
    try:
        df = pd.read_csv(tcga_file, low_memory=False)
        
        # Normalize column names
        df.columns = [c.lower().strip() for c in df.columns]
        
        # Filter for gene
        if "gene_symbol" in df.columns:
            gene_df = df[df["gene_symbol"].str.upper() == gene_symbol.upper()]
            
            if gene_df.empty:
                return None
            
            # Calculate logFC: tumor vs normal
            # Look for sample_type or similar column
            if "sample_type" in gene_df.columns:
                tumor_df = gene_df[gene_df["sample_type"].str.contains("Tumor", case=False, na=False)]
                normal_df = gene_df[gene_df["sample_type"].str.contains("Normal", case=False, na=False)]
                
                # Get expression values
                expr_col = None
                for col in ["log2_fpkm", "log2_expression", "expression_value", "fpkm"]:
                    if col in gene_df.columns:
                        expr_col = col
                        break
                
                if expr_col:
                    tumor_expr = tumor_df[expr_col].dropna()
                    normal_expr = normal_df[expr_col].dropna()
                    
                    if len(tumor_expr) > 0 and len(normal_expr) > 0:
                        tumor_mean = tumor_expr.mean()
                        normal_mean = normal_expr.mean()
                        
                        # Calculate logFC
                        if normal_mean > 0:
                            logfc = tumor_mean - normal_mean
                            return float(logfc)
                        elif tumor_mean > 0:
                            # If normal is 0, use tumor value as logFC
                            return float(tumor_mean)
            
            # Fallback: if no sample_type, try to calculate from available columns
            if "log2_fpkm" in gene_df.columns:
                # Assume all samples are tumors if no sample_type
                mean_expr = gene_df["log2_fpkm"].dropna().mean()
                if not pd.isna(mean_expr):
                    return float(mean_expr)
    
    except Exception as e:
        print(f"⚠️  Error loading TCGA data: {e}")
    
    return None

