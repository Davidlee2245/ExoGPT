"""
HPA Pathology Data Converter
-----------------------------

Converts HPA pathology API data into format expected by Step 1 pipeline.
"""

from __future__ import annotations
import os
import pandas as pd
from typing import Optional, Dict, Union, Tuple


def convert_hpa_pathology_to_expression_format(
    hpa_csv_path: str,
    cancer_name: str,
    output_dir: str = "data/hpa",
    return_diagnostics: bool = False
) -> Union[Optional[str], Tuple[Optional[str], Dict]]:
    """
    Convert HPA pathology CSV to expression format expected by pipeline.
    
    HPA pathology data has columns like:
    - Gene
    - Ensembl
    - RNA tissue category
    - Pathology score (<cancer_name>)
    - Protein class
    etc.
    
    Converts to format with:
    - gene_symbol
    - tissue (disease vs normal)
    - expression_value (pathology score)
    
    Args:
        hpa_csv_path: Path to HPA pathology CSV file
        cancer_name: Cancer name used in search (e.g., "colorectal cancer")
        output_dir: Output directory for converted file
        return_diagnostics: If True, return tuple (path, diagnostics_dict) instead of just path
    
    Returns:
        Path to converted CSV file, or None if conversion failed
        If return_diagnostics=True, returns (path, diagnostics) or (None, diagnostics)
    """
    diagnostics = {
        "csv_path": hpa_csv_path,
        "cancer_name": cancer_name,
        "errors": [],
        "warnings": [],
        "info": []
    }
    
    if not os.path.exists(hpa_csv_path):
        diagnostics["errors"].append(f"CSV file does not exist: {hpa_csv_path}")
        if return_diagnostics:
            return None, diagnostics
        return None
    
    try:
        # Load HPA pathology data
        df = pd.read_csv(hpa_csv_path)
        
        if df.empty:
            error_msg = f"HPA CSV file is empty: {hpa_csv_path}"
            diagnostics["errors"].append(error_msg)
            if return_diagnostics:
                return None, diagnostics
            return None
        
        diagnostics["info"].append(f"Loaded HPA data: {len(df)} rows, {len(df.columns)} columns")
        
        # Normalize column names
        original_columns = list(df.columns)
        df.columns = [c.lower().strip() for c in df.columns]
        normalized_columns = list(df.columns)
        
        diagnostics["info"].append(f"Available columns: {', '.join(normalized_columns)}")
        
        # Show first few rows for debugging
        if len(df) > 0:
            sample_data = {}
            for col in normalized_columns[:5]:  # Show first 5 columns
                if col in df.columns:
                    sample_data[col] = str(df.iloc[0][col])
            diagnostics["info"].append(f"Sample first row: {sample_data}")
        
        # Find gene column (could be "gene", "gene_symbol", "gene name", etc.)
        gene_col = None
        for col in ["gene", "gene_symbol", "gene name", "gene_name", "gene symbol"]:
            if col in df.columns:
                gene_col = col
                break
        
        if not gene_col:
            error_msg = f"Could not find gene column. Available columns: {normalized_columns}"
            diagnostics["errors"].append(error_msg)
            if return_diagnostics:
                return None, diagnostics
            return None
        
        diagnostics["info"].append(f"Found gene column: '{gene_col}'")
        
        # Find pathology score column (e.g., "pathology score (colorectal cancer)")
        pathology_col = None
        cancer_lower = cancer_name.lower()
        for col in df.columns:
            if "pathology" in col.lower() and "score" in col.lower():
                # Check if it matches the cancer name
                if cancer_lower in col.lower() or any(word in col.lower() for word in cancer_lower.split()):
                    pathology_col = col
                    break
        
        # If no exact match, try to find any pathology score column
        if not pathology_col:
            for col in df.columns:
                if "pathology" in col.lower() and "score" in col.lower():
                    pathology_col = col
                    break
        
        if not pathology_col:
            # Try to find any numeric column that might be a score
            numeric_cols = df.select_dtypes(include=['number']).columns.tolist()
            if len(numeric_cols) > 0:
                pathology_col = numeric_cols[0]
                diagnostics["warnings"].append(f"Using first numeric column as pathology score: '{pathology_col}'")
            else:
                error_msg = f"Could not find pathology score column. Available columns: {normalized_columns}. Numeric columns: {numeric_cols}"
                diagnostics["errors"].append(error_msg)
                if return_diagnostics:
                    return None, diagnostics
                return None
        
        diagnostics["info"].append(f"Found pathology score column: '{pathology_col}'")
        
        # Extract gene symbols and pathology scores
        converted_data = []
        valid_scores = []
        
        # First pass: collect all valid scores to compute baseline
        for _, row in df.iterrows():
            gene = row[gene_col]
            if pd.isna(gene):
                continue
            
            score = row[pathology_col]
            if pd.isna(score):
                continue
            
            try:
                score_val = float(score)
                if score_val > 0:  # Only positive scores
                    valid_scores.append(score_val)
            except (ValueError, TypeError):
                continue
        
        # Compute baseline (median of all scores, or 1.0 if no valid scores)
        baseline = pd.Series(valid_scores).median() if valid_scores else 1.0
        if baseline <= 0:
            baseline = 1.0
        
        # Second pass: create expression data
        for _, row in df.iterrows():
            gene = row[gene_col]
            if pd.isna(gene):
                continue
            
            score = row[pathology_col]
            if pd.isna(score):
                continue
            
            # Convert to numeric
            try:
                score_val = float(score)
                if score_val <= 0:
                    continue
            except (ValueError, TypeError):
                continue
            
            # Create rows for disease and normal
            # Disease: use pathology score (higher = more associated with cancer)
            # Normal: use baseline (represents normal tissue expression)
            converted_data.append({
                "gene_symbol": str(gene).strip(),
                "tissue": f"{cancer_name} tissue",
                "expression_value": score_val
            })
            
            converted_data.append({
                "gene_symbol": str(gene).strip(),
                "tissue": "Normal tissue",
                "expression_value": baseline
            })
        
        if not converted_data:
            error_msg = f"No valid data extracted. Found {len(df)} rows but no valid gene-score pairs. Gene column: '{gene_col}', Pathology column: '{pathology_col}', Valid scores: {len(valid_scores)}"
            diagnostics["errors"].append(error_msg)
            if return_diagnostics:
                return None, diagnostics
            return None
        
        diagnostics["info"].append(f"Converted {len(converted_data) // 2} genes to expression format")
        
        # Create output DataFrame
        converted_df = pd.DataFrame(converted_data)
        
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        # Generate output filename
        from exo_gpt.hpa_downloader import slugify
        slug = slugify(cancer_name)
        output_path = os.path.join(output_dir, f"{slug}_expression.csv")
        
        # Save converted file
        converted_df.to_csv(output_path, index=False)
        
        diagnostics["info"].append(f"Saved converted file to: {output_path}")
        
        if return_diagnostics:
            return output_path, diagnostics
        return output_path
        
    except Exception as e:
        import traceback
        error_msg = f"Exception during conversion: {str(e)}"
        diagnostics["errors"].append(error_msg)
        diagnostics["errors"].append(f"Traceback: {traceback.format_exc()}")
        if return_diagnostics:
            return None, diagnostics
        return None


def map_disease_to_hpa_search_term(disease: str) -> Optional[str]:
    """
    Map disease name to HPA API search term.
    
    Args:
        disease: Disease name (e.g., "CRC", "Metastatic melanoma")
    
    Returns:
        HPA search term (e.g., "colorectal cancer", "melanoma") or None
    """
    disease_lower = disease.lower().strip()
    
    # Mapping dictionary - using HPA API exact terminology
    # HPA uses specific cancer type names - these must match exactly
    # Based on HPA pathology database terminology
    mappings = {
        "crc": "Colorectal cancer",
        "colorectal": "Colorectal cancer",
        "colorectal cancer": "Colorectal cancer",
        "melanoma": "Skin Cutaneous Melanoma",  # HPA exact term
        "metastatic melanoma": "Skin Cutaneous Melanoma",
        "malignant melanoma": "Skin Cutaneous Melanoma",
        "skin cancer": "Skin Cutaneous Melanoma",
        "breast": "Breast cancer",
        "breast cancer": "Breast cancer",
        "lung": "Lung cancer",
        "lung cancer": "Lung cancer",
        "prostate": "Prostate cancer",
        "prostate cancer": "Prostate cancer",
        "pancreatic": "Pancreatic cancer",
        "pancreatic cancer": "Pancreatic cancer",
        "gastric": "Gastric cancer",
        "gastric cancer": "Gastric cancer",
        "liver": "Liver cancer",
        "liver cancer": "Liver cancer",
        "ovarian": "Ovarian cancer",
        "ovarian cancer": "Ovarian cancer",
        "cervical": "Cervical cancer",
        "cervical cancer": "Cervical cancer",
        "endometrial": "Endometrial cancer",
        "endometrial cancer": "Endometrial cancer",
        "renal": "Renal cancer",
        "renal cancer": "Renal cancer",
        "kidney": "Renal cancer",
        "kidney cancer": "Renal cancer",
    }
    
    # Try exact match first
    if disease_lower in mappings:
        return mappings[disease_lower]
    
    # Try partial match
    for key, value in mappings.items():
        if key in disease_lower or disease_lower in key:
            return value
    
    # If no mapping found, try to extract cancer type from disease name
    # Common patterns: "X cancer", "X carcinoma", etc.
    if "cancer" in disease_lower:
        # Extract the part before "cancer"
        parts = disease_lower.split("cancer")
        if parts[0].strip():
            return f"{parts[0].strip()} cancer"
    
    if "carcinoma" in disease_lower:
        parts = disease_lower.split("carcinoma")
        if parts[0].strip():
            return f"{parts[0].strip()} cancer"
    
    # Default: return None (will use legacy pipeline)
    return None

