"""
Step 1 — EV Biomarker Finder (Enhanced with Multi-Source Querying)
-------------------------------------------------------------------

Automatically searches across:
- ./data/publications/ - Extracted from published articles
- ./data/databases/ - Pre-packaged from Vesiclepedia, ExoCarta, EVpedia, EVmiRNA, exoRBase
- ./data/experiments/ - User-generated experimental data

Prioritizes local files for speed, with optional database download capability.
"""

from __future__ import annotations

import argparse
import json
import os
from dataclasses import dataclass, asdict
from typing import List, Dict, Any, Optional, Tuple

import numpy as np
import pandas as pd

DATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "data"))


@dataclass
class NormalizedDisease:
    name: str
    doid: Optional[str] = None
    mesh: Optional[str] = None


@dataclass
class BiomarkerEvidence:
    gene_symbol: str
    uniprot: Optional[str]
    analyte_type: str
    num_studies: int
    mean_logfc: Optional[float]
    min_p_value: Optional[float]
    min_fdr: Optional[float]
    evidence_level: str
    surface_likelihood: Optional[float]
    score: float
    supporting_studies: List[Dict[str, Any]]


def _load_disease_mapping(mapping_path: str) -> pd.DataFrame:
    """Load disease name normalization table."""
    if not os.path.exists(mapping_path):
        return pd.DataFrame(columns=["disease_name", "doid", "mesh", "synonyms"])
    return pd.read_csv(mapping_path, sep=None, engine="python")


def normalize_disease_name(disease: str, mapping_df: pd.DataFrame) -> NormalizedDisease:
    """Normalize disease name using mapping table."""
    if mapping_df.empty:
        return NormalizedDisease(name=disease)

    d = disease.lower().strip()
    for _, row in mapping_df.iterrows():
        names = [str(row.get("disease_name", "")).lower()]
        syn = str(row.get("synonyms", "")).lower()
        if syn:
            names.extend([s.strip() for s in syn.split("|") if s.strip()])
        if any(d == n or d in n for n in names):
            return NormalizedDisease(
                name=row.get("disease_name", disease),
                doid=row.get("doid") if pd.notna(row.get("doid")) else None,
                mesh=row.get("mesh") if pd.notna(row.get("mesh")) else None,
            )
    return NormalizedDisease(name=disease)


def _load_ev_tables(data_dir: str, search_subdirs: bool = True) -> List[Tuple[str, pd.DataFrame]]:
    """
    Load EV-like tables from data_dir and all subdirectories.
    
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


def _filter_table_for_disease_biofluid(
    df: pd.DataFrame, disease: str, biofluid: str
) -> pd.DataFrame:
    """Filter DataFrame by disease and biofluid (case-insensitive)."""
    disease_lower = disease.lower().strip()
    biofluid_lower = biofluid.lower().strip()
    
    # Convert to string and handle NaN
    disease_col = df["disease"].astype(str).str.lower()
    biofluid_col = df["biofluid"].astype(str).str.lower()
    
    # Create boolean masks
    disease_mask = disease_col.str.contains(disease_lower, na=False, case=False)
    biofluid_mask = biofluid_col.str.contains(biofluid_lower, na=False, case=False)
    
    # Combine masks
    mask = disease_mask & biofluid_mask
    return df[mask].copy()


def _estimate_surface_likelihood(row: pd.Series) -> Optional[float]:
    """Estimate if protein is likely surface-exposed based on annotations."""
    if "protein" not in str(row.get("analyte_type", "")).lower():
        return None
    
    score = 0.0
    text_fields = []
    for col in ["notes", "method", "gene_symbol", "uniprot"]:
        if col in row.index and pd.notna(row[col]):
            text_fields.append(str(row[col]).lower())
    
    combined_text = " ".join(text_fields)
    
    # Surface protein keywords
    surface_keywords = [
        "surface", "membrane", "extracellular", "ecd", "transmembrane",
        "cd274", "pd-l1", "pdcd1", "pd-1", "ctla4", "cd63", "cd81", "cd9",
        "receptor", "ligand", "adhesion", "integrin", "tetraspanin"
    ]
    for keyword in surface_keywords:
        if keyword in combined_text:
            score += 0.1
    
    # Exosome marker proteins (typically surface)
    exosome_markers = ["cd63", "cd81", "cd9", "tsg101", "alix", "flot"]
    for marker in exosome_markers:
        if marker in combined_text:
            score += 0.15
    
    return min(score, 1.0) if score > 0 else 0.2  # Default 0.2 if no keywords


def _compute_biomarker_score(
    num_studies: int,
    mean_logfc: Optional[float],
    min_p: Optional[float],
    min_fdr: Optional[float],
    surface_likelihood: Optional[float],
) -> float:
    """Compute composite score for ranking biomarkers."""
    score = 0.0
    
    # Number of studies (more = better)
    score += min(num_studies * 2.0, 10.0)
    
    # Effect size (logFC)
    if mean_logfc is not None:
        score += abs(mean_logfc) * 3.0
    
    # Statistical significance
    if min_p is not None and min_p > 0:
        score += -np.log10(min_p) * 0.5
    
    if min_fdr is not None and min_fdr > 0:
        score += -np.log10(min_fdr) * 0.3
    
    # Surface likelihood bonus (for nanobinder targets)
    if surface_likelihood is not None:
        score += surface_likelihood * 5.0
    
    return float(score)


def aggregate_biomarkers(
    disease: str,
    biofluid: str,
    data_dir: str = DATA_DIR,
    search_subdirs: bool = True,
) -> Dict[str, Any]:
    """
    Aggregate and rank EV biomarkers across all data sources.
    
    Automatically searches:
    - publications/ subdirectory
    - databases/ subdirectory  
    - experiments/ subdirectory
    - root data_dir (legacy)
    
    Returns JSON-serializable dict with ranked biomarkers.
    """
    disease_map_path = os.path.join(data_dir, "disease_mapping.csv")
    mapping_df = _load_disease_mapping(disease_map_path)
    norm = normalize_disease_name(disease, mapping_df)

    # Load all tables from all sources
    tables = _load_ev_tables(data_dir, search_subdirs=search_subdirs)
    
    if not tables:
        return {
            "disease": disease,
            "biofluid": biofluid,
            "normalized_disease_ids": {
                "doid": norm.doid,
                "mesh": norm.mesh,
            },
            "ev_biomarkers": [],
            "data_sources": [],
        }

    records: List[Dict[str, Any]] = []
    data_sources = set()
    
    for dataset_name, df in tables:
        sub = _filter_table_for_disease_biofluid(df, norm.name, biofluid)
        if sub.empty:
            continue
        sub = sub.copy()
        sub["dataset"] = dataset_name
        records.append(sub)
        data_sources.add(dataset_name)

    if not records:
        return {
            "disease": norm.name,
            "biofluid": biofluid,
            "normalized_disease_ids": {
                "doid": norm.doid,
                "mesh": norm.mesh,
            },
            "ev_biomarkers": [],
            "data_sources": list(data_sources),
        }

    all_df = pd.concat(records, ignore_index=True)

    # Normalize column existence
    for col in ["uniprot", "logfc", "p_value", "fdr", "method", "pmid", "study_id"]:
        if col not in all_df.columns:
            all_df[col] = pd.NA

    grouped = all_df.groupby(["gene_symbol", "analyte_type"], dropna=False)

    biomarkers: List[BiomarkerEvidence] = []
    for (gene_symbol, analyte_type), g in grouped:
        uniprot_vals = g["uniprot"].dropna().unique()
        uniprot = uniprot_vals[0] if len(uniprot_vals) > 0 else None

        num_studies = g["study_id"].nunique() if "study_id" in g.columns and g["study_id"].notna().any() else g["dataset"].nunique()
        mean_logfc = float(g["logfc"].astype(float).mean()) if "logfc" in g.columns and g["logfc"].notna().any() else None

        min_p = float(g["p_value"].astype(float).min()) if "p_value" in g.columns and g["p_value"].notna().any() else None
        min_fdr = float(g["fdr"].astype(float).min()) if "fdr" in g.columns and g["fdr"].notna().any() else None

        surface_scores = g.apply(_estimate_surface_likelihood, axis=1)
        if surface_scores.notna().any():
            surface_likelihood = float(surface_scores.mean())
        else:
            surface_likelihood = None

        if num_studies >= 5:
            evidence_level = "high"
        elif num_studies >= 2:
            evidence_level = "moderate"
        else:
            evidence_level = "low"

        score = _compute_biomarker_score(
            num_studies=num_studies,
            mean_logfc=mean_logfc,
            min_p=min_p,
            min_fdr=min_fdr,
            surface_likelihood=surface_likelihood,
        )

        supporting_studies = []
        for _, row in g.iterrows():
            # Convert NaN to None for JSON serialization
            def clean_value(val):
                if val is None:
                    return None
                try:
                    if pd.isna(val):
                        return None
                except (TypeError, ValueError):
                    pass
                if isinstance(val, float):
                    import math
                    if math.isnan(val):
                        return None
                if isinstance(val, str) and val.lower() in ['nan', 'none', '']:
                    return None
                return val
            
            supporting_studies.append({
                "study_id": clean_value(row.get("study_id")),
                "dataset": clean_value(row.get("dataset")),
                "logFC": clean_value(row.get("logfc")),
                "p_value": clean_value(row.get("p_value")),
                "fdr": clean_value(row.get("fdr")),
                "method": clean_value(row.get("method")),
                "biofluid": clean_value(row.get("biofluid")),
                "pmid": clean_value(row.get("pmid")),
            })

        # Clean NaN values before creating dataclass
        def clean_float(val):
            if val is None:
                return None
            if pd.isna(val) if hasattr(pd, 'isna') else (isinstance(val, float) and str(val) == 'nan'):
                return None
            return float(val) if isinstance(val, (int, float)) else val
        
        biomarkers.append(
            BiomarkerEvidence(
                gene_symbol=str(gene_symbol),
                uniprot=str(uniprot) if uniprot and not pd.isna(uniprot) else None,
                analyte_type=str(analyte_type),
                num_studies=int(num_studies),
                mean_logfc=clean_float(mean_logfc),
                min_p_value=clean_float(min_p),
                min_fdr=clean_float(min_fdr),
                evidence_level=evidence_level,
                surface_likelihood=clean_float(surface_likelihood),
                score=float(score),
                supporting_studies=supporting_studies,
            )
        )

    # Sort: proteins first, then by score
    biomarkers_sorted = sorted(
        biomarkers,
        key=lambda b: (0 if "protein" not in b.analyte_type.lower() else 1, b.score),
        reverse=True,
    )

    # Convert dataclasses to dicts and clean NaN values
    def clean_value_for_json(val):
        """Convert NaN/NA values to None for JSON serialization"""
        if pd.isna(val) if hasattr(pd, 'isna') else (isinstance(val, float) and str(val) == 'nan'):
            return None
        return val
    
    biomarkers_dicts = []
    for b in biomarkers_sorted:
        b_dict = asdict(b)
        # Clean all values in the biomarker dict
        cleaned_dict = {}
        for k, v in b_dict.items():
            if isinstance(v, list):
                cleaned_dict[k] = [
                    {k2: clean_value_for_json(v2) for k2, v2 in item.items()} if isinstance(item, dict) else clean_value_for_json(item)
                    for item in v
                ]
            else:
                cleaned_dict[k] = clean_value_for_json(v)
        biomarkers_dicts.append(cleaned_dict)
    
    result = {
        "disease": norm.name,
        "biofluid": biofluid,
        "normalized_disease_ids": {
            "doid": clean_value_for_json(norm.doid),
            "mesh": clean_value_for_json(norm.mesh),
        },
        "ev_biomarkers": biomarkers_dicts,
        "data_sources": sorted(list(data_sources)),
    }
    return result


def biomarkers_to_markdown(ev_biomarkers: List[Dict[str, Any]], top_n: int = 20) -> str:
    """Generate markdown table of top biomarkers."""
    lines = []
    lines.append("|Rank|Gene|UniProt|Analyte|NumStudies|MeanLogFC|MinP|MinFDR|Evidence|SurfaceLikelihood|Score|")
    lines.append("|---|---|---|---|---|---|---|---|---|---|---|")
    
    for idx, bm in enumerate(ev_biomarkers[:top_n], 1):
        rank = idx
        gene = bm.get("gene_symbol", "—")
        uniprot = bm.get("uniprot") or "—"
        analyte = bm.get("analyte_type", "—")
        num_studies = bm.get("num_studies", 0)
        mean_logfc = f"{bm.get('mean_logfc', 0):.2f}" if bm.get("mean_logfc") is not None else "—"
        min_p = f"{bm.get('min_p_value', 0):.2e}" if bm.get("min_p_value") is not None else "—"
        min_fdr = f"{bm.get('min_fdr', 0):.2e}" if bm.get("min_fdr") is not None else "—"
        evidence = bm.get("evidence_level", "—")
        surface = f"{bm.get('surface_likelihood', 0):.2f}" if bm.get("surface_likelihood") is not None else "—"
        score = f"{bm.get('score', 0):.2f}"
        
        lines.append(f"|{rank}|{gene}|{uniprot}|{analyte}|{num_studies}|{mean_logfc}|{min_p}|{min_fdr}|{evidence}|{surface}|{score}|")
    
    return "\n".join(lines)


def run_cli() -> None:
    """CLI entrypoint for Step 1."""
    parser = argparse.ArgumentParser(description="Step 1 — EV Biomarker Finder")
    parser.add_argument("--disease", required=True, help="Disease name, e.g. 'Metastatic melanoma'")
    parser.add_argument("--biofluid", required=True, help="Biofluid, e.g. 'plasma'")
    parser.add_argument("--data_dir", default=DATA_DIR, help="Directory with EV CSV/TSV tables")
    parser.add_argument("--out_json", required=True, help="Path to write JSON result")
    parser.add_argument("--out_md", required=True, help="Path to write markdown summary table")

    args = parser.parse_args()

    res = aggregate_biomarkers(args.disease, args.biofluid, data_dir=args.data_dir)
    os.makedirs(os.path.dirname(os.path.abspath(args.out_json)), exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(args.out_md)), exist_ok=True)

    # Helper function to convert NaN to None for JSON serialization
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


if __name__ == "__main__":
    run_cli()

