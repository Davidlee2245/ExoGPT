"""
Step 2 — Target & Epitope Curator
-----------------------------------

Selects nanobinder-suitable targets from Step 1 biomarker results.
Maps targets to structures (PDB/AlphaFold) and defines epitope patches
for binder design.

Input: Step 1 JSON output (biomarker list)
Output: Curated targets with epitope definitions (JSON + markdown table)
"""

from __future__ import annotations

import argparse
import json
import os
from dataclasses import dataclass, asdict
from statistics import mean
from typing import List, Dict, Any, Optional, Tuple

import numpy as np
import pandas as pd
import sys

try:
    import requests
    REQUESTS_AVAILABLE = True
except ImportError:  # pragma: no cover - env without requests
    requests = None
    REQUESTS_AVAILABLE = False

try:
    from Bio import PDB
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False
    PDB = None

try:
    from .structure_preprocessor import (
        StructureInfo,
        PreprocessedStructure,
        preprocess_structure,
        save_structure,
    )
    from .residue_features import (
        compute_residue_features,
        StructureFeatures,
    )
except ImportError:  # pragma: no cover - fallback for script execution
    sys.path.append(os.path.dirname(__file__))
    from structure_preprocessor import (
        StructureInfo,
        PreprocessedStructure,
        preprocess_structure,
        save_structure,
    )
    from residue_features import (
        compute_residue_features,
        StructureFeatures,
    )

DATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "data"))
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
STRUCTURE_CACHE_DIR = os.path.join(PROJECT_ROOT, "structures", "processed_step2")
AUTO_STRUCTURE_DIR = os.path.join(PROJECT_ROOT, "structures", "auto_downloaded")
ALPHAFOLD_URL_TEMPLATE = "https://alphafold.ebi.ac.uk/files/AF-{uniprot}-F1-model_v4.pdb"


@dataclass
class EpitopePatch:
    """Defines an epitope patch on a target structure."""
    epitope_id: str
    epitope_residues: List[int]
    epitope_description: str
    selection_method: str  # e.g., "surface_loop", "sasa_based", "manual"
    notes: Optional[str] = None


@dataclass
class CuratedTarget:
    """A curated target protein with structure and epitope information."""
    target_name: str  # e.g., "PD-L1"
    gene_symbol: str
    uniprot: Optional[str]
    structure_type: str  # "pdb" or "alphafold"
    model_path: str  # Path to PDB file or AlphaFold model
    chain_id: Optional[str] = None
    pdb_id: Optional[str] = None  # If structure_type is "pdb"
    domain: str = "ectodomain"  # e.g., "ectodomain", "full_length", "extracellular_region"
    domain_residues: Optional[Tuple[int, int]] = None  # (start, end) if domain is extracted
    epitopes: List[EpitopePatch] = None
    selection_rationale: str = ""
    surface_likelihood: Optional[float] = None
    biomarker_score: Optional[float] = None
    cleaned_structure_path: Optional[str] = None
    preprocessing_summary: Optional[str] = None
    surface_patch_stats: Optional[List[Dict[str, Any]]] = None
    
    def __post_init__(self):
        if self.epitopes is None:
            self.epitopes = []
        if self.surface_patch_stats is None:
            self.surface_patch_stats = []


def load_structure_mapping(mapping_path: str) -> pd.DataFrame:
    """
    Load gene_symbol/UniProt → structure mapping file.
    
    Expected columns:
    - gene_symbol: Gene name (e.g., "CD274")
    - uniprot: UniProt ID (e.g., "Q9NZQ7")
    - structure_type: "pdb" or "alphafold"
    - model_path: Path to structure file
    - pdb_id: PDB ID if structure_type is "pdb" (optional)
    - domain: Domain description (e.g., "ectodomain")
    - domain_residues: "start-end" format (e.g., "19-238") or empty
    - notes: Additional notes (optional)
    """
    if not os.path.exists(mapping_path):
        return pd.DataFrame(columns=[
            "gene_symbol", "uniprot", "structure_type", "model_path",
            "pdb_id", "domain", "domain_residues", "notes"
        ])
    
    try:
        # Read CSV with proper quote handling to handle commas in notes field
        df = pd.read_csv(mapping_path, sep=",", quotechar='"', skipinitialspace=True, engine="python")
    except Exception:
        # Fallback to auto-detect separator
        try:
            df = pd.read_csv(mapping_path, sep=None, engine="python", quotechar='"')
        except Exception:
            return pd.DataFrame(columns=[
                "gene_symbol", "uniprot", "structure_type", "model_path",
                "pdb_id", "domain", "domain_residues", "notes"
            ])
    
    # Remove empty rows
    df = df.dropna(how="all")
    
    # Normalize column names
    df.columns = [c.lower().strip() for c in df.columns]
    
    # Ensure required columns exist
    required = ["gene_symbol", "structure_type", "model_path"]
    for col in required:
        if col not in df.columns:
            df[col] = pd.NA
    
    # Filter out rows where required columns are missing
    if "gene_symbol" in df.columns:
        # Filter out rows where gene_symbol is NaN or empty string
        df = df[df["gene_symbol"].notna()]
        df = df[df["gene_symbol"].astype(str).str.strip() != ""]
    
    # Reset index to ensure clean indexing
    df = df.reset_index(drop=True)
    
    return df


def find_structure_for_target(
    gene_symbol: str,
    uniprot: Optional[str],
    mapping_df: pd.DataFrame
) -> Optional[Dict[str, Any]]:
    """
    Find structure mapping for a target protein.
    
    Returns dict with structure info or None if not found.
    """
    if mapping_df.empty:
        return None
    
    # Check if required columns exist
    if "gene_symbol" not in mapping_df.columns:
        return None
    
    # Normalize inputs
    gene_symbol_normalized = str(gene_symbol).strip().upper() if gene_symbol else ""
    uniprot_normalized = str(uniprot).strip().upper() if uniprot else ""
    
    # Debug: Print what we're looking for and what's in the DataFrame
    # This will help us see what's actually happening
    debug_info = {
        "searching_for_gene": gene_symbol_normalized,
        "searching_for_uniprot": uniprot_normalized,
        "mapping_df_columns": list(mapping_df.columns),
        "mapping_df_shape": mapping_df.shape,
        "mapping_genes": [],
        "mapping_uniprots": []
    }
    
    # Try to match by gene_symbol first - use simple iteration for reliability
    if gene_symbol_normalized:
        for idx in range(len(mapping_df)):
            row = mapping_df.iloc[idx]
            # Get gene_symbol value directly from Series
            try:
                mapping_gene_val = row["gene_symbol"]
                debug_info["mapping_genes"].append(str(mapping_gene_val))
                
                # Normalize and compare
                if pd.notna(mapping_gene_val):
                    mapping_gene = str(mapping_gene_val).strip().upper()
                    if mapping_gene == gene_symbol_normalized:
                        # Build result - use .iloc for reliable access
                        result = {
                            "gene_symbol": str(row["gene_symbol"]),
                            "uniprot": str(row["uniprot"]) if pd.notna(row["uniprot"]) else (uniprot or ""),
                            "structure_type": str(row["structure_type"]).lower(),
                            "model_path": str(row["model_path"]),
                            "pdb_id": str(row["pdb_id"]) if pd.notna(row["pdb_id"]) and str(row["pdb_id"]).strip() else None,
                            "domain": str(row["domain"]) if pd.notna(row["domain"]) else "ectodomain",
                            "domain_residues": _parse_domain_residues(row["domain_residues"] if pd.notna(row["domain_residues"]) else None),
                            "notes": str(row["notes"]) if pd.notna(row["notes"]) else None,
                        }
                        return result
            except (KeyError, IndexError) as e:
                # If column doesn't exist, skip this row
                continue
    
    # Try to match by UniProt if available
    if uniprot_normalized and "uniprot" in mapping_df.columns:
        for idx in range(len(mapping_df)):
            row = mapping_df.iloc[idx]
            # Get uniprot value directly from Series
            try:
                mapping_uniprot_val = row["uniprot"]
                debug_info["mapping_uniprots"].append(str(mapping_uniprot_val))
                
                # Normalize and compare
                if pd.notna(mapping_uniprot_val):
                    mapping_uniprot = str(mapping_uniprot_val).strip().upper()
                    if mapping_uniprot == uniprot_normalized:
                        # Build result - use .iloc for reliable access
                        result = {
                            "gene_symbol": str(row["gene_symbol"]),
                            "uniprot": str(row["uniprot"]),
                            "structure_type": str(row["structure_type"]).lower(),
                            "model_path": str(row["model_path"]),
                            "pdb_id": str(row["pdb_id"]) if pd.notna(row["pdb_id"]) and str(row["pdb_id"]).strip() else None,
                            "domain": str(row["domain"]) if pd.notna(row["domain"]) else "ectodomain",
                            "domain_residues": _parse_domain_residues(row["domain_residues"] if pd.notna(row["domain_residues"]) else None),
                            "notes": str(row["notes"]) if pd.notna(row["notes"]) else None,
                        }
                        return result
            except (KeyError, IndexError) as e:
                # If column doesn't exist, skip this row
                continue
    
    # If we get here, no match was found - store debug info for later
    # We'll return this in the skipped_biomarkers debug info
    return None


def _parse_domain_residues(domain_str: Any) -> Optional[Tuple[int, int]]:
    """Parse domain_residues string like '19-238' into (start, end) tuple."""
    if pd.isna(domain_str) or not domain_str:
        return None
    
    try:
        parts = str(domain_str).strip().split("-")
        if len(parts) == 2:
            return (int(parts[0]), int(parts[1]))
    except (ValueError, AttributeError):
        pass
    
    return None


def _resolve_model_path(model_path: str, mapping_file: str) -> str:
    """
    Resolve structure path that may be relative to different anchors.
    """
    candidate_paths = []
    
    if model_path:
        candidate_paths.append(os.path.abspath(model_path))
    
    mapping_dir = os.path.dirname(os.path.abspath(mapping_file)) if mapping_file else None
    if mapping_dir and model_path:
        candidate_paths.append(os.path.abspath(os.path.join(mapping_dir, model_path)))
        candidate_paths.append(os.path.abspath(os.path.join(mapping_dir, "..", model_path)))
    
    # Relative to project root
    if model_path:
        candidate_paths.append(os.path.abspath(os.path.join(PROJECT_ROOT, model_path.lstrip("./"))))
    
    for path in candidate_paths:
        if os.path.exists(path):
            return path
    
    # Fallback to absolute interpretation
    return os.path.abspath(model_path)


def _download_alphafold_model(
    gene_symbol: str,
    uniprot: str,
) -> Optional[Dict[str, Any]]:
    """
    Try to download an AlphaFold model for a UniProt ID.
    
    This uses the public AlphaFold DB URL template and saves the PDB under
    structures/auto_downloaded/. If download fails or requests is unavailable,
    returns None.
    """
    if not REQUESTS_AVAILABLE or not uniprot:
        return None
    
    os.makedirs(AUTO_STRUCTURE_DIR, exist_ok=True)
    uniprot_norm = str(uniprot).strip().upper()
    outfile = os.path.join(AUTO_STRUCTURE_DIR, f"{uniprot_norm}_alphafold.pdb")
    
    # Reuse if already downloaded
    if os.path.exists(outfile):
        return {
            "gene_symbol": gene_symbol,
            "uniprot": uniprot_norm,
            "structure_type": "alphafold",
            "model_path": outfile,
            "pdb_id": None,
            "domain": "full_length",
            "domain_residues": None,
            "notes": "Auto-downloaded AlphaFold model",
        }
    
    url = ALPHAFOLD_URL_TEMPLATE.format(uniprot=uniprot_norm)
    try:
        resp = requests.get(url, timeout=30)
        if resp.status_code != 200 or "ATOM" not in resp.text:
            return None
        with open(outfile, "w") as f:
            f.write(resp.text)
    except Exception:
        return None
    
    return {
        "gene_symbol": gene_symbol,
        "uniprot": uniprot_norm,
        "structure_type": "alphafold",
        "model_path": outfile,
        "pdb_id": None,
        "domain": "full_length",
        "domain_residues": None,
        "notes": "Auto-downloaded AlphaFold model",
    }


def _infer_chain_id(structure_file: str) -> str:
    """
    Infer a reasonable default chain ID by selecting the first polypeptide chain.
    """
    if not BIOPYTHON_AVAILABLE:
        return "A"
    
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("infer", structure_file)
    
    for chain in structure.get_chains():
        for residue in chain:
            res_id = residue.id
            if isinstance(res_id, tuple) and len(res_id) >= 1 and res_id[0] == " ":
                return chain.id
    
    return "A"


def _remove_non_protein_residues(structure: Any, chain_id: str) -> Tuple[Any, int]:
    """
    Remove ligands/waters from the selected chain to keep only protein residues.
    
    Returns:
        (cleaned_structure, removed_count)
    """
    if not BIOPYTHON_AVAILABLE:
        return structure, 0
    
    removed = 0
    for chain in structure.get_chains():
        if chain.id != chain_id:
            continue
        to_remove = []
        for residue in chain:
            hetflag = residue.id[0] if isinstance(residue.id, tuple) else residue.id
            if hetflag is None:
                continue
            if hetflag.strip():  # Non-empty = hetero residue (ligand, water, etc.)
                to_remove.append(residue.id)
        for residue_id in to_remove:
            chain.detach_child(residue_id)
            removed += 1
    
    return structure, removed


def _extract_ca_coordinates(structure: Any, chain_id: str) -> Dict[int, np.ndarray]:
    """Extract CA coordinates for residues in chain."""
    coords = {}
    if not BIOPYTHON_AVAILABLE:
        return coords
    
    for chain in structure.get_chains():
        if chain.id != chain_id:
            continue
        for residue in chain:
            res_id = residue.id
            if isinstance(res_id, tuple) and len(res_id) >= 2:
                resnum = res_id[1]
                if isinstance(resnum, int):
                    try:
                        ca = residue["CA"]
                        coords[resnum] = np.array(ca.get_coord(), dtype=float)
                    except KeyError:
                        continue
    return coords


def _cluster_surface_residues(
    surface_residues: List[int],
    ca_coords: Dict[int, np.ndarray],
    distance_cutoff: float = 7.0
) -> List[List[int]]:
    """
    Cluster surface residues using a simple Cα distance graph.
    """
    if not surface_residues:
        return []
    
    nodes = [res for res in surface_residues if res in ca_coords]
    if len(nodes) < 2:
        return [nodes] if nodes else []
    
    adjacency: Dict[int, List[int]] = {res: [] for res in nodes}
    for i in range(len(nodes)):
        for j in range(i + 1, len(nodes)):
            res_i = nodes[i]
            res_j = nodes[j]
            coord_i = ca_coords[res_i]
            coord_j = ca_coords[res_j]
            distance = float(np.linalg.norm(coord_i - coord_j))
            if distance <= distance_cutoff:
                adjacency[res_i].append(res_j)
                adjacency[res_j].append(res_i)
    
    visited = set()
    components: List[List[int]] = []
    
    for node in nodes:
        if node in visited:
            continue
        stack = [node]
        component = []
        while stack:
            current = stack.pop()
            if current in visited:
                continue
            visited.add(current)
            component.append(current)
            stack.extend(adjacency.get(current, []))
        if component:
            components.append(sorted(component))
    
    return components


def _select_surface_patches(
    components: List[List[int]],
    residue_features: Dict[int, Any],
    size_range: Tuple[int, int] = (5, 25)
) -> List[Dict[str, Any]]:
    """Select candidate epitope patches from connected components."""
    min_size, max_size = size_range
    patch_stats: List[Dict[str, Any]] = []
    
    for comp in components:
        feature_subset = [residue_features[res] for res in comp if res in residue_features]
        if len(feature_subset) < min_size or len(feature_subset) > max_size:
            continue
        avg_sasa = mean([float(rf.sasa) for rf in feature_subset]) if feature_subset else 0.0
        mean_b = mean([float(rf.bfactor) for rf in feature_subset]) if feature_subset else 0.0
        loop_fraction = (
            sum(1 for rf in feature_subset if rf.secondary_structure_category == "loop") / len(feature_subset)
            if feature_subset else 0.0
        )
        patch_stats.append({
            "residues": comp,
            "size": len(feature_subset),
            "avg_sasa": avg_sasa,
            "mean_bfactor": mean_b,
            "loop_fraction": loop_fraction,
        })
    
    patch_stats.sort(
        key=lambda p: (p["avg_sasa"], p["loop_fraction"], -p["mean_bfactor"]),
        reverse=True
    )
    return patch_stats


def _analyze_structure_for_target(
    target_id: str,
    structure_info: Dict[str, Any],
    structure_mapping_path: str,
    sasa_threshold: float = 30.0,
    graph_distance_cutoff: float = 7.0
) -> Optional[Dict[str, Any]]:
    """
    Run structure preprocessing + SASA/patch identification for a target.
    """
    if not BIOPYTHON_AVAILABLE:
        raise ImportError("BioPython is required for structure analysis but is not installed.")
    
    resolved_model_path = _resolve_model_path(structure_info["model_path"], structure_mapping_path)
    if not os.path.exists(resolved_model_path):
        raise FileNotFoundError(f"Structure file not found at {resolved_model_path}")
    
    chain_id = structure_info.get("chain_id")
    if not chain_id or chain_id == "nan":
        chain_id = _infer_chain_id(resolved_model_path)
    
    os.makedirs(STRUCTURE_CACHE_DIR, exist_ok=True)
    cleaned_filename = f"{target_id.lower()}_{chain_id.lower()}_{structure_info['structure_type']}_cleaned.pdb"
    cleaned_path = os.path.join(STRUCTURE_CACHE_DIR, cleaned_filename)
    
    structure_meta = StructureInfo(
        target_id=target_id,
        structure_file=resolved_model_path,
        structure_type=structure_info["structure_type"],
        chain_id=chain_id,
        domain_type=structure_info.get("domain", "ectodomain"),
        domain_residues=structure_info.get("domain_residues"),
        is_transmembrane="full" in str(structure_info.get("domain", "")).lower()
    )
    
    preprocessed: PreprocessedStructure = preprocess_structure(
        structure_meta,
        output_file=cleaned_path
    )
    
    protein_only_structure, removed_residues = _remove_non_protein_residues(preprocessed.structure, chain_id)
    save_structure(protein_only_structure, cleaned_path)
    
    structure_features: StructureFeatures = compute_residue_features(
        structure_file=cleaned_path,
        structure=protein_only_structure,
        chain_id=chain_id,
        sasa_threshold=sasa_threshold,
        structure_type=structure_info["structure_type"]
    )
    
    residue_feature_map = {rf.residue_number: rf for rf in structure_features.residues}
    surface_residues = [
        rf.residue_number
        for rf in structure_features.residues
        if rf.is_surface and not rf.is_low_confidence
    ]
    
    ca_coords = _extract_ca_coordinates(protein_only_structure, chain_id)
    components = _cluster_surface_residues(surface_residues, ca_coords, distance_cutoff=graph_distance_cutoff)
    patch_stats = _select_surface_patches(components, residue_feature_map)
    
    analysis = {
        "chain_id": chain_id,
        "cleaned_structure_path": cleaned_path,
        "structure_features": {
            "sasa_method": structure_features.sasa_method,
            "total_residues": structure_features.total_residues,
            "surface_residues": structure_features.surface_residues,
            "low_confidence_residues": structure_features.low_confidence_residues,
        },
        "surface_residue_count": len(surface_residues),
        "removed_ligands": removed_residues,
        "patch_stats": patch_stats,
        "candidate_patch_count": len(components),
        "selected_patch_count": len(patch_stats),
        "sasa_threshold": sasa_threshold,
        "graph_distance_cutoff": graph_distance_cutoff,
    }
    
    return analysis


def is_nanobinder_suitable(biomarker: Dict[str, Any]) -> bool:
    """
    Determine if a biomarker is suitable for nanobinder targeting.
    
    Criteria:
    - Must be a protein (not miRNA, lipid, etc.)
    - Should have surface_likelihood > 0.3 OR be a known surface protein
    - Should have reasonable evidence (at least 1 study)
    - Prefer targets with UniProt ID (for structure mapping)
    """
    analyte_type = str(biomarker.get("analyte_type", "")).lower()
    
    # Must be a protein
    if "protein" not in analyte_type:
        return False
    
    # Should have at least some evidence
    num_studies = biomarker.get("num_studies", 0)
    if num_studies < 1:
        return False
    
    # Surface likelihood check
    surface_likelihood = biomarker.get("surface_likelihood")
    if surface_likelihood is not None and surface_likelihood > 0.3:
        return True
    
    # Known surface protein keywords
    gene_symbol = str(biomarker.get("gene_symbol", "")).upper()
    surface_keywords = [
        "CD274", "PD-L1", "PDCD1", "PD-1", "CTLA4", "CD63", "CD81", "CD9",
        "CD44", "CD151", "ITGA", "ITGB", "ICAM", "VCAM", "TETRASPANIN"
    ]
    for keyword in surface_keywords:
        if keyword in gene_symbol:
            return True
    
    # If has UniProt and decent score, include it
    if biomarker.get("uniprot") and biomarker.get("score", 0) > 5.0:
        return True
    
    return False


def define_epitope_patches(
    target: CuratedTarget,
    structure_info: Dict[str, Any],
    structure_analysis: Optional[Dict[str, Any]] = None
) -> List[EpitopePatch]:
    """
    Define epitope patches for a target using SASA-derived clusters when available.
    """
    if structure_analysis and structure_analysis.get("patch_stats"):
        patches = []
        for idx, patch in enumerate(structure_analysis["patch_stats"], start=1):
            epitope_id = f"{target.gene_symbol.lower()}_patch_{idx:02d}"
            desc = (
                f"SASA mean {patch['avg_sasa']:.1f} Å², "
                f"{patch['loop_fraction'] * 100:.0f}% loop residues"
            )
            notes = (
                f"Surface patch detected by SASA>{structure_analysis['sasa_threshold']} Å² "
                f"and Cα graph cutoff {structure_analysis['graph_distance_cutoff']} Å; "
                f"method={structure_analysis['structure_features']['sasa_method']}"
            )
            patches.append(EpitopePatch(
                epitope_id=epitope_id,
                epitope_residues=patch["residues"],
                epitope_description=desc,
                selection_method="sasa_surface_patch",
                notes=notes
            ))
        return patches
    
    return _fallback_epitope_templates(target, structure_info)


def _fallback_epitope_templates(
    target: CuratedTarget,
    structure_info: Dict[str, Any]
) -> List[EpitopePatch]:
    """Legacy placeholder epitopes when structure analysis is unavailable."""
    epitopes = []
    gene_upper = target.gene_symbol.upper()
    
    if gene_upper in ["CD274", "PD-L1"]:
        epitopes.append(EpitopePatch(
            epitope_id=f"{target.gene_symbol.lower()}_epitope_01",
            epitope_residues=[52, 54, 55, 76, 77, 120, 121, 122],
            epitope_description="Surface-exposed loop patch on IgV domain, away from PD-1 binding interface",
            selection_method="literature_based",
            notes="Based on known antibody epitopes. In practice, would use SASA > 30 Å² and loop region analysis."
        ))
        epitopes.append(EpitopePatch(
            epitope_id=f"{target.gene_symbol.lower()}_epitope_02",
            epitope_residues=[19, 20, 21, 22, 25, 26, 27, 28],
            epitope_description="N-terminal region, highly solvent-exposed",
            selection_method="literature_based",
            notes="Alternative epitope site. Would verify with SASA calculation."
        ))
    elif gene_upper in ["CD63", "CD81", "CD9"]:
        epitopes.append(EpitopePatch(
            epitope_id=f"{target.gene_symbol.lower()}_epitope_01",
            epitope_residues=[100, 101, 102, 103, 104, 105, 106, 107, 108],
            epitope_description="Large extracellular loop (ECL2), highly exposed",
            selection_method="domain_based",
            notes="Tetraspanins have characteristic large ECLs. Would use domain annotation + SASA to refine."
        ))
    elif gene_upper in ["CD44"]:
        epitopes.append(EpitopePatch(
            epitope_id=f"{target.gene_symbol.lower()}_epitope_01",
            epitope_residues=[50, 51, 52, 53, 54, 55, 56, 57],
            epitope_description="Extracellular domain surface patch",
            selection_method="domain_based",
            notes="Would use SASA and loop identification to select optimal patch."
        ))
    else:
        domain_res = structure_info.get("domain_residues")
        if domain_res:
            start, end = domain_res
            mid_point = (start + end) // 2
            epitopes.append(EpitopePatch(
                epitope_id=f"{target.gene_symbol.lower()}_epitope_01",
                epitope_residues=list(range(mid_point - 3, mid_point + 5)),
                epitope_description="Surface-exposed region (placeholder - requires SASA analysis)",
                selection_method="template",
                notes=f"PLACEHOLDER: Would calculate SASA for residues {start}-{end} and select highest exposure patch."
            ))
        else:
            epitopes.append(EpitopePatch(
                epitope_id=f"{target.gene_symbol.lower()}_epitope_01",
                epitope_residues=[50, 51, 52, 53, 54, 55, 56, 57],
                epitope_description="Surface-exposed region (placeholder - requires structure analysis)",
                selection_method="template",
                notes="PLACEHOLDER: Requires loading structure file and computing SASA."
            ))
    
    return epitopes


def curate_targets(
    step1_json_path: str,
    structure_mapping_path: str,
    user_specified_target: Optional[str] = None,
    max_targets: int = 5
) -> Dict[str, Any]:
    """
    Main curation function: Select targets and define epitopes.
    
    Args:
        step1_json_path: Path to Step 1 JSON output
        structure_mapping_path: Path to structure mapping CSV
        user_specified_target: Optional specific target (e.g., "CD274" or "PD-L1")
        max_targets: Maximum number of targets to curate
    
    Returns:
        JSON-serializable dict with curated targets
    """
    # Load Step 1 results
    with open(step1_json_path, "r") as f:
        step1_data = json.load(f)
    
    biomarkers = step1_data.get("ev_biomarkers", [])
    disease = step1_data.get("disease", "Unknown")
    biofluid = step1_data.get("biofluid", "Unknown")
    
    # Load structure mapping
    mapping_df = load_structure_mapping(structure_mapping_path)
    
    # If user specified a target, search ALL biomarkers first (not just suitable ones)
    # This ensures we can find and curate the specified target even if it doesn't pass
    # all the suitability filters
    if user_specified_target:
        target_upper = user_specified_target.upper()
        
        # Handle common aliases
        target_aliases = [target_upper]
        if target_upper in ["CD274", "PD-L1", "PDL1"]:
            target_aliases.extend(["CD274", "PD-L1", "PDL1"])
        elif target_upper in ["PDCD1", "PD-1", "PD1"]:
            target_aliases.extend(["PDCD1", "PD-1", "PD1"])
        
        # Search all biomarkers for the specified target (check all aliases)
        specified_biomarkers = [
            bm for bm in biomarkers
            if any(
                alias in str(bm.get("gene_symbol", "")).upper()
                or alias in str(bm.get("uniprot", "")).upper()
                for alias in target_aliases
            )
        ]
        
        # Filter for nanobinder-suitable targets
        suitable_biomarkers = [
            bm for bm in biomarkers
            if is_nanobinder_suitable(bm)
        ]
        
        # Prioritize specified targets (even if not in suitable list)
        # But still check if they're proteins
        specified_proteins = [
            bm for bm in specified_biomarkers
            if "protein" in str(bm.get("analyte_type", "")).lower()
        ]
        
        # Combine: specified proteins first, then other suitable biomarkers
        if specified_proteins:
            suitable_biomarkers = specified_proteins + [
                bm for bm in suitable_biomarkers
                if bm not in specified_proteins
            ]
    else:
        # Filter for nanobinder-suitable targets
        suitable_biomarkers = [
            bm for bm in biomarkers
            if is_nanobinder_suitable(bm)
        ]
    
    # Sort by score (highest first)
    suitable_biomarkers.sort(key=lambda b: b.get("score", 0), reverse=True)
    
    # Limit to max_targets
    selected_biomarkers = suitable_biomarkers[:max_targets]
    
    # Curate each target
    curated_targets: List[CuratedTarget] = []
    
    # Track which biomarkers were skipped and why
    skipped_biomarkers = []
    
    # Debug: Check what's actually in the mapping DataFrame
    mapping_genes = []
    mapping_uniprots = []
    if not mapping_df.empty and "gene_symbol" in mapping_df.columns:
        mapping_genes = [str(g).strip().upper() for g in mapping_df["gene_symbol"].dropna().tolist()]
    if not mapping_df.empty and "uniprot" in mapping_df.columns:
        mapping_uniprots = [str(u).strip().upper() for u in mapping_df["uniprot"].dropna().tolist()]
    
    structure_processing_debug = []
    
    for bm in selected_biomarkers:
        gene_symbol = str(bm.get("gene_symbol", ""))
        uniprot = bm.get("uniprot")
        
        # Find structure mapping
        structure_info = find_structure_for_target(gene_symbol, uniprot, mapping_df)
        
        if not structure_info:
            # Optionally try to auto-download an AlphaFold model if UniProt is available
            auto_download_attempted = False
            auto_download_succeeded = False
            auto_structure_info: Optional[Dict[str, Any]] = None
            
            if uniprot:
                auto_download_attempted = True
                auto_structure_info = _download_alphafold_model(gene_symbol, uniprot)
                if auto_structure_info:
                    auto_download_succeeded = True
                    structure_info = auto_structure_info
            
        if not structure_info:
            # Track why this was skipped with more detail (no mapping and/or auto-download failed)
            gene_normalized = str(gene_symbol).strip().upper()
            uniprot_normalized = str(uniprot).strip().upper() if uniprot else ""
            in_genes = gene_normalized in mapping_genes
            in_uniprots = uniprot_normalized in mapping_uniprots if uniprot_normalized else False
            
            # Get actual values from DataFrame for debugging
            actual_genes_in_df = []
            actual_uniprots_in_df = []
            if not mapping_df.empty:
                try:
                    for idx in range(len(mapping_df)):
                        row = mapping_df.iloc[idx]
                        actual_genes_in_df.append(str(row["gene_symbol"]).strip().upper() if pd.notna(row["gene_symbol"]) else "N/A")
                        if "uniprot" in mapping_df.columns:
                            actual_uniprots_in_df.append(str(row["uniprot"]).strip().upper() if pd.notna(row["uniprot"]) else "N/A")
                except Exception:
                    pass
            
            skipped_biomarkers.append({
                "gene_symbol": gene_symbol,
                "uniprot": uniprot,
                "gene_normalized": gene_normalized,
                "uniprot_normalized": uniprot_normalized,
                "reason": "no_structure_mapping",
                "mapping_file_exists": os.path.exists(structure_mapping_path) if structure_mapping_path else False,
                "mapping_df_empty": mapping_df.empty,
                "mapping_df_rows": len(mapping_df),
                "mapping_df_columns": list(mapping_df.columns) if not mapping_df.empty else [],
                "gene_in_mapping": in_genes,
                "uniprot_in_mapping": in_uniprots,
                "mapping_genes_sample": mapping_genes[:5] if mapping_genes else [],
                "mapping_uniprots_sample": mapping_uniprots[:5] if mapping_uniprots else [],
                "actual_genes_in_df": actual_genes_in_df[:10],
                "actual_uniprots_in_df": actual_uniprots_in_df[:10],
                "auto_download_attempted": auto_download_attempted,
                "auto_download_succeeded": auto_download_succeeded,
            })
            continue
        
        # Analyze structure (cleaning + SASA + clustering)
        analysis_error = None
        structure_analysis: Optional[Dict[str, Any]] = None
        try:
            structure_analysis = _analyze_structure_for_target(
                gene_symbol or structure_info.get("gene_symbol") or structure_info.get("uniprot") or "target",
                structure_info,
                structure_mapping_path
            )
        except Exception as exc:
            analysis_error = str(exc)
        
        if analysis_error:
            skipped_biomarkers.append({
                "gene_symbol": gene_symbol,
                "uniprot": uniprot,
                "reason": "structure_processing_failed",
                "error": analysis_error
            })
            continue
        
        # Create curated target
        # Format surface likelihood for display
        surf_likelihood = bm.get('surface_likelihood')
        surf_likelihood_str = f"{surf_likelihood:.2f}" if surf_likelihood is not None else "N/A"
        
        target = CuratedTarget(
            target_name=gene_symbol,  # Could be enhanced with protein name lookup
            gene_symbol=gene_symbol,
            uniprot=structure_info.get("uniprot") or uniprot,
            structure_type=structure_info["structure_type"],
            model_path=structure_info["model_path"],
            pdb_id=structure_info.get("pdb_id"),
            domain=structure_info.get("domain", "ectodomain"),
            domain_residues=structure_info.get("domain_residues"),
            selection_rationale=f"Selected from Step 1 biomarkers. Score: {bm.get('score', 0):.2f}, "
                               f"Surface likelihood: {surf_likelihood_str}, "
                               f"Evidence: {bm.get('evidence_level', 'unknown')}",
            surface_likelihood=bm.get("surface_likelihood"),
            biomarker_score=bm.get("score")
        )
        
        if structure_analysis:
            target.chain_id = structure_analysis["chain_id"]
            target.cleaned_structure_path = structure_analysis["cleaned_structure_path"]
            sf = structure_analysis["structure_features"]
            target.preprocessing_summary = (
                f"{sf['total_residues']} residues | {sf['surface_residues']} surface "
                f"| SASA via {sf['sasa_method']}"
            )
            target.surface_patch_stats = structure_analysis["patch_stats"]
        
        # Define epitopes
        epitopes = define_epitope_patches(target, structure_info, structure_analysis)
        if not epitopes:
            epitopes = define_epitope_patches(target, structure_info, None)
        target.epitopes = epitopes
        
        curated_targets.append(target)
        
        structure_processing_debug.append({
            "gene_symbol": gene_symbol,
            "chain_id": target.chain_id,
            "cleaned_structure_path": target.cleaned_structure_path,
            "sasa_method": structure_analysis["structure_features"]["sasa_method"] if structure_analysis else None,
            "surface_residues": structure_analysis["surface_residue_count"] if structure_analysis else None,
            "removed_ligands": structure_analysis["removed_ligands"] if structure_analysis else None,
            "candidate_patch_count": structure_analysis["candidate_patch_count"] if structure_analysis else 0,
            "selected_patches": structure_analysis["selected_patch_count"] if structure_analysis else 0,
            "epitopes_emitted": len(epitopes),
            "sasa_epitopes": len(structure_analysis["patch_stats"]) if structure_analysis else 0,
        })
    
    # Convert to dict for JSON serialization
    result = {
        "disease": disease,
        "biofluid": biofluid,
        "user_specified_target": user_specified_target,
        "curated_targets": [_curated_target_to_dict(t) for t in curated_targets],
        "structure_mapping_file": structure_mapping_path,
        "num_targets": len(curated_targets),
        "methodology_notes": _get_methodology_notes(),
        "debug_info": {
            "total_biomarkers": len(biomarkers),
            "suitable_biomarkers": len(suitable_biomarkers),
            "selected_biomarkers": len(selected_biomarkers),
            "structure_mapping_file": structure_mapping_path,
            "structure_mapping_exists": os.path.exists(structure_mapping_path),
            "structure_mapping_loaded": not mapping_df.empty,
            "structure_mapping_rows": len(mapping_df),
            "curated_targets": len(curated_targets),
            "skipped_no_structure": len(selected_biomarkers) - len(curated_targets),
            "skipped_biomarkers": skipped_biomarkers,
            "structure_processing_details": structure_processing_debug,
            "selected_biomarker_details": [
                {
                    "gene_symbol": str(bm.get("gene_symbol", "")),
                    "uniprot": bm.get("uniprot"),
                    "analyte_type": bm.get("analyte_type"),
                    "score": bm.get("score")
                }
                for bm in selected_biomarkers
            ]
        }
    }
    
    return result


def _curated_target_to_dict(target: CuratedTarget) -> Dict[str, Any]:
    """Convert CuratedTarget dataclass to dict for JSON serialization."""
    d = asdict(target)
    # Convert domain_residues tuple to list or None
    if d.get("domain_residues"):
        d["domain_residues"] = list(d["domain_residues"])
    else:
        d["domain_residues"] = None
    # Convert epitopes
    d["epitopes"] = [asdict(ep) for ep in target.epitopes]
    return d


def _get_methodology_notes() -> str:
    """Return notes on how epitope selection would work with 3D tools."""
    return """
Epitope Selection Methodology (for future implementation with 3D analysis tools):

1. Structure Loading:
   - Load PDB file or AlphaFold model using BioPython or PyMOL
   - Extract chain of interest (if multi-chain structure)

2. Domain Extraction (if needed):
   - If full-length structure, extract extracellular domain using:
     - Signal peptide prediction (e.g., SignalP)
     - Transmembrane region prediction (e.g., TMHMM)
     - Domain annotation from UniProt or Pfam

3. Solvent Accessible Surface Area (SASA) Calculation:
   - Use tools like:
     - BioPython's DSSP or FreeSASA
     - PyMOL's get_area command
     - NACCESS
   - Calculate SASA for each residue
   - Filter residues with SASA > threshold (e.g., 0.5 or 20 Å²)

4. Epitope Patch Selection:
   - Group solvent-exposed residues into contiguous patches
   - Prefer patches in loop regions (high B-factor or secondary structure annotation)
   - Avoid:
     - Glycosylation sites (check UniProt PTM annotations)
     - Known functional interfaces (e.g., ligand binding sites)
     - Disulfide bonds (if known)
   - Select 1-3 patches per target, each with 5-15 residues

5. Validation:
   - Visualize patches in PyMOL or ChimeraX
   - Check for steric accessibility
   - Verify patch is not occluded by other domains/chains

Current implementation uses literature-based and template-based epitopes.
Replace with computed epitopes once 3D analysis tools are integrated.
"""


def curated_targets_to_markdown(curated_targets: List[Dict[str, Any]]) -> str:
    """Generate markdown table summarizing curated targets and epitopes."""
    lines = []
    lines.append("# Curated Targets & Epitopes\n")
    lines.append("|Target|Gene|UniProt|Structure|Domain|Epitopes|Epitope Residues|")
    lines.append("|---|---|---|---|---|---|---|")
    
    for target in curated_targets:
        target_name = target.get("target_name", "—")
        gene = target.get("gene_symbol", "—")
        uniprot = target.get("uniprot") or "—"
        structure_type = target.get("structure_type", "—")
        pdb_id = target.get("pdb_id", "")
        structure = f"{structure_type.upper()}"
        if pdb_id:
            structure += f" ({pdb_id})"
        domain = target.get("domain", "—")
        
        epitopes = target.get("epitopes", [])
        if epitopes:
            for idx, ep in enumerate(epitopes, 1):
                ep_id = ep.get("epitope_id", f"epitope_{idx}")
                ep_res = ", ".join(map(str, ep.get("epitope_residues", [])))
                ep_desc = ep.get("epitope_description", "")
                
                if idx == 1:
                    # First row includes target info
                    lines.append(
                        f"|{target_name}|{gene}|{uniprot}|{structure}|{domain}|"
                        f"{ep_id}<br/>{ep_desc}|{ep_res}|"
                    )
                else:
                    # Subsequent rows for additional epitopes
                    lines.append(
                        f"|||{uniprot}|||{ep_id}<br/>{ep_desc}|{ep_res}|"
                    )
        else:
            lines.append(
                f"|{target_name}|{gene}|{uniprot}|{structure}|{domain}|None|—|"
            )
    
    return "\n".join(lines)


def run_cli() -> None:
    """CLI entrypoint for Step 2."""
    parser = argparse.ArgumentParser(description="Step 2 — Target & Epitope Curator")
    parser.add_argument("--step1_json", required=True, help="Path to Step 1 JSON output")
    parser.add_argument(
        "--structure_mapping",
        default=os.path.join(DATA_DIR, "structure_mapping.csv"),
        help="Path to structure mapping CSV file"
    )
    parser.add_argument(
        "--target",
        default=None,
        help="Optional: Specify a particular target (e.g., 'CD274' or 'PD-L1')"
    )
    parser.add_argument(
        "--max_targets",
        type=int,
        default=5,
        help="Maximum number of targets to curate (default: 5)"
    )
    parser.add_argument("--out_json", required=True, help="Path to write JSON result")
    parser.add_argument("--out_md", required=True, help="Path to write markdown summary table")
    
    args = parser.parse_args()
    
    result = curate_targets(
        step1_json_path=args.step1_json,
        structure_mapping_path=args.structure_mapping,
        user_specified_target=args.target,
        max_targets=args.max_targets
    )
    
    # Ensure output directories exist
    os.makedirs(os.path.dirname(os.path.abspath(args.out_json)), exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(args.out_md)), exist_ok=True)
    
    # Write JSON
    with open(args.out_json, "w") as f:
        json.dump(result, f, indent=2)
    
    # Write markdown
    md = curated_targets_to_markdown(result["curated_targets"])
    with open(args.out_md, "w") as f:
        f.write(md + "\n")
    
    print(f"✓ Curated {result['num_targets']} targets")
    print(f"✓ Wrote JSON to: {args.out_json}")
    print(f"✓ Wrote markdown to: {args.out_md}")


if __name__ == "__main__":
    run_cli()

