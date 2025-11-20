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
from typing import List, Dict, Any, Optional, Tuple

import pandas as pd

DATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "data"))


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
    pdb_id: Optional[str] = None  # If structure_type is "pdb"
    domain: str = "ectodomain"  # e.g., "ectodomain", "full_length", "extracellular_region"
    domain_residues: Optional[Tuple[int, int]] = None  # (start, end) if domain is extracted
    epitopes: List[EpitopePatch] = None
    selection_rationale: str = ""
    surface_likelihood: Optional[float] = None
    biomarker_score: Optional[float] = None
    
    def __post_init__(self):
        if self.epitopes is None:
            self.epitopes = []


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
    structure_info: Dict[str, Any]
) -> List[EpitopePatch]:
    """
    Define epitope patches for a target.
    
    NOTE: This is a template/placeholder function. In a real implementation,
    you would:
    1. Load the structure file (PDB or AlphaFold)
    2. Calculate Solvent Accessible Surface Area (SASA) for each residue
    3. Identify surface-exposed loops/regions
    4. Avoid glycosylation sites (if known)
    5. Select patches that are:
       - Solvent-exposed (high SASA)
       - In flexible loop regions (if available from B-factors)
       - Away from known functional sites (if annotated)
    
    For now, we provide example epitopes based on known structures and
    document the methodology for future implementation.
    """
    epitopes = []
    
    # Example epitope definitions for common targets
    # These are based on literature and known structures
    gene_upper = target.gene_symbol.upper()
    
    if gene_upper in ["CD274", "PD-L1"]:
        # PD-L1 ectodomain epitope (based on PDB 3BIK, 5JDR, etc.)
        # This is a surface-exposed loop region on the IgV domain
        epitopes.append(EpitopePatch(
            epitope_id=f"{target.gene_symbol.lower()}_epitope_01",
            epitope_residues=[52, 54, 55, 76, 77, 120, 121, 122],
            epitope_description="Surface-exposed loop patch on IgV domain, away from PD-1 binding interface",
            selection_method="literature_based",
            notes="Based on known antibody epitopes. In practice, would use SASA > 0.5 and loop region analysis."
        ))
        # Alternative epitope
        epitopes.append(EpitopePatch(
            epitope_id=f"{target.gene_symbol.lower()}_epitope_02",
            epitope_residues=[19, 20, 21, 22, 25, 26, 27, 28],
            epitope_description="N-terminal region, highly solvent-exposed",
            selection_method="literature_based",
            notes="Alternative epitope site. Would verify with SASA calculation."
        ))
    
    elif gene_upper in ["CD63", "CD81", "CD9"]:
        # Tetraspanin proteins - typically have large extracellular loops
        # These are placeholder residues; real analysis would identify the largest ECL
        epitopes.append(EpitopePatch(
            epitope_id=f"{target.gene_symbol.lower()}_epitope_01",
            epitope_residues=[100, 101, 102, 103, 104, 105, 106, 107, 108],
            epitope_description="Large extracellular loop (ECL2), highly exposed",
            selection_method="domain_based",
            notes="Tetraspanins have characteristic large ECLs. Would use domain annotation + SASA to refine."
        ))
    
    elif gene_upper in ["CD44"]:
        # CD44 has a large extracellular domain
        epitopes.append(EpitopePatch(
            epitope_id=f"{target.gene_symbol.lower()}_epitope_01",
            epitope_residues=[50, 51, 52, 53, 54, 55, 56, 57],
            epitope_description="Extracellular domain surface patch",
            selection_method="domain_based",
            notes="Would use SASA and loop identification to select optimal patch."
        ))
    
    else:
        # Generic epitope template for unknown targets
        # In practice, this would be computed from structure analysis
        domain_res = structure_info.get("domain_residues")
        if domain_res:
            start, end = domain_res
            # Placeholder: select middle region (would be refined with SASA)
            mid_point = (start + end) // 2
            epitopes.append(EpitopePatch(
                epitope_id=f"{target.gene_symbol.lower()}_epitope_01",
                epitope_residues=list(range(mid_point - 3, mid_point + 5)),
                epitope_description="Surface-exposed region (placeholder - requires SASA analysis)",
                selection_method="template",
                notes=f"PLACEHOLDER: Would calculate SASA for residues {start}-{end} and select highest exposure patch. "
                      f"Should avoid glycosylation sites and functional interfaces."
            ))
        else:
            # No domain info - very generic placeholder
            epitopes.append(EpitopePatch(
                epitope_id=f"{target.gene_symbol.lower()}_epitope_01",
                epitope_residues=[50, 51, 52, 53, 54, 55, 56, 57],
                epitope_description="Surface-exposed region (placeholder - requires structure analysis)",
                selection_method="template",
                notes="PLACEHOLDER: Requires loading structure file and computing SASA. "
                      "Would select residues with SASA > 0.5 in loop regions, avoiding known functional sites."
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
    
    for bm in selected_biomarkers:
        gene_symbol = str(bm.get("gene_symbol", ""))
        uniprot = bm.get("uniprot")
        
        # Find structure mapping
        structure_info = find_structure_for_target(gene_symbol, uniprot, mapping_df)
        
        if not structure_info:
            # Track why this was skipped with more detail
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
                "actual_uniprots_in_df": actual_uniprots_in_df[:10]
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
        
        # Define epitopes
        epitopes = define_epitope_patches(target, structure_info)
        target.epitopes = epitopes
        
        curated_targets.append(target)
    
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

