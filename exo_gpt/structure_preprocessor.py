"""
Structure Preprocessing Module
==============================

Pre-processes protein structures (PDB/AlphaFold) for nanobinder design:
1. Loads structure files
2. Selects target chain(s)
3. Extracts ectodomain (removes signal peptide, TM helix, cytosolic tail)
4. For AlphaFold: trims low-confidence termini (pLDDT < 70)
5. Outputs clean 3D model ready for analysis

Inputs:
- Target ID: gene symbol / UniProt (e.g. CD274 / Q9NZQ7)
- Structure file: PDB or AlphaFold model
- Chain ID and topology info (which side is extracellular, if TM protein)

Output:
- Clean 3D model (ectodomain only) ready for analysis
"""

from __future__ import annotations

import os
from typing import Optional, List, Tuple, Dict, Any
from dataclasses import dataclass
from pathlib import Path

try:
    from Bio import PDB
    from Bio.PDB import PDBIO, Select
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False
    PDB = None
    PDBIO = None
    Select = None

import numpy as np


@dataclass
class StructureInfo:
    """Information about a protein structure"""
    target_id: str  # Gene symbol or UniProt ID
    structure_file: str  # Path to PDB file
    structure_type: str  # "pdb" or "alphafold"
    chain_id: str  # Target chain ID (e.g., "A")
    domain_type: str  # "ectodomain", "full_length", etc.
    domain_residues: Optional[Tuple[int, int]] = None  # (start, end) if domain is extracted
    is_transmembrane: bool = False  # Whether this is a transmembrane protein
    extracellular_side: Optional[str] = None  # "N-terminal" or "C-terminal" for TM proteins
    notes: Optional[str] = None


@dataclass
class PreprocessedStructure:
    """Result of structure preprocessing"""
    structure: Any  # BioPython Structure object
    chain_id: str  # Selected chain ID
    residues: List[int]  # Residue numbers in the cleaned structure
    domain_start: int  # First residue number in cleaned structure
    domain_end: int  # Last residue number in cleaned structure
    trimmed_n_term: int  # Number of residues trimmed from N-terminus
    trimmed_c_term: int  # Number of residues trimmed from C-terminus
    low_confidence_trimmed: bool = False  # Whether low-confidence regions were trimmed
    output_file: Optional[str] = None  # Path to saved cleaned structure


class ChainSelector(Select):
    """BioPython Select class to select specific chain(s)"""
    
    def __init__(self, chain_ids: List[str]):
        self.chain_ids = set(chain_ids)
    
    def accept_chain(self, chain):
        return chain.id in self.chain_ids


class EctodomainSelector(Select):
    """BioPython Select class to select ectodomain residues"""
    
    def __init__(self, chain_id: str, start_res: int, end_res: int):
        self.chain_id = chain_id
        self.start_res = start_res
        self.end_res = end_res
    
    def accept_chain(self, chain):
        return chain.id == self.chain_id
    
    def accept_residue(self, residue):
        """Accept residues within the domain range"""
        res_id = residue.id
        # Handle insertion codes (e.g., (' ', 100, ' '))
        if isinstance(res_id, tuple) and len(res_id) >= 2:
            resnum = res_id[1]
            if isinstance(resnum, int):
                return self.start_res <= resnum <= self.end_res
        return False


class ConfidenceSelector(Select):
    """BioPython Select class to filter by confidence (pLDDT for AlphaFold)"""
    
    def __init__(self, chain_id: str, min_confidence: float = 70.0):
        self.chain_id = chain_id
        self.min_confidence = min_confidence
    
    def accept_chain(self, chain):
        return chain.id == self.chain_id
    
    def accept_residue(self, residue):
        """Accept residues with confidence >= min_confidence"""
        # Check B-factor of CA atom (pLDDT in AlphaFold)
        try:
            ca_atom = residue['CA']
            bfactor = ca_atom.get_bfactor()
            return bfactor >= self.min_confidence
        except (KeyError, AttributeError):
            # If no CA atom or no B-factor, accept by default
            return True


def load_structure(structure_file: str) -> Any:
    """
    Load a structure file (PDB or AlphaFold) using BioPython.
    
    Args:
        structure_file: Path to PDB file
        
    Returns:
        BioPython Structure object
        
    Raises:
        ImportError: If BioPython is not available
        FileNotFoundError: If structure file doesn't exist
        Exception: If structure parsing fails
    """
    if not BIOPYTHON_AVAILABLE:
        raise ImportError(
            "BioPython is required for structure preprocessing. "
            "Install with: pip install biopython"
        )
    
    if not os.path.exists(structure_file):
        raise FileNotFoundError(f"Structure file not found: {structure_file}")
    
    parser = PDB.PDBParser(QUIET=True)  # QUIET=True suppresses warnings
    
    try:
        structure = parser.get_structure("structure", structure_file)
        return structure
    except Exception as e:
        raise Exception(f"Failed to parse structure file {structure_file}: {str(e)}")


def get_chain_residue_range(structure: Any, chain_id: str) -> Tuple[int, int]:
    """
    Get the residue number range for a chain.
    
    Args:
        structure: BioPython Structure object
        chain_id: Chain ID to check
        
    Returns:
        (min_residue, max_residue) tuple
    """
    residues = []
    for chain in structure.get_chains():
        if chain.id == chain_id:
            for residue in chain:
                res_id = residue.id
                if isinstance(res_id, tuple) and len(res_id) >= 2:
                    resnum = res_id[1]
                    if isinstance(resnum, int):
                        residues.append(resnum)
    
    if not residues:
        raise ValueError(f"No residues found in chain {chain_id}")
    
    return (min(residues), max(residues))


def extract_ectodomain(
    structure: Any,
    chain_id: str,
    domain_residues: Optional[Tuple[int, int]] = None,
    signal_peptide_end: Optional[int] = None,
    tm_helix_start: Optional[int] = None,
    tm_helix_end: Optional[int] = None,
    cytosolic_tail_start: Optional[int] = None
) -> Tuple[Any, int, int, int, int]:
    """
    Extract ectodomain from structure by removing signal peptide, TM helix, and cytosolic tail.
    
    Args:
        structure: BioPython Structure object
        chain_id: Target chain ID
        domain_residues: Optional (start, end) tuple specifying domain directly
        signal_peptide_end: Last residue of signal peptide (if known)
        tm_helix_start: First residue of transmembrane helix (if known)
        tm_helix_end: Last residue of transmembrane helix (if known)
        cytosolic_tail_start: First residue of cytosolic tail (if known)
        
    Returns:
        (cleaned_structure, start_res, end_res, trimmed_n, trimmed_c)
        - cleaned_structure: New Structure object with only ectodomain
        - start_res: First residue number in cleaned structure
        - end_res: Last residue number in cleaned structure
        - trimmed_n: Number of residues trimmed from N-terminus
        - trimmed_c: Number of residues trimmed from C-terminus
    """
    if not BIOPYTHON_AVAILABLE:
        raise ImportError("BioPython is required")
    
    # Get full chain range
    full_start, full_end = get_chain_residue_range(structure, chain_id)
    
    # Determine domain boundaries
    if domain_residues:
        # Use provided domain residues directly
        start_res, end_res = domain_residues
        trimmed_n = start_res - full_start if start_res > full_start else 0
        trimmed_c = full_end - end_res if end_res < full_end else 0
    else:
        # Calculate from topology information
        start_res = full_start
        end_res = full_end
        
        # Remove signal peptide
        if signal_peptide_end:
            start_res = signal_peptide_end + 1
            trimmed_n = signal_peptide_end - full_start + 1
        else:
            trimmed_n = 0
        
        # Remove transmembrane helix
        if tm_helix_start and tm_helix_end:
            # For ectodomain, we want everything before TM helix
            end_res = tm_helix_start - 1
            trimmed_c = full_end - end_res
        elif cytosolic_tail_start:
            # Remove cytosolic tail
            end_res = cytosolic_tail_start - 1
            trimmed_c = full_end - end_res
        else:
            trimmed_c = 0
    
    # Create selector for ectodomain
    selector = EctodomainSelector(chain_id, start_res, end_res)
    
    # Create new structure with only selected residues
    io = PDBIO()
    io.set_structure(structure)
    
    # Create a new structure
    new_structure = PDB.StructureBuilder.Structure("ectodomain")
    new_model = PDB.Model.Model(0)
    new_structure.add(new_model)
    
    # Copy selected chain and residues
    for chain in structure.get_chains():
        if chain.id == chain_id:
            new_chain = PDB.Chain.Chain(chain_id)
            new_model.add(new_chain)
            
            for residue in chain:
                res_id = residue.id
                if isinstance(res_id, tuple) and len(res_id) >= 2:
                    resnum = res_id[1]
                    if isinstance(resnum, int) and start_res <= resnum <= end_res:
                        new_chain.add(residue.copy())
    
    return new_structure, start_res, end_res, trimmed_n, trimmed_c


def trim_low_confidence_alphafold(
    structure: Any,
    chain_id: str,
    min_plddt: float = 70.0
) -> Tuple[Any, bool, int, int]:
    """
    Trim low-confidence termini from AlphaFold structure.
    
    In AlphaFold PDB files, the B-factor column contains pLDDT scores.
    This function removes residues at N- and C-termini with pLDDT < min_plddt.
    
    Args:
        structure: BioPython Structure object (AlphaFold model)
        chain_id: Target chain ID
        min_plddt: Minimum pLDDT score to keep (default: 70.0)
        
    Returns:
        (cleaned_structure, was_trimmed, trimmed_n, trimmed_c)
        - cleaned_structure: New Structure object with low-confidence regions removed
        - was_trimmed: Whether any trimming occurred
        - trimmed_n: Number of residues trimmed from N-terminus
        - trimmed_c: Number of residues trimmed from C-terminus
    """
    if not BIOPYTHON_AVAILABLE:
        raise ImportError("BioPython is required")
    
    # Get chain and collect residues with their confidence scores
    chain = None
    for c in structure.get_chains():
        if c.id == chain_id:
            chain = c
            break
    
    if chain is None:
        raise ValueError(f"Chain {chain_id} not found in structure")
    
    # Collect residues with their confidence scores
    residue_confidences = []
    for residue in chain:
        try:
            ca_atom = residue['CA']
            bfactor = ca_atom.get_bfactor()
            res_id = residue.id
            if isinstance(res_id, tuple) and len(res_id) >= 2:
                resnum = res_id[1]
                if isinstance(resnum, int):
                    residue_confidences.append((resnum, bfactor, residue))
        except (KeyError, AttributeError):
            continue
    
    if not residue_confidences:
        # No confidence data available, return structure unchanged
        return structure, False, 0, 0
    
    # Sort by residue number
    residue_confidences.sort(key=lambda x: x[0])
    
    # Find first and last residues with confidence >= min_plddt
    start_idx = None
    end_idx = None
    
    for i, (resnum, conf, _) in enumerate(residue_confidences):
        if conf >= min_plddt:
            if start_idx is None:
                start_idx = i
            end_idx = i
    
    if start_idx is None:
        # No residues meet confidence threshold - return all (better than nothing)
        return structure, False, 0, 0
    
    # Calculate trimmed residues
    trimmed_n = start_idx
    trimmed_c = len(residue_confidences) - end_idx - 1
    
    # Get residue number range to keep
    start_res = residue_confidences[start_idx][0]
    end_res = residue_confidences[end_idx][0]
    
    # Create new structure with only high-confidence residues
    new_structure = PDB.StructureBuilder.Structure("trimmed")
    new_model = PDB.Model.Model(0)
    new_structure.add(new_model)
    
    new_chain = PDB.Chain.Chain(chain_id)
    new_model.add(new_chain)
    
    for residue in chain:
        res_id = residue.id
        if isinstance(res_id, tuple) and len(res_id) >= 2:
            resnum = res_id[1]
            if isinstance(resnum, int) and start_res <= resnum <= end_res:
                # Check confidence before adding
                try:
                    ca_atom = residue['CA']
                    bfactor = ca_atom.get_bfactor()
                    if bfactor >= min_plddt:
                        new_chain.add(residue.copy())
                except (KeyError, AttributeError):
                    # If no CA or B-factor, include it
                    new_chain.add(residue.copy())
    
    was_trimmed = trimmed_n > 0 or trimmed_c > 0
    return new_structure, was_trimmed, trimmed_n, trimmed_c


def select_chain(structure: Any, chain_id: str) -> Any:
    """
    Select a specific chain from structure.
    
    Args:
        structure: BioPython Structure object
        chain_id: Chain ID to select
        
    Returns:
        New Structure object with only the selected chain
    """
    if not BIOPYTHON_AVAILABLE:
        raise ImportError("BioPython is required")
    
    selector = ChainSelector([chain_id])
    
    # Create new structure with only selected chain
    new_structure = PDB.StructureBuilder.Structure("selected_chain")
    new_model = PDB.Model.Model(0)
    new_structure.add(new_model)
    
    for chain in structure.get_chains():
        if chain.id == chain_id:
            new_chain = PDB.Chain.Chain(chain_id)
            new_model.add(new_chain)
            for residue in chain:
                new_chain.add(residue.copy())
            break
    
    return new_structure


def save_structure(structure: Any, output_file: str) -> None:
    """
    Save a BioPython Structure object to a PDB file.
    
    Args:
        structure: BioPython Structure object
        output_file: Path to output PDB file
    """
    if not BIOPYTHON_AVAILABLE:
        raise ImportError("BioPython is required")
    
    # Ensure output directory exists
    os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)
    
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_file)


def preprocess_structure(
    structure_info: StructureInfo,
    output_file: Optional[str] = None,
    trim_alphafold_confidence: bool = True,
    min_plddt: float = 70.0
) -> PreprocessedStructure:
    """
    Main function: Pre-process a structure for nanobinder design.
    
    Steps:
    1. Load structure (PDB or AlphaFold)
    2. Select target chain(s)
    3. Extract ectodomain (if TM protein) - remove signal peptide, TM helix, cytosolic tail
    4. For AlphaFold: trim low-confidence termini (pLDDT < min_plddt)
    5. Save cleaned structure
    
    Args:
        structure_info: StructureInfo object with structure details
        output_file: Optional path to save cleaned structure (if None, auto-generates)
        trim_alphafold_confidence: Whether to trim low-confidence regions for AlphaFold (default: True)
        min_plddt: Minimum pLDDT score for AlphaFold trimming (default: 70.0)
        
    Returns:
        PreprocessedStructure object with cleaned structure and metadata
    """
    if not BIOPYTHON_AVAILABLE:
        raise ImportError(
            "BioPython is required for structure preprocessing. "
            "Install with: pip install biopython"
        )
    
    # Step 1: Load structure
    structure = load_structure(structure_info.structure_file)
    
    # Step 2: Select target chain
    structure = select_chain(structure, structure_info.chain_id)
    
    # Step 3: Extract ectodomain if needed
    trimmed_n_term = 0
    trimmed_c_term = 0
    domain_start = None
    domain_end = None
    
    if structure_info.domain_type == "ectodomain" and structure_info.is_transmembrane:
        # Extract ectodomain using domain_residues or topology info
        structure, domain_start, domain_end, trimmed_n_term, trimmed_c_term = extract_ectodomain(
            structure,
            structure_info.chain_id,
            domain_residues=structure_info.domain_residues
        )
    elif structure_info.domain_residues:
        # Extract specific domain range
        structure, domain_start, domain_end, trimmed_n_term, trimmed_c_term = extract_ectodomain(
            structure,
            structure_info.chain_id,
            domain_residues=structure_info.domain_residues
        )
    else:
        # Get full chain range
        domain_start, domain_end = get_chain_residue_range(structure, structure_info.chain_id)
    
    # Step 4: Trim low-confidence regions for AlphaFold
    low_confidence_trimmed = False
    if structure_info.structure_type == "alphafold" and trim_alphafold_confidence:
        structure, was_trimmed, trim_n, trim_c = trim_low_confidence_alphafold(
            structure,
            structure_info.chain_id,
            min_plddt=min_plddt
        )
        if was_trimmed:
            low_confidence_trimmed = True
            trimmed_n_term += trim_n
            trimmed_c_term += trim_c
            # Update domain range after trimming
            domain_start, domain_end = get_chain_residue_range(structure, structure_info.chain_id)
    
    # Step 5: Get final residue list
    residues = []
    for chain in structure.get_chains():
        if chain.id == structure_info.chain_id:
            for residue in chain:
                res_id = residue.id
                if isinstance(res_id, tuple) and len(res_id) >= 2:
                    resnum = res_id[1]
                    if isinstance(resnum, int):
                        residues.append(resnum)
            break
    
    # Step 6: Save cleaned structure
    if output_file is None:
        # Auto-generate output filename
        base_name = os.path.splitext(os.path.basename(structure_info.structure_file))[0]
        output_dir = os.path.dirname(os.path.abspath(structure_info.structure_file))
        output_file = os.path.join(output_dir, f"{base_name}_cleaned.pdb")
    
    save_structure(structure, output_file)
    
    return PreprocessedStructure(
        structure=structure,
        chain_id=structure_info.chain_id,
        residues=sorted(residues),
        domain_start=domain_start,
        domain_end=domain_end,
        trimmed_n_term=trimmed_n_term,
        trimmed_c_term=trimmed_c_term,
        low_confidence_trimmed=low_confidence_trimmed,
        output_file=output_file
    )



