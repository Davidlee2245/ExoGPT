"""
Per-Residue Feature Computation
================================

Computes per-residue features for structure analysis:
1. SASA (Solvent Accessible Surface Area) - surface exposure
2. Secondary structure - helix/sheet/loop classification
3. Flexibility/confidence - B-factor (PDB) or pLDDT (AlphaFold)

For each residue in the cleaned model:
- SASA: Run FreeSASA/NACCESS or use approximation
- Secondary structure: Run DSSP → classify helix/sheet/loop
- Flexibility: Use B-factor (PDB) or pLDDT (AlphaFold)
- Flag surface residues (SASA > threshold, e.g. 30 Å²)
- Mark low-confidence/flexible residues (avoid as cores of epitopes)
"""

from __future__ import annotations

import os
import subprocess
import tempfile
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass
import numpy as np

try:
    from Bio import PDB
    from Bio.PDB import DSSP
    from Bio.PDB.DSSP import dssp_dict_from_pdb_file
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False
    PDB = None
    DSSP = None

# Try to import FreeSASA (optional)
try:
    import freesasa
    FREESASA_AVAILABLE = True
except ImportError:
    FREESASA_AVAILABLE = False
    freesasa = None


@dataclass
class ResidueFeatures:
    """Per-residue features for structure analysis"""
    residue_number: int
    residue_name: str  # 3-letter code (e.g., "ALA")
    chain_id: str
    
    # SASA (Solvent Accessible Surface Area)
    sasa: float  # Å²
    is_surface: bool  # True if SASA > threshold
    
    # Secondary structure
    secondary_structure: str  # H (helix), E (sheet), C (coil/loop), etc.
    secondary_structure_category: str  # "helix", "sheet", "loop"
    
    # Flexibility/confidence
    bfactor: float  # B-factor (PDB) or pLDDT (AlphaFold)
    is_low_confidence: bool  # True if confidence is too low
    is_highly_flexible: bool  # True if B-factor is very high
    
    # Additional info
    phi: Optional[float] = None  # Phi angle (degrees)
    psi: Optional[float] = None  # Psi angle (degrees)


@dataclass
class StructureFeatures:
    """Complete feature set for a structure"""
    structure_file: str
    chain_id: str
    residues: List[ResidueFeatures]
    sasa_method: str  # "freesasa", "naccess", "approximation", "unavailable"
    dssp_available: bool
    total_residues: int
    surface_residues: int
    low_confidence_residues: int
    highly_flexible_residues: int


def check_dssp_available() -> bool:
    """Check if DSSP executable is available"""
    try:
        result = subprocess.run(
            ["dssp", "--version"],
            capture_output=True,
            timeout=5
        )
        return result.returncode == 0
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False


def check_freesasa_available() -> bool:
    """Check if FreeSASA Python library is available"""
    return FREESASA_AVAILABLE


def check_naccess_available() -> bool:
    """Check if NACCESS executable is available"""
    try:
        result = subprocess.run(
            ["naccess", "-h"],
            capture_output=True,
            timeout=5
        )
        return result.returncode == 0
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False


def calculate_sasa_freesasa(
    structure_file: str,
    chain_id: str
) -> Optional[Dict[int, float]]:
    """
    Calculate SASA using FreeSASA Python library.
    
    Args:
        structure_file: Path to PDB file
        chain_id: Chain ID to analyze
        
    Returns:
        Dictionary mapping residue number to SASA (Å²), or None if failed
    """
    if not FREESASA_AVAILABLE:
        return None
    
    try:
        # Load structure with FreeSASA
        structure = freesasa.Structure(structure_file)
        result = freesasa.calc(structure)
        
        # Extract per-residue SASA
        sasa_dict = {}
        for i, atom in enumerate(structure.atoms()):
            if atom.chainId() == chain_id:
                resnum = atom.residueNumber()
                sasa_value = result.atomArea(i)
                
                # Sum SASA for all atoms in the same residue
                if resnum not in sasa_dict:
                    sasa_dict[resnum] = 0.0
                sasa_dict[resnum] += sasa_value
        
        return sasa_dict
    except Exception as e:
        print(f"FreeSASA calculation failed: {e}")
        return None


def calculate_sasa_naccess(
    structure_file: str,
    chain_id: str,
    temp_dir: Optional[str] = None
) -> Optional[Dict[int, float]]:
    """
    Calculate SASA using NACCESS executable.
    
    Args:
        structure_file: Path to PDB file
        chain_id: Chain ID to analyze
        temp_dir: Temporary directory for NACCESS output
        
    Returns:
        Dictionary mapping residue number to SASA (Å²), or None if failed
    """
    if not check_naccess_available():
        return None
    
    try:
        if temp_dir is None:
            temp_dir = tempfile.mkdtemp()
        
        # Run NACCESS
        base_name = os.path.splitext(os.path.basename(structure_file))[0]
        output_base = os.path.join(temp_dir, base_name)
        
        result = subprocess.run(
            ["naccess", structure_file, "-h"],
            capture_output=True,
            timeout=60,
            cwd=temp_dir
        )
        
        if result.returncode != 0:
            return None
        
        # Parse NACCESS output (.rsa file)
        rsa_file = f"{output_base}.rsa"
        if not os.path.exists(rsa_file):
            return None
        
        sasa_dict = {}
        with open(rsa_file, 'r') as f:
            for line in f:
                if line.startswith("RES"):
                    # Parse NACCESS format
                    # RES NUM      All-atoms   Total-Side   Main-Chain    Non-polar    All polar
                    parts = line.split()
                    if len(parts) >= 3:
                        try:
                            resnum = int(parts[1])
                            chain = parts[2] if len(parts) > 2 else ""
                            if chain == chain_id:
                                sasa = float(parts[3])  # All-atoms SASA
                                sasa_dict[resnum] = sasa
                        except (ValueError, IndexError):
                            continue
        
        return sasa_dict
    except Exception as e:
        print(f"NACCESS calculation failed: {e}")
        return None


def calculate_sasa_approximation(
    structure: Any,
    chain_id: str
) -> Dict[int, float]:
    """
    Approximate SASA using a simple method based on atom exposure.
    
    This is a fallback when FreeSASA/NACCESS are not available.
    Uses a simple sphere-based approximation.
    
    Args:
        structure: BioPython Structure object
        chain_id: Chain ID to analyze
        
    Returns:
        Dictionary mapping residue number to approximate SASA (Å²)
    """
    if not BIOPYTHON_AVAILABLE:
        return {}
    
    # Get all CA atoms and their coordinates
    ca_atoms = []
    residue_map = {}
    
    for chain in structure.get_chains():
        if chain.id == chain_id:
            for residue in chain:
                res_id = residue.id
                if isinstance(res_id, tuple) and len(res_id) >= 2:
                    resnum = res_id[1]
                    if isinstance(resnum, int):
                        try:
                            ca = residue['CA']
                            ca_atoms.append(ca.get_coord())
                            residue_map[len(ca_atoms) - 1] = resnum
                        except KeyError:
                            continue
    
    if not ca_atoms:
        return {}
    
    ca_coords = np.array(ca_atoms)
    
    # Simple approximation: calculate exposure based on local neighbor density
    sasa_dict = {}
    cutoff = 9.0  # Å
    neighbor_capacity = 18.0  # approximate number of neighbors for buried residues
    
    for i, coord in enumerate(ca_coords):
        resnum = residue_map[i]
        distances = np.linalg.norm(ca_coords - coord, axis=1)
        neighbor_count = np.sum((distances <= cutoff) & (distances > 0.1))
        
        exposure = max(0.0, 1.0 - (neighbor_count / neighbor_capacity))
        approximate_sasa = 80.0 * (exposure ** 1.5)
        
        sasa_dict[resnum] = float(approximate_sasa)
    
    return sasa_dict


def calculate_sasa(
    structure_file: str,
    structure: Any,
    chain_id: str,
    method_preference: List[str] = None
) -> Tuple[Dict[int, float], str]:
    """
    Calculate SASA using the best available method.
    
    Args:
        structure_file: Path to PDB file
        structure: BioPython Structure object
        chain_id: Chain ID to analyze
        method_preference: Preferred methods in order (default: ["freesasa", "naccess", "approximation"])
        
    Returns:
        (sasa_dict, method_used) tuple
    """
    if method_preference is None:
        method_preference = ["freesasa", "naccess", "approximation"]
    
    sasa_dict = None
    method_used = "unavailable"
    
    for method in method_preference:
        if method == "freesasa" and check_freesasa_available():
            sasa_dict = calculate_sasa_freesasa(structure_file, chain_id)
            if sasa_dict:
                method_used = "freesasa"
                break
        elif method == "naccess" and check_naccess_available():
            sasa_dict = calculate_sasa_naccess(structure_file, chain_id)
            if sasa_dict:
                method_used = "naccess"
                break
        elif method == "approximation":
            sasa_dict = calculate_sasa_approximation(structure, chain_id)
            if sasa_dict:
                method_used = "approximation"
                break
    
    if sasa_dict is None:
        sasa_dict = {}
        method_used = "unavailable"
    
    return sasa_dict, method_used


def calculate_secondary_structure(
    structure_file: str,
    structure: Any,
    chain_id: str
) -> Tuple[Dict[int, Dict[str, Any]], bool]:
    """
    Calculate secondary structure using DSSP.
    
    Args:
        structure_file: Path to PDB file
        structure: BioPython Structure object
        chain_id: Chain ID to analyze
        
    Returns:
        (dssp_dict, dssp_available) tuple
        - dssp_dict: Dictionary mapping (chain_id, res_id) to DSSP features
        - dssp_available: Whether DSSP was successfully run
    """
    if not BIOPYTHON_AVAILABLE:
        return {}, False
    
    dssp_available = check_dssp_available()
    
    if not dssp_available:
        # Try to use BioPython's built-in DSSP (requires DSSP executable)
        return {}, False
    
    try:
        # Use BioPython's DSSP
        dssp_dict = dssp_dict_from_pdb_file(structure_file)
        
        # Filter for our chain
        filtered_dict = {}
        for key, value in dssp_dict.items():
            if key[0] == chain_id:
                filtered_dict[key] = value
        
        return filtered_dict, True
    except Exception as e:
        print(f"DSSP calculation failed: {e}")
        # Fallback: use simple secondary structure from structure
        return {}, False


def get_secondary_structure_simple(
    structure: Any,
    chain_id: str
) -> Dict[int, str]:
    """
    Simple secondary structure assignment based on phi/psi angles.
    This is a fallback when DSSP is not available.
    
    Args:
        structure: BioPython Structure object
        chain_id: Chain ID to analyze
        
    Returns:
        Dictionary mapping residue number to secondary structure code
    """
    if not BIOPYTHON_AVAILABLE:
        return {}
    
    ss_dict = {}
    
    # Get CA coordinates for simple helix/sheet detection
    ca_coords = []
    residue_numbers = []
    
    for chain in structure.get_chains():
        if chain.id == chain_id:
            for residue in chain:
                res_id = residue.id
                if isinstance(res_id, tuple) and len(res_id) >= 2:
                    resnum = res_id[1]
                    if isinstance(resnum, int):
                        try:
                            ca = residue['CA']
                            ca_coords.append(ca.get_coord())
                            residue_numbers.append(resnum)
                        except KeyError:
                            continue
    
    if len(ca_coords) < 3:
        return {}
    
    ca_coords = np.array(ca_coords)
    
    # Simple classification based on CA-CA distances and angles
    for i, resnum in enumerate(residue_numbers):
        ss_dict[resnum] = "C"  # Default to coil
        
        if i < 2 or i >= len(ca_coords) - 2:
            continue
        
        # Calculate CA-CA distances
        d1 = np.linalg.norm(ca_coords[i] - ca_coords[i-1])
        d2 = np.linalg.norm(ca_coords[i+1] - ca_coords[i])
        
        # Helix: CA-CA distance ~5.5 Å
        if 4.5 < d1 < 6.5 and 4.5 < d2 < 6.5:
            ss_dict[resnum] = "H"
        # Sheet: CA-CA distance ~3.8 Å (in strand)
        elif 3.0 < d1 < 4.5 and 3.0 < d2 < 4.5:
            ss_dict[resnum] = "E"
    
    return ss_dict


def extract_bfactor(
    structure: Any,
    chain_id: str
) -> Dict[int, float]:
    """
    Extract B-factor (PDB) or pLDDT (AlphaFold) for each residue.
    
    Args:
        structure: BioPython Structure object
        chain_id: Chain ID to analyze
        
    Returns:
        Dictionary mapping residue number to B-factor/pLDDT
    """
    if not BIOPYTHON_AVAILABLE:
        return {}
    
    bfactor_dict = {}
    
    for chain in structure.get_chains():
        if chain.id == chain_id:
            for residue in chain:
                res_id = residue.id
                if isinstance(res_id, tuple) and len(res_id) >= 2:
                    resnum = res_id[1]
                    if isinstance(resnum, int):
                        try:
                            # Use CA atom B-factor as representative
                            ca_atom = residue['CA']
                            bfactor = ca_atom.get_bfactor()
                            bfactor_dict[resnum] = bfactor
                        except KeyError:
                            # If no CA, try to get average of all atoms
                            bfactors = []
                            for atom in residue:
                                bfactors.append(atom.get_bfactor())
                            if bfactors:
                                bfactor_dict[resnum] = np.mean(bfactors)
    
    return bfactor_dict


def ss_code_to_category(ss_code: str) -> str:
    """Convert DSSP secondary structure code to category"""
    # DSSP codes: H=helix, E=sheet, C=coil/loop
    if ss_code in ['H', 'G', 'I']:  # Helix types
        return "helix"
    elif ss_code in ['E', 'B']:  # Sheet types
        return "sheet"
    else:  # Coil, turn, etc.
        return "loop"


def compute_residue_features(
    structure_file: str,
    structure: Any,
    chain_id: str,
    sasa_threshold: float = 30.0,
    low_confidence_threshold: float = 50.0,
    high_flexibility_threshold: float = 80.0,
    structure_type: str = "pdb"
) -> StructureFeatures:
    """
    Main function: Compute per-residue features for a structure.
    
    Args:
        structure_file: Path to PDB file
        structure: BioPython Structure object
        chain_id: Chain ID to analyze
        sasa_threshold: SASA threshold for surface residues (Å², default: 30.0)
        low_confidence_threshold: B-factor/pLDDT threshold for low confidence (default: 50.0)
        high_flexibility_threshold: B-factor threshold for high flexibility (default: 80.0)
        structure_type: "pdb" or "alphafold"
        
    Returns:
        StructureFeatures object with all per-residue features
    """
    if not BIOPYTHON_AVAILABLE:
        raise ImportError("BioPython is required for feature computation")
    
    # 1. Calculate SASA
    sasa_dict, sasa_method = calculate_sasa(structure_file, structure, chain_id)
    
    # 2. Calculate secondary structure
    dssp_dict, dssp_available = calculate_secondary_structure(structure_file, structure, chain_id)
    
    # Fallback to simple method if DSSP not available
    if not dssp_available:
        ss_dict = get_secondary_structure_simple(structure, chain_id)
    else:
        # Convert DSSP dict to simple residue number -> SS code
        ss_dict = {}
        for (chain, res_id), dssp_data in dssp_dict.items():
            if chain == chain_id:
                if isinstance(res_id, tuple) and len(res_id) >= 2:
                    resnum = res_id[1]
                    if isinstance(resnum, int):
                        ss_code = dssp_data[1]  # Secondary structure code
                        ss_dict[resnum] = ss_code
    
    # 3. Extract B-factor/pLDDT
    bfactor_dict = extract_bfactor(structure, chain_id)
    
    # 4. Compile per-residue features
    residue_features_list = []
    
    for chain in structure.get_chains():
        if chain.id == chain_id:
            for residue in chain:
                res_id = residue.id
                if isinstance(res_id, tuple) and len(res_id) >= 2:
                    resnum = res_id[1]
                    if isinstance(resnum, int):
                        # Get residue name
                        resname = residue.get_resname()
                        
                        # Get SASA
                        sasa = sasa_dict.get(resnum, 0.0)
                        is_surface = sasa > sasa_threshold
                        
                        # Get secondary structure
                        ss_code = ss_dict.get(resnum, "C")
                        ss_category = ss_code_to_category(ss_code)
                        
                        # Get B-factor/pLDDT
                        bfactor = bfactor_dict.get(resnum, 0.0)
                        
                        # For AlphaFold, pLDDT is in B-factor column
                        # Low confidence = low pLDDT
                        # For PDB, high B-factor = high flexibility
                        if structure_type == "alphafold":
                            is_low_confidence = bfactor < low_confidence_threshold
                            is_highly_flexible = False  # pLDDT doesn't indicate flexibility
                        else:
                            is_low_confidence = bfactor > high_flexibility_threshold  # High B-factor = low confidence
                            is_highly_flexible = bfactor > high_flexibility_threshold
                        
                        # Get phi/psi if available from DSSP
                        phi = None
                        psi = None
                        if dssp_available:
                            dssp_key = (chain_id, res_id)
                            if dssp_key in dssp_dict:
                                dssp_data = dssp_dict[dssp_key]
                                phi = dssp_data[4] if len(dssp_data) > 4 else None
                                psi = dssp_data[5] if len(dssp_data) > 5 else None
                        
                        residue_features = ResidueFeatures(
                            residue_number=resnum,
                            residue_name=resname,
                            chain_id=chain_id,
                            sasa=sasa,
                            is_surface=is_surface,
                            secondary_structure=ss_code,
                            secondary_structure_category=ss_category,
                            bfactor=bfactor,
                            is_low_confidence=is_low_confidence,
                            is_highly_flexible=is_highly_flexible,
                            phi=phi,
                            psi=psi
                        )
                        
                        residue_features_list.append(residue_features)
    
    # Calculate summary statistics
    total_residues = len(residue_features_list)
    surface_residues = sum(1 for rf in residue_features_list if rf.is_surface)
    low_confidence_residues = sum(1 for rf in residue_features_list if rf.is_low_confidence)
    highly_flexible_residues = sum(1 for rf in residue_features_list if rf.is_highly_flexible)
    
    return StructureFeatures(
        structure_file=structure_file,
        chain_id=chain_id,
        residues=residue_features_list,
        sasa_method=sasa_method,
        dssp_available=dssp_available,
        total_residues=total_residues,
        surface_residues=surface_residues,
        low_confidence_residues=low_confidence_residues,
        highly_flexible_residues=highly_flexible_residues
    )



