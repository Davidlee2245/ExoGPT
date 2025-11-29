#!/usr/bin/env python
"""
Helper script to download PDB structure files referenced in structure_mapping.csv
"""
import os
import sys
import pandas as pd
import requests
from pathlib import Path

def download_pdb(pdb_id: str, output_path: str) -> bool:
    """Download a PDB file from RCSB PDB."""
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        with open(output_path, 'w') as f:
            f.write(response.text)
        print(f"✓ Downloaded {pdb_id} to {output_path}")
        return True
    except Exception as e:
        print(f"✗ Failed to download {pdb_id}: {e}")
        return False

def download_alphafold(uniprot_id: str, output_path: str) -> bool:
    """Download an AlphaFold model from AlphaFold DB."""
    # AlphaFold DB URL format: AF-{uniprot}-F1-model_v{version}.pdb
    # Try different versions (v4, v3, v2)
    for version in ["v4", "v3", "v2"]:
        url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_{version}.pdb"
        try:
            response = requests.get(url, timeout=30)
            if response.status_code == 200:
                os.makedirs(os.path.dirname(output_path), exist_ok=True)
                with open(output_path, 'w') as f:
                    f.write(response.text)
                print(f"✓ Downloaded AlphaFold model {version} for {uniprot_id} to {output_path}")
                return True
        except Exception:
            continue
    
    # If all versions fail, try alternative URL format
    alt_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    try:
        response = requests.get(alt_url, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if data and len(data) > 0:
                pdb_url = data[0].get("cifUrl", "").replace(".cif", ".pdb")
                if pdb_url:
                    pdb_response = requests.get(pdb_url, timeout=30)
                    if pdb_response.status_code == 200:
                        os.makedirs(os.path.dirname(output_path), exist_ok=True)
                        with open(output_path, 'w') as f:
                            f.write(pdb_response.text)
                        print(f"✓ Downloaded AlphaFold model for {uniprot_id} to {output_path}")
                        return True
    except Exception:
        pass
    
    print(f"✗ Failed to download AlphaFold model for {uniprot_id}")
    print(f"  Try downloading manually from: https://alphafold.ebi.ac.uk/entry/{uniprot_id}")
    return False

def main():
    project_root = Path(__file__).parent
    mapping_path = project_root / "data" / "structure_mapping.csv"
    
    if not mapping_path.exists():
        print(f"Error: structure_mapping.csv not found at {mapping_path}")
        sys.exit(1)
    
    # Load structure mapping
    df = pd.read_csv(mapping_path, quotechar='"', skipinitialspace=True)
    
    # Create structures directory
    structures_dir = project_root / "structures"
    structures_dir.mkdir(exist_ok=True)
    
    print(f"Downloading structure files to {structures_dir}...")
    print("=" * 60)
    
    downloaded = 0
    failed = 0
    
    for _, row in df.iterrows():
        model_path = str(row.get("model_path", ""))
        structure_type = str(row.get("structure_type", "")).lower()
        pdb_id = str(row.get("pdb_id", "")) if pd.notna(row.get("pdb_id")) else ""
        uniprot = str(row.get("uniprot", "")) if pd.notna(row.get("uniprot")) else ""
        
        if not model_path:
            continue
        
        # Make path absolute
        if not os.path.isabs(model_path):
            model_path = os.path.join(str(project_root), model_path.lstrip('./'))
        model_path = os.path.normpath(model_path)
        
        # Skip if file already exists
        if os.path.exists(model_path):
            print(f"⊘ {os.path.basename(model_path)} already exists, skipping")
            continue
        
        # Download based on structure type
        if structure_type == "pdb" and pdb_id:
            if download_pdb(pdb_id, model_path):
                downloaded += 1
            else:
                failed += 1
        elif structure_type == "alphafold" and uniprot:
            if download_alphafold(uniprot, model_path):
                downloaded += 1
            else:
                failed += 1
        else:
            print(f"⚠ Skipping {model_path}: missing pdb_id or uniprot")
            failed += 1
    
    print("=" * 60)
    print(f"Download complete: {downloaded} files downloaded, {failed} failed")
    
    if failed > 0:
        print("\nNote: Some files failed to download. You may need to:")
        print("1. Check internet connection")
        print("2. Verify PDB IDs and UniProt IDs in structure_mapping.csv")
        print("3. Download files manually from:")
        print("   - RCSB PDB: https://www.rcsb.org/")
        print("   - AlphaFold DB: https://alphafold.ebi.ac.uk/")

if __name__ == "__main__":
    main()

