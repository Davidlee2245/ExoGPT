"""
Step 3 — Nanobinder Design Orchestrator
-----------------------------------------

Orchestrates the nanobinder design pipeline:
1. RFdiffusion (backbone generation)
2. ProteinMPNN (sequence design)
3. AlphaFold-Multimer / ColabFold (complex prediction)

Generates reusable scripts and configs for execution on external machines/HPC.

Input: Step 2 JSON output (curated targets with epitopes)
Output: Design plan, config files, and execution scripts
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import List, Dict, Any, Optional, Callable

import yaml


@dataclass
class DesignConfig:
    """Configuration for a single design run (one target + one epitope)."""
    target_name: str
    gene_symbol: str
    uniprot: Optional[str]
    structure_type: str  # "pdb" or "alphafold"
    model_path: str  # Path to target PDB file
    pdb_id: Optional[str]
    target_chain_id: str  # Chain ID in PDB (default: "A")
    epitope_id: str
    epitope_residues: List[int]
    binder_length: int  # Binder length in amino acids (default: 80-100)
    binder_style: str  # e.g., "helical_mini_domain", "beta_sheet", "mixed"
    num_backbones: int  # Number of RFdiffusion designs (default: 10)
    num_sequences_per_backbone: int  # Number of ProteinMPNN sequences per backbone (default: 8)
    num_af_models: int  # Number of AlphaFold models per sequence (default: 5)
    num_af_recycles: int  # Number of AlphaFold recycles (default: 3)
    output_base_dir: str  # Base directory for all outputs
    rfdiffusion_path: Optional[str] = None  # Path to RFdiffusion executable/script
    proteinmpnn_path: Optional[str] = None  # Path to ProteinMPNN executable/script
    colabfold_path: Optional[str] = None  # Path to ColabFold executable/script
    alphafold_path: Optional[str] = None  # Path to AlphaFold-Multimer executable/script
    use_colabfold: bool = True  # Use ColabFold (True) or AlphaFold-Multimer (False)
    omitted_aa: List[str] = None  # Amino acids to omit in ProteinMPNN (default: ["C", "X"])
    binder_chain_id: str = "B"  # Chain ID for binder in complex (default: "B")
    
    def __post_init__(self):
        if self.omitted_aa is None:
            self.omitted_aa = ["C", "X"]


def load_step2_json(step2_json_path: str) -> Dict[str, Any]:
    """Load Step 2 JSON output."""
    with open(step2_json_path, "r") as f:
        return json.load(f)


def create_design_configs(
    step2_data: Dict[str, Any],
    output_base_dir: str = "./workflows",
    binder_length: int = 90,
    num_backbones: int = 10,
    num_sequences_per_backbone: int = 8,
    num_af_models: int = 5,
    num_af_recycles: int = 3,
    use_colabfold: bool = True,
    rfdiffusion_path: Optional[str] = None,
    proteinmpnn_path: Optional[str] = None,
    colabfold_path: Optional[str] = None,
    alphafold_path: Optional[str] = None,
) -> List[DesignConfig]:
    """
    Create design configurations for each target/epitope combination.
    
    Args:
        step2_data: Step 2 JSON output
        output_base_dir: Base directory for all outputs
        binder_length: Length of binder in amino acids
        num_backbones: Number of RFdiffusion backbone designs
        num_sequences_per_backbone: Number of ProteinMPNN sequences per backbone
        num_af_models: Number of AlphaFold models per sequence
        num_af_recycles: Number of AlphaFold recycles
        use_colabfold: Use ColabFold (True) or AlphaFold-Multimer (False)
        rfdiffusion_path: Path to RFdiffusion executable
        proteinmpnn_path: Path to ProteinMPNN executable
        colabfold_path: Path to ColabFold executable
        alphafold_path: Path to AlphaFold-Multimer executable
    
    Returns:
        List of DesignConfig objects
    """
    configs = []
    curated_targets = step2_data.get("curated_targets", [])
    
    for target in curated_targets:
        target_name = target.get("target_name", "")
        gene_symbol = target.get("gene_symbol", "")
        uniprot = target.get("uniprot")
        structure_type = target.get("structure_type", "pdb")
        model_path = target.get("model_path", "")
        pdb_id = target.get("pdb_id")
        epitopes = target.get("epitopes", [])
        
        # Determine target chain ID (default: "A")
        # In practice, this could be extracted from PDB or specified in config
        target_chain_id = "A"
        
        for epitope in epitopes:
            epitope_id = epitope.get("epitope_id", "")
            epitope_residues = epitope.get("epitope_residues", [])
            
            # Create unique output directory for this target/epitope
            safe_target = target_name.replace(" ", "_").replace("/", "_")
            safe_epitope = epitope_id.replace(" ", "_").replace("/", "_")
            design_output_dir = os.path.join(
                output_base_dir,
                "designs",
                safe_target,
                safe_epitope
            )
            
            config = DesignConfig(
                target_name=target_name,
                gene_symbol=gene_symbol,
                uniprot=uniprot,
                structure_type=structure_type,
                model_path=model_path,
                pdb_id=pdb_id,
                target_chain_id=target_chain_id,
                epitope_id=epitope_id,
                epitope_residues=epitope_residues,
                binder_length=binder_length,
                binder_style="helical_mini_domain",  # Default style
                num_backbones=num_backbones,
                num_sequences_per_backbone=num_sequences_per_backbone,
                num_af_models=num_af_models,
                num_af_recycles=num_af_recycles,
                output_base_dir=output_base_dir,
                rfdiffusion_path=rfdiffusion_path,
                proteinmpnn_path=proteinmpnn_path,
                colabfold_path=colabfold_path,
                alphafold_path=alphafold_path,
                use_colabfold=use_colabfold,
            )
            
            configs.append(config)
    
    return configs


def generate_rfdiffusion_script(config: DesignConfig) -> str:
    """
    Generate bash script to run RFdiffusion for backbone generation.
    
    This script:
    - Takes target PDB and epitope residues as input
    - Generates binder backbones using RFdiffusion
    - Outputs complex PDBs (target + binder backbone)
    """
    script_lines = [
        "#!/bin/bash",
        "# RFdiffusion Backbone Generation Script",
        f"# Target: {config.target_name} ({config.gene_symbol})",
        f"# Epitope: {config.epitope_id}",
        f"# Residues: {', '.join(map(str, config.epitope_residues))}",
        "",
        "set -e  # Exit on error",
        "",
        "# Configuration",
        f'TARGET_PDB="{config.model_path}"',
        f'TARGET_CHAIN="{config.target_chain_id}"',
        f'EPITOPE_RESIDUES="{",".join(map(str, config.epitope_residues))}"',
        f'BINDER_LENGTH={config.binder_length}',
        f'NUM_DESIGNS={config.num_backbones}',
        f'OUTPUT_DIR="{os.path.join(config.output_base_dir, "designs", config.target_name.replace(" ", "_"), config.epitope_id.replace(" ", "_"), "rfdiffusion_output")}"',
        "",
        "# Create output directory",
        "mkdir -p \"$OUTPUT_DIR\"",
        "",
        "# RFdiffusion command",
        "# NOTE: Adjust the command below based on your RFdiffusion installation",
        "# Example for RFdiffusion v1.0+:",
        f'if [ -z "$RFDIFFUSION_PATH" ]; then',
        '    RFDIFFUSION_CMD="python scripts/run_inference.py"',
        "else",
        '    RFDIFFUSION_CMD="$RFDIFFUSION_PATH"',
        "fi",
        "",
        "# Run RFdiffusion",
        "# Format: RFdiffusion with hotspot-guided design",
        "# The exact command depends on RFdiffusion version and API",
        "# Example command structure (adjust as needed):",
        "",
        "for i in $(seq 1 $NUM_DESIGNS); do",
        '    OUTPUT_PDB="$OUTPUT_DIR/complex_${i}.pdb"',
        "",
        "    # RFdiffusion inference command",
        "    # This is a template - adjust based on actual RFdiffusion API",
        '    $RFDIFFUSION_CMD \\',
        f'        --target_pdb "$TARGET_PDB" \\',
        f'        --target_chain "$TARGET_CHAIN" \\',
        f'        --hotspot_residues "$EPITOPE_RESIDUES" \\',
        f'        --binder_length $BINDER_LENGTH \\',
        f'        --output_pdb "$OUTPUT_PDB" \\',
        '        --num_samples 1',
        "",
        "    if [ $? -eq 0 ]; then",
        '        echo "✓ Generated backbone $i: $OUTPUT_PDB"',
        "    else",
        '        echo "✗ Failed to generate backbone $i"',
        "        exit 1",
        "    fi",
        "done",
        "",
        "echo \"",
        f"RFdiffusion completed: {config.num_backbones} backbones generated",
        f"Output directory: $OUTPUT_DIR",
        "\"",
    ]
    
    return "\n".join(script_lines)


def generate_proteinmpnn_script(config: DesignConfig) -> str:
    """
    Generate bash script to run ProteinMPNN for sequence design.
    
    This script:
    - Takes RFdiffusion output PDBs (complex with binder backbone)
    - Designs sequences for the binder chain
    - Outputs designed sequences and updated PDBs
    """
    rfdiffusion_output = os.path.join(
        config.output_base_dir,
        "designs",
        config.target_name.replace(" ", "_"),
        config.epitope_id.replace(" ", "_"),
        "rfdiffusion_output"
    )
    proteinmpnn_output = os.path.join(
        config.output_base_dir,
        "designs",
        config.target_name.replace(" ", "_"),
        config.epitope_id.replace(" ", "_"),
        "proteinmpnn_output"
    )
    
    omitted_aa_str = ",".join(config.omitted_aa)
    
    script_lines = [
        "#!/bin/bash",
        "# ProteinMPNN Sequence Design Script",
        f"# Target: {config.target_name} ({config.gene_symbol})",
        f"# Epitope: {config.epitope_id}",
        "",
        "set -e  # Exit on error",
        "",
        "# Configuration",
        f'RFDIFFUSION_OUTPUT="{rfdiffusion_output}"',
        f'PROTEINMPNN_OUTPUT="{proteinmpnn_output}"',
        f'BINDER_CHAIN="{config.binder_chain_id}"',
        f'OMITTED_AA="{omitted_aa_str}"',
        f'NUM_SEQUENCES={config.num_sequences_per_backbone}',
        "",
        "# Create output directory",
        "mkdir -p \"$PROTEINMPNN_OUTPUT\"",
        "",
        "# ProteinMPNN command",
        "# NOTE: Adjust the command below based on your ProteinMPNN installation",
        f'if [ -z "$PROTEINMPNN_PATH" ]; then',
        '    PROTEINMPNN_CMD="python protein_mpnn_run.py"',
        "else",
        '    PROTEINMPNN_CMD="$PROTEINMPNN_PATH"',
        "fi",
        "",
        "# Process each RFdiffusion output",
        "# RFdiffusion generates files with pattern: backbone_*.pdb (e.g., backbone_1_0.pdb, backbone_2_0.pdb)",
        "for complex_pdb in \"$RFDIFFUSION_OUTPUT\"/backbone_*.pdb; do",
        "    if [ ! -f \"$complex_pdb\" ]; then",
        "        echo \"No RFdiffusion outputs found in $RFDIFFUSION_OUTPUT\"",
        "        exit 1",
        "    fi",
        "",
        "    # Extract base name (e.g., complex_1.pdb -> complex_1)",
        '    BASE_NAME=$(basename "$complex_pdb" .pdb)',
        "",
        "    # Run ProteinMPNN for this backbone",
        "    # Format: ProteinMPNN with specified chain and omitted amino acids",
        '    $PROTEINMPNN_CMD \\',
        '        --pdb_path "$complex_pdb" \\',
        f'        --chain_id "$BINDER_CHAIN" \\',
        f'        --num_seq_per_target $NUM_SEQUENCES \\',
        f'        --omit_AAs "$OMITTED_AA" \\',
        '        --out_path "$PROTEINMPNN_OUTPUT/${BASE_NAME}_sequences.fasta"',
        "",
        "    if [ $? -eq 0 ]; then",
        '        echo "✓ Designed sequences for $BASE_NAME"',
        "    else",
        '        echo "✗ Failed to design sequences for $BASE_NAME"',
        "        exit 1",
        "    fi",
        "done",
        "",
        "echo \"",
        f"ProteinMPNN completed: sequences designed for all backbones",
        f"Output directory: $PROTEINMPNN_OUTPUT",
        "\"",
    ]
    
    return "\n".join(script_lines)


def generate_colabfold_script(config: DesignConfig) -> str:
    """
    Generate bash script to run ColabFold for complex structure prediction.
    
    This script:
    - Takes designed sequences from ProteinMPNN
    - Creates FASTA files with target + binder sequences
    - Runs ColabFold batch prediction
    - Organizes outputs per design
    """
    proteinmpnn_output = os.path.join(
        config.output_base_dir,
        "designs",
        config.target_name.replace(" ", "_"),
        config.epitope_id.replace(" ", "_"),
        "proteinmpnn_output"
    )
    colabfold_output = os.path.join(
        config.output_base_dir,
        "designs",
        config.target_name.replace(" ", "_"),
        config.epitope_id.replace(" ", "_"),
        "colabfold_output"
    )
    
    # Get target sequence (would need to extract from PDB in practice)
    # For now, use placeholder
    target_sequence_placeholder = "TARGET_SEQUENCE_FROM_PDB"
    
    script_lines = [
        "#!/bin/bash",
        "# ColabFold Complex Prediction Script",
        f"# Target: {config.target_name} ({config.gene_symbol})",
        f"# Epitope: {config.epitope_id}",
        "",
        "set -e  # Exit on error",
        "",
        "# Configuration",
        f'PROTEINMPNN_OUTPUT="{proteinmpnn_output}"',
        f'COLABFOLD_OUTPUT="{colabfold_output}"',
        f'TARGET_PDB="{config.model_path}"',
        f'TARGET_CHAIN="{config.target_chain_id}"',
        f'BINDER_CHAIN="{config.binder_chain_id}"',
        f'NUM_MODELS={config.num_af_models}',
        f'NUM_RECYCLES={config.num_af_recycles}',
        "",
        "# Create output directory",
        "mkdir -p \"$COLABFOLD_OUTPUT\"",
        "",
        "# Extract target sequence from PDB",
        "# NOTE: In practice, use a tool like BioPython or PyMOL to extract sequence",
        "# For now, this is a placeholder - you'll need to implement sequence extraction",
        f'TARGET_SEQ="{target_sequence_placeholder}"  # TODO: Extract from $TARGET_PDB chain $TARGET_CHAIN',
        "",
        "# ColabFold command",
        "# NOTE: Adjust the command below based on your ColabFold installation",
        f'if [ -z "$COLABFOLD_PATH" ]; then',
        '    COLABFOLD_CMD="colabfold_batch"',
        "else",
        '    COLABFOLD_CMD="$COLABFOLD_PATH"',
        "fi",
        "",
        "# Process each ProteinMPNN output",
        "for seq_file in \"$PROTEINMPNN_OUTPUT\"/*_sequences.fasta; do",
        "    if [ ! -f \"$seq_file\" ]; then",
        "        echo \"No ProteinMPNN outputs found in $PROTEINMPNN_OUTPUT\"",
        "        exit 1",
        "    fi",
        "",
        "    # Extract base name",
        '    BASE_NAME=$(basename "$seq_file" _sequences.fasta)',
        "",
        "    # Create batch FASTA file with target + binder pairs",
        "    # Format: For each binder sequence, create a FASTA entry with target and binder",
        '    BATCH_FASTA="$COLABFOLD_OUTPUT/${BASE_NAME}_batch.fasta"',
        "    > \"$BATCH_FASTA\"  # Clear/create file",
        "",
        "    # Read binder sequences from ProteinMPNN output",
        "    # ProteinMPNN typically outputs sequences in FASTA format",
        "    # Format: >seq_0\\nSEQUENCE\\n>seq_1\\nSEQUENCE...",
        "    seq_count=0",
        "    while IFS= read -r line; do",
        "        if [[ $line =~ ^\\> ]]; then",
        "            # This is a header line, skip or use as identifier",
        "            continue",
        "        elif [[ -n \"$line\" && ! $line =~ ^# ]]; then",
        "            # This is a sequence line",
        "            BINDER_SEQ=\"$line\"",
        "            seq_count=$((seq_count + 1))",
        "",
        "            # Write target:binder pair to batch FASTA",
        "            echo \">${BASE_NAME}_seq${seq_count}_target\" >> \"$BATCH_FASTA\"",
        "            echo \"$TARGET_SEQ\" >> \"$BATCH_FASTA\"",
        "            echo \">${BASE_NAME}_seq${seq_count}_binder\" >> \"$BATCH_FASTA\"",
        "            echo \"$BINDER_SEQ\" >> \"$BATCH_FASTA\"",
        "        fi",
        "    done < \"$seq_file\"",
        "",
        "    if [ $seq_count -eq 0 ]; then",
        "        echo \"Warning: No sequences found in $seq_file\"",
        "        continue",
        "    fi",
        "",
        "    # Run ColabFold batch prediction",
        "    OUTPUT_DIR_DESIGN=\"$COLABFOLD_OUTPUT/${BASE_NAME}\"",
        "    mkdir -p \"$OUTPUT_DIR_DESIGN\"",
        "",
        "    $COLABFOLD_CMD \\",
        "        \"$BATCH_FASTA\" \\",
        "        \"$OUTPUT_DIR_DESIGN\" \\",
        f'        --num-models $NUM_MODELS \\',
        f'        --num-recycles $NUM_RECYCLES \\',
        "        --model-type alphafold2_multimer_v3",
        "",
        "    if [ $? -eq 0 ]; then",
        '        echo "✓ Predicted structures for $BASE_NAME ($seq_count sequences)"',
        "    else",
        '        echo "✗ Failed to predict structures for $BASE_NAME"',
        "        exit 1",
        "    fi",
        "done",
        "",
        "echo \"",
        f"ColabFold completed: structures predicted for all designs",
        f"Output directory: $COLABFOLD_OUTPUT",
        "\"",
    ]
    
    return "\n".join(script_lines)


def generate_alphafold_script(config: DesignConfig) -> str:
    """
    Generate bash script to run AlphaFold-Multimer for complex structure prediction.
    
    Alternative to ColabFold if user prefers AlphaFold-Multimer.
    """
    proteinmpnn_output = os.path.join(
        config.output_base_dir,
        "designs",
        config.target_name.replace(" ", "_"),
        config.epitope_id.replace(" ", "_"),
        "proteinmpnn_output"
    )
    alphafold_output = os.path.join(
        config.output_base_dir,
        "designs",
        config.target_name.replace(" ", "_"),
        config.epitope_id.replace(" ", "_"),
        "alphafold_output"
    )
    
    script_lines = [
        "#!/bin/bash",
        "# AlphaFold-Multimer Complex Prediction Script",
        f"# Target: {config.target_name} ({config.gene_symbol})",
        f"# Epitope: {config.epitope_id}",
        "",
        "set -e  # Exit on error",
        "",
        "# Configuration",
        f'PROTEINMPNN_OUTPUT="{proteinmpnn_output}"',
        f'ALPHAFOLD_OUTPUT="{alphafold_output}"',
        f'TARGET_PDB="{config.model_path}"',
        f'TARGET_CHAIN="{config.target_chain_id}"',
        f'BINDER_CHAIN="{config.binder_chain_id}"',
        f'NUM_MODELS={config.num_af_models}',
        f'NUM_RECYCLES={config.num_af_recycles}',
        "",
        "# Create output directory",
        "mkdir -p \"$ALPHAFOLD_OUTPUT\"",
        "",
        "# AlphaFold-Multimer command",
        "# NOTE: Adjust the command below based on your AlphaFold installation",
        f'if [ -z "$ALPHAFOLD_PATH" ]; then',
        '    ALPHAFOLD_CMD="python run_alphafold.py"',
        "else",
        '    ALPHAFOLD_CMD="$ALPHAFOLD_PATH"',
        "fi",
        "",
        "# Process each ProteinMPNN output",
        "# Similar structure to ColabFold script but using AlphaFold-Multimer API",
        "# This is a template - adjust based on actual AlphaFold-Multimer API",
        "",
        "echo \"",
        "AlphaFold-Multimer script template generated.",
        "Adjust commands based on your AlphaFold-Multimer installation and API.",
        "\"",
    ]
    
    return "\n".join(script_lines)


def generate_master_script(config: DesignConfig) -> str:
    """
    Generate master script that runs the entire pipeline sequentially.
    """
    rfdiffusion_script = os.path.join(
        config.output_base_dir,
        "designs",
        config.target_name.replace(" ", "_"),
        config.epitope_id.replace(" ", "_"),
        "run_rfdiffusion.sh"
    )
    proteinmpnn_script = os.path.join(
        config.output_base_dir,
        "designs",
        config.target_name.replace(" ", "_"),
        config.epitope_id.replace(" ", "_"),
        "run_proteinmpnn.sh"
    )
    if config.use_colabfold:
        prediction_script = os.path.join(
            config.output_base_dir,
            "designs",
            config.target_name.replace(" ", "_"),
            config.epitope_id.replace(" ", "_"),
            "run_colabfold.sh"
        )
    else:
        prediction_script = os.path.join(
            config.output_base_dir,
            "designs",
            config.target_name.replace(" ", "_"),
            config.epitope_id.replace(" ", "_"),
            "run_alphafold.sh"
        )
    
    script_lines = [
        "#!/bin/bash",
        "# Master Pipeline Script",
        f"# Target: {config.target_name} ({config.gene_symbol})",
        f"# Epitope: {config.epitope_id}",
        "",
        "set -e  # Exit on error",
        "",
        "echo \"========================================\"",
        "echo \"Nanobinder Design Pipeline\"",
        f"echo \"Target: {config.target_name}\"",
        f"echo \"Epitope: {config.epitope_id}\"",
        "echo \"========================================\"",
        "",
        "# Step 1: RFdiffusion backbone generation",
        "echo \"",
        "Step 1/3: RFdiffusion Backbone Generation",
        "\"",
        f'bash "{rfdiffusion_script}"',
        "",
        "# Step 2: ProteinMPNN sequence design",
        "echo \"",
        "Step 2/3: ProteinMPNN Sequence Design",
        "\"",
        f'bash "{proteinmpnn_script}"',
        "",
        "# Step 3: Structure prediction",
        "echo \"",
        "Step 3/3: Complex Structure Prediction",
        "\"",
        f'bash "{prediction_script}"',
        "",
        "echo \"",
        "==========================================",
        "Pipeline completed successfully!",
        "==========================================",
        "\"",
    ]
    
    return "\n".join(script_lines)


def save_config_yaml(config: DesignConfig, output_path: str) -> None:
    """Save design configuration as YAML file."""
    config_dict = asdict(config)
    # Convert None to null for YAML
    for key, value in config_dict.items():
        if value is None:
            config_dict[key] = None
    
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w") as f:
        yaml.dump(config_dict, f, default_flow_style=False, sort_keys=False)


def save_config_json(config: DesignConfig, output_path: str) -> None:
    """Save design configuration as JSON file."""
    config_dict = asdict(config)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(config_dict, f, indent=2)


# ============================================================================
# Tool Detection Functions
# ============================================================================

def check_rfdiffusion_available(rfdiffusion_path: Optional[str] = None) -> Tuple[bool, Optional[str], Optional[str]]:
    """
    Check if RFdiffusion is available.
    
    Returns:
        (is_available, executable_path, error_message)
    """
    # Check if path is provided
    if rfdiffusion_path:
        if os.path.exists(rfdiffusion_path):
            return True, rfdiffusion_path, None
        else:
            return False, None, f"RFdiffusion path does not exist: {rfdiffusion_path}"
    
    # Get project root (where exo_gpt module is located)
    # This file is in exo_gpt/, so parent is project root
    project_root = Path(__file__).parent.parent
    cwd = os.getcwd()
    
    # Try common RFdiffusion installation paths (relative and absolute)
    # Prefer script paths over module imports
    common_paths = [
        os.path.join(str(project_root), "RFdiffusion", "scripts", "run_inference.py"),
        os.path.join(str(project_root), "rfdiffusion", "scripts", "run_inference.py"),
        os.path.join(cwd, "RFdiffusion", "scripts", "run_inference.py"),
        os.path.join(cwd, "rfdiffusion", "scripts", "run_inference.py"),
        "RFdiffusion/scripts/run_inference.py",
        "rfdiffusion/scripts/run_inference.py",
        "scripts/run_inference.py",
    ]
    
    for path in common_paths:
        try:
            if os.path.exists(path):
                # Make path absolute for consistency
                abs_path = os.path.abspath(path)
                return True, abs_path, None
        except Exception:
            continue
    
    # Fallback: try module import
    try:
        __import__("rfdiffusion")
        return True, "python -m rfdiffusion", None
    except ImportError:
        pass
    
    # Try checking if it's in PATH
    try:
        result = subprocess.run(
            ["which", "rfdiffusion"],
            capture_output=True,
            text=True,
            timeout=5
        )
        if result.returncode == 0:
            return True, result.stdout.strip(), None
    except Exception:
        pass
    
    return False, None, (
        "RFdiffusion not found. Install from: https://github.com/RosettaCommons/RFdiffusion\n"
        "Installation Instructions (in ExoGPT environment):\n"
        "1. cd RFdiffusion/env/SE3Transformer && pip install -e . && cd ../..\n"
        "2. cd RFdiffusion && pip install -e . --no-deps\n"
        "3. pip install omegaconf hydra-core e3nn wandb dgl\n"
        "4. Download model weights: mkdir -p models && cd models && wget http://files.ipd.uw.edu/pub/RFdiffusion/.../Base_ckpt.pt\n"
        "See RFdiffusion/README.md for full instructions."
    )


def check_proteinmpnn_available(proteinmpnn_path: Optional[str] = None) -> Tuple[bool, Optional[str], Optional[str]]:
    """
    Check if ProteinMPNN is available.
    
    Returns:
        (is_available, executable_path, error_message)
    """
    if proteinmpnn_path:
        if os.path.exists(proteinmpnn_path):
            return True, proteinmpnn_path, None
        else:
            return False, None, f"ProteinMPNN path does not exist: {proteinmpnn_path}"
    
    # Get project root (where exo_gpt module is located)
    # This file is in exo_gpt/, so parent is project root
    project_root = Path(__file__).parent.parent
    cwd = os.getcwd()
    
    # Try common ProteinMPNN installation paths (relative and absolute)
    common_paths = [
        "protein_mpnn_run.py",
        "ProteinMPNN/protein_mpnn_run.py",
        os.path.join(str(project_root), "ProteinMPNN", "protein_mpnn_run.py"),
        os.path.join(str(project_root), "protein_mpnn_run.py"),
        os.path.join(cwd, "ProteinMPNN", "protein_mpnn_run.py"),
        os.path.join(cwd, "protein_mpnn_run.py"),
        "python -m protein_mpnn",
    ]
    
    for path in common_paths:
        try:
            if path.startswith("python"):
                module_name = path.split()[-1]
                try:
                    __import__(module_name)
                    return True, path, None
                except ImportError:
                    continue
            elif os.path.exists(path):
                # Make path absolute for consistency
                abs_path = os.path.abspath(path)
                return True, abs_path, None
        except Exception:
            continue
    
    # Try importing protein_mpnn
    try:
        import protein_mpnn
        return True, "python -m protein_mpnn", None
    except ImportError:
        pass
    
    return False, None, "ProteinMPNN not found. Install from: https://github.com/dauparas/ProteinMPNN"


def check_colabfold_available(colabfold_path: Optional[str] = None) -> Tuple[bool, Optional[str], Optional[str]]:
    """
    Check if ColabFold is available.
    
    Returns:
        (is_available, executable_path, error_message)
    """
    if colabfold_path:
        if os.path.exists(colabfold_path) or colabfold_path == "colabfold_batch":
            return True, colabfold_path, None
        else:
            return False, None, f"ColabFold path does not exist: {colabfold_path}"
    
    # Try common ColabFold commands
    try:
        result = subprocess.run(
            ["colabfold_batch", "--help"],
            capture_output=True,
            text=True,
            timeout=5
        )
        if result.returncode == 0 or "colabfold" in result.stderr.lower() or "colabfold" in result.stdout.lower():
            return True, "colabfold_batch", None
    except FileNotFoundError:
        pass
    except Exception:
        pass
    
    # Try importing colabfold
    try:
        import colabfold
        return True, "colabfold_batch", None
    except ImportError:
        pass
    
    return False, None, "ColabFold not found. Install with: pip install colabfold"


def check_alphafold_available(alphafold_path: Optional[str] = None) -> Tuple[bool, Optional[str], Optional[str]]:
    """
    Check if AlphaFold-Multimer is available.
    
    Returns:
        (is_available, executable_path, error_message)
    """
    if alphafold_path:
        if os.path.exists(alphafold_path):
            return True, alphafold_path, None
        else:
            return False, None, f"AlphaFold path does not exist: {alphafold_path}"
    
    # Try common AlphaFold paths
    common_paths = [
        "run_alphafold.py",
        "alphafold/run_alphafold.py",
    ]
    
    for path in common_paths:
        if os.path.exists(path):
            return True, path, None
    
    return False, None, "AlphaFold-Multimer not found. Install from: https://github.com/deepmind/alphafold"


# ============================================================================
# Tool Execution Functions
# ============================================================================

def map_epitope_residues_to_pdb(
    epitope_residues: List[int],
    pdb_path: str,
    chain_id: str,
    domain_residues: Optional[Tuple[int, int]] = None
) -> List[int]:
    """
    Map epitope residues (UniProt numbering) to PDB residue numbers.
    
    For PDB fragments, PDB residue numbers typically correspond directly to 
    UniProt residue numbers (PDB residue 33 = UniProt residue 33).
    This function validates that epitope residues exist in the PDB file.
    
    Args:
        epitope_residues: List of epitope residue numbers (UniProt numbering)
        pdb_path: Path to PDB file
        chain_id: Chain ID in PDB (e.g., "A")
        domain_residues: Optional (start, end) tuple for domain in UniProt numbering
    
    Returns:
        List of PDB residue numbers that exist in the PDB file
    """
    try:
        from Bio.PDB import PDBParser
        parser = PDBParser(QUIET=True)
        
        # Make path absolute if needed
        if not os.path.isabs(pdb_path):
            project_root = Path(__file__).parent.parent
            pdb_path = os.path.join(str(project_root), pdb_path.lstrip('./'))
            pdb_path = os.path.normpath(pdb_path)
        
        if not os.path.exists(pdb_path):
            # If PDB doesn't exist, return original residues (will fail later with better error)
            return epitope_residues
        
        structure = parser.get_structure("target", pdb_path)
        chain = structure[0][chain_id]
        residues = list(chain.get_residues())
        
        if not residues:
            return epitope_residues
        
        # Get actual PDB residue numbers
        pdb_residue_numbers = set()
        for r in residues:
            if r.id[0] == ' ':  # Only regular residues, not heteroatoms
                pdb_residue_numbers.add(r.id[1])
        
        # For PDB fragments, PDB residue numbers = UniProt residue numbers
        # So we can directly check if epitope residues exist in PDB
        valid_residues = []
        invalid_residues = []
        
        for res in epitope_residues:
            if res in pdb_residue_numbers:
                valid_residues.append(res)
            else:
                invalid_residues.append(res)
        
        if invalid_residues:
            # Log warning but continue with valid residues
            import warnings
            warnings.warn(
                f"Some epitope residues not found in PDB: {invalid_residues}. "
                f"PDB has residues {min(pdb_residue_numbers)}-{max(pdb_residue_numbers)}. "
                f"Using valid residues: {valid_residues}"
            )
        
        if not valid_residues:
            # If no valid residues, return original (will fail with better error message)
            return epitope_residues
        
        return valid_residues
        
    except Exception as e:
        # If parsing fails, return original residues (will fail later with better error)
        import warnings
        warnings.warn(f"Could not parse PDB to validate epitope residues: {str(e)}")
        return epitope_residues


def execute_rfdiffusion(
    config: DesignConfig,
    rfdiffusion_cmd: str,
    progress_callback: Optional[Callable[[str, float], None]] = None
) -> Dict[str, Any]:
    """
    Execute RFdiffusion to generate binder backbones.
    
    Args:
        config: Design configuration
        rfdiffusion_cmd: Command to run RFdiffusion
        progress_callback: Optional callback for progress updates (message, percent)
    
    Returns:
        Dictionary with execution results
    """
    output_dir = os.path.join(
        config.output_base_dir,
        "designs",
        config.target_name.replace(" ", "_"),
        config.epitope_id.replace(" ", "_"),
        "rfdiffusion_output"
    )
    os.makedirs(output_dir, exist_ok=True)
    
    if progress_callback:
        progress_callback(f"RFdiffusion: Now processing target {config.target_name} / {config.epitope_id}", 0.0)
    
    # Validate PDB file exists before starting
    model_path = config.model_path
    if not os.path.isabs(model_path):
        project_root = Path(__file__).parent.parent
        model_path = os.path.join(str(project_root), model_path.lstrip('./'))
        model_path = os.path.normpath(model_path)
    
    if not os.path.exists(model_path):
        error_msg = (
            f"RFdiffusion: PDB file not found: {model_path}\n"
            f"Please download the structure file. For PDB ID {config.pdb_id or 'N/A'}, "
            f"download from: https://www.rcsb.org/structure/{config.pdb_id if config.pdb_id else 'UNKNOWN'}\n"
            f"Or for AlphaFold models, download from: https://alphafold.ebi.ac.uk/entry/{config.uniprot or 'UNKNOWN'}"
        )
        if progress_callback:
            progress_callback(error_msg, 0.0)
        return {
            "success": False,
            "generated_files": [],
            "errors": [error_msg],
            "error_summary": f"PDB file not found: {model_path}"
        }
    
    generated_files = []
    errors = []
    
    # Note: RFdiffusion API varies by version. This is a template that needs adjustment.
    # Common patterns:
    # - Some versions use command-line args
    # - Some versions use Python API
    # - Some versions need specific file formats
    
    for i in range(1, config.num_backbones + 1):
        output_pdb = os.path.join(output_dir, f"complex_{i}.pdb")
        
        if progress_callback:
            progress = (i / config.num_backbones) * 100
            progress_callback(f"RFdiffusion: Generating backbone {i}/{config.num_backbones} for {config.target_name} / {config.epitope_id}", progress)
        
        # RFdiffusion uses Hydra config syntax
        # Format: 'contigmap.contigs=[TARGET_CHAIN/0 BINDER_LENGTH]' 'ppi.hotspot_res=[CHAIN_RES,...]'
        # Example: 'contigmap.contigs=[A1-150/0 70-100]' 'ppi.hotspot_res=[A59,A83,A91]'
        
        # Build hotspot residues string (format: CHAIN_RES,CHAIN_RES,...)
        # Epitope residues are in UniProt numbering, which directly maps to PDB residue numbers for fragments
        # Format from examples: 'ppi.hotspot_res=[A59,A83,A91]'
        # Note: The brackets will be added in the command, so we just need the comma-separated list
        
        # Map epitope residues (UniProt numbering) to PDB residue numbers
        # For PDB fragments, UniProt residue numbers = PDB residue numbers
        # We need to validate that the residues exist in the PDB file
        pdb_epitope_residues = map_epitope_residues_to_pdb(
            config.epitope_residues,
            model_path,
            config.target_chain_id,
            getattr(config, 'domain_residues', None)
        )
        
        if not pdb_epitope_residues:
            error_msg = (
                f"RFdiffusion: No valid epitope residues found in PDB file.\n"
                f"Epitope residues (UniProt): {config.epitope_residues}\n"
                f"Please check that epitope residues correspond to residues in the PDB file."
            )
            errors.append(error_msg)
            if progress_callback:
                progress_callback(error_msg, (i / config.num_backbones) * 100)
            continue  # Skip this backbone
        
        # Build hotspot residue list using validated PDB residue numbers
        # pdb_epitope_residues already contains only residues that exist in the PDB
        hotspot_res_list = [f"{config.target_chain_id}{res}" for res in pdb_epitope_residues]
        
        hotspot_res_str = ','.join(hotspot_res_list)  # Just the values, brackets added in command
        
        # Build contig string: [TARGET_REGION/0 BINDER_LENGTH]
        # Format: [CHAIN_START-END/0 BINDER_MIN-BINDER_MAX]
        # We'll use a range for binder length to allow some variation
        binder_min = max(70, config.binder_length - 10)
        binder_max = config.binder_length + 10
        
        # Get actual residue numbers from PDB file (to handle gaps in numbering)
        # PDB files may have gaps (e.g., residues 60, 62, 63 exist but 61 doesn't)
        # We need to use only residues that actually exist
        target_start = 1
        target_end = 150
        actual_residue_numbers = []
        try:
            from Bio.PDB import PDBParser
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure("target", model_path)
            chain = structure[0][config.target_chain_id]
            residues = list(chain.get_residues())
            if residues:
                # Get actual residue numbers (PDB numbering, not sequential)
                actual_residue_numbers = sorted([r.id[1] for r in residues if r.id[0] == ' '])  # Only regular residues, not heteroatoms
                if actual_residue_numbers:
                    target_start = min(actual_residue_numbers)
                    target_end = max(actual_residue_numbers)
                    # Log for debugging
                    if progress_callback:
                        progress_callback(f"RFdiffusion: PDB chain {config.target_chain_id} residue range: {target_start}-{target_end} (actual residues: {len(actual_residue_numbers)} residues)", progress)
        except Exception as e:
            # Fallback: use a reasonable default
            # If we can't parse, assume starting at 1
            target_start = 1
            target_end = 150
            if progress_callback:
                progress_callback(f"RFdiffusion: Warning - Could not parse PDB residue numbers, using default range 1-150: {str(e)}", progress)
        
        # Check if there are gaps in residue numbering
        # If there are gaps, we need to list individual residues instead of using a range
        if actual_residue_numbers:
            expected_residues = set(range(target_start, target_end + 1))
            actual_residues_set = set(actual_residue_numbers)
            missing_residues = expected_residues - actual_residues_set
            
            if missing_residues:
                # There are gaps - need to build contig with individual residue ranges
                # Group consecutive residues into ranges
                residue_ranges = []
                start = actual_residue_numbers[0]
                end = actual_residue_numbers[0]
                
                for i in range(1, len(actual_residue_numbers)):
                    if actual_residue_numbers[i] == end + 1:
                        # Consecutive, extend range
                        end = actual_residue_numbers[i]
                    else:
                        # Gap found, save current range and start new one
                        if start == end:
                            residue_ranges.append(f"{config.target_chain_id}{start}")
                        else:
                            residue_ranges.append(f"{config.target_chain_id}{start}-{end}")
                        start = actual_residue_numbers[i]
                        end = actual_residue_numbers[i]
                
                # Add last range
                if start == end:
                    residue_ranges.append(f"{config.target_chain_id}{start}")
                else:
                    residue_ranges.append(f"{config.target_chain_id}{start}-{end}")
                
                # Build contig with "/" separated ranges (RFdiffusion format for multiple receptor fragments)
                # Format: A33-60/A62-149/0 80-100 (multiple fragments separated by /)
                contig_value = f"{'/'.join(residue_ranges)}/0 {binder_min}-{binder_max}"
                if progress_callback:
                    progress_callback(f"RFdiffusion: PDB has gaps, using contig: {contig_value}", progress)
            else:
                # No gaps, use simple range
                contig_value = f"{config.target_chain_id}{target_start}-{target_end}/0 {binder_min}-{binder_max}"
        else:
            # Fallback: use simple range
            contig_value = f"{config.target_chain_id}{target_start}-{target_end}/0 {binder_min}-{binder_max}"
        
        # Determine script path
        is_python_script = rfdiffusion_cmd.endswith('.py') or (os.path.isfile(rfdiffusion_cmd) and rfdiffusion_cmd.endswith('.py'))
        cmd_parts = rfdiffusion_cmd.split()
        
        if is_python_script:
            script_path = rfdiffusion_cmd
        elif len(cmd_parts) > 1 and cmd_parts[0] == "python":
            script_path = cmd_parts[1]
        else:
            # Default to RFdiffusion scripts path
            project_root = Path(__file__).parent.parent
            script_path = os.path.join(str(project_root), "RFdiffusion", "scripts", "run_inference.py")
        
        # Build Hydra command
        # Output prefix (without .pdb extension - RFdiffusion adds it)
        output_prefix = os.path.join(output_dir, f"backbone_{i}")
        
        # Hydra syntax: contigmap.contigs=[VALUE] where VALUE contains the contig specification
        # According to RFdiffusion README: "the entire argument MUST be enclosed in '' so that 
        # the commandline does not attempt to parse any of the special characters"
        # However, when using subprocess.run() with a list, we don't use shell quotes.
        # The brackets are part of the value that Hydra will parse.
        # Format: contigmap.contigs=[CHAIN_START-END/0 BINDER_MIN-BINDER_MAX]
        contig_arg = f"contigmap.contigs=[{contig_value}]"
        hotspot_arg = f"ppi.hotspot_res=[{hotspot_res_str}]"
        
        # Make model_path absolute if it's relative (needed when cwd is RFdiffusion root)
        model_path = config.model_path
        if not os.path.isabs(model_path):
            # If relative, make it absolute relative to project root
            project_root = Path(__file__).parent.parent
            model_path = os.path.join(str(project_root), model_path.lstrip('./'))
            model_path = os.path.normpath(model_path)
        
        # Validate that the PDB file exists
        if not os.path.exists(model_path):
            error_msg = (
                f"RFdiffusion: PDB file not found: {model_path}\n"
                f"Please ensure the structure file exists. For PDB files, you may need to download them from RCSB PDB.\n"
                f"For AlphaFold models, you may need to download them from AlphaFold DB."
            )
            errors.append(error_msg)
            if progress_callback:
                progress_callback(error_msg, (i / config.num_backbones) * 100)
            continue  # Skip this backbone
        
        cmd = [
            "python", script_path,
            f"inference.input_pdb={model_path}",
            f"inference.output_prefix={output_prefix}",
            f"inference.num_designs=1",  # One design per iteration
            contig_arg,
            hotspot_arg
        ]
        
        try:
            # Set working directory and environment for RFdiffusion
            cwd = None
            env = os.environ.copy()
            
            # DGL with CUDA 11.6 support (dgl-cuda11.6) is installed, so graphs will use GPU
            # No patching needed - DGL will automatically use CUDA when tensors are on GPU
            
            # Add conda lib directory to LD_LIBRARY_PATH for DGL to find CUDA libraries
            conda_prefix = os.environ.get("CONDA_PREFIX", "")
            if conda_prefix:
                lib_path = os.path.join(conda_prefix, "lib")
                current_ld_path = env.get("LD_LIBRARY_PATH", "")
                if current_ld_path:
                    env["LD_LIBRARY_PATH"] = f"{lib_path}:{current_ld_path}"
                else:
                    env["LD_LIBRARY_PATH"] = lib_path
            
            # Note: PyTorch 2.1.0 + CUDA 11.8 is installed and working
            # DGL is forced to CPU mode to avoid CUDA compatibility issues
            # PyTorch will still use GPU for most operations
            
            # Set working directory to RFdiffusion root (required for RFdiffusion)
            if is_python_script:
                script_dir = os.path.dirname(os.path.abspath(script_path))
                # RFdiffusion root is typically the parent of scripts/ directory
                rfdiffusion_root = os.path.dirname(script_dir)  # Go up from scripts/ to RFdiffusion/
                
                if os.path.exists(rfdiffusion_root) and os.path.basename(rfdiffusion_root) in ["RFdiffusion", "rfdiffusion"]:
                    # Set working directory to RFdiffusion root
                    cwd = rfdiffusion_root
                    
                    # Add RFdiffusion directory to PYTHONPATH so 'rfdiffusion' module can be found
                    current_pythonpath = env.get("PYTHONPATH", "")
                    if current_pythonpath:
                        env["PYTHONPATH"] = f"{rfdiffusion_root}:{current_pythonpath}"
                    else:
                        env["PYTHONPATH"] = rfdiffusion_root
                else:
                    # Fallback: use script directory
                    cwd = script_dir
            
            # Use the current Python environment (ExoGPT) - no need for conda run
            # RFdiffusion and dependencies are installed in the same environment as the GUI
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=3600,  # 1 hour timeout
                check=False,
                cwd=cwd,
                env=env,
                shell=False  # Explicitly set to False to avoid shell interpretation
            )
            
            # RFdiffusion outputs files with pattern: {output_prefix}_*.pdb
            # Find the generated PDB file
            output_files = []
            if result.returncode == 0:
                # Look for generated PDB files
                import glob
                pattern = os.path.join(output_dir, f"backbone_{i}_*.pdb")
                output_files = glob.glob(pattern)
                if not output_files:
                    # Try alternative pattern
                    pattern = os.path.join(output_dir, f"*{i}*.pdb")
                    output_files = glob.glob(pattern)
                
                if output_files:
                    generated_files.extend(output_files)
                else:
                    # Check if output_pdb exists (fallback)
                    if os.path.exists(output_pdb):
                        generated_files.append(output_pdb)
            else:
                # Include both stdout and stderr in error message for debugging
                error_details = []
                if result.stdout:
                    # Show last 1000 chars of stdout (most recent messages)
                    stdout_msg = result.stdout[-1000:] if len(result.stdout) > 1000 else result.stdout
                    error_details.append(f"stdout (last 1000 chars): {stdout_msg}")
                if result.stderr:
                    # Show full stderr to see complete traceback
                    stderr_msg = result.stderr
                    error_details.append(f"stderr: {stderr_msg}")
                    
                    # Check for common dependency errors and provide helpful suggestions
                    if "ModuleNotFoundError" in stderr_msg or "No module named" in stderr_msg:
                        # Extract the missing module name from the error
                        missing_module = None
                        if "No module named '" in stderr_msg:
                            try:
                                missing_module = stderr_msg.split("No module named '")[1].split("'")[0]
                            except:
                                pass
                        elif "No module named " in stderr_msg:
                            try:
                                missing_module = stderr_msg.split("No module named ")[1].split()[0].strip()
                            except:
                                pass
                        
                        if missing_module:
                            if missing_module == "torch":
                                error_details.append(f"SUGGESTION: Missing module '{missing_module}'. Install with: pip install torch")
                            elif missing_module in ["omegaconf", "hydra"]:
                                error_details.append(f"SUGGESTION: Missing module '{missing_module}'. Install with: pip install omegaconf hydra-core")
                            elif missing_module == "rfdiffusion":
                                error_details.append(f"SUGGESTION: Missing module '{missing_module}'. Install with: cd RFdiffusion && pip install -e . --no-deps")
                                error_details.append("NOTE: RFdiffusion should be installed in the ExoGPT environment. See INSTALL_RFDIFFUSION.md for setup instructions.")
                            elif missing_module == "pyrsistent":
                                error_details.append(f"SUGGESTION: Missing module '{missing_module}'. Install with: pip install pyrsistent")
                            elif missing_module in ["dgl", "e3nn", "wandb"]:
                                error_details.append(f"SUGGESTION: Missing module '{missing_module}'. Install with: pip install {missing_module}")
                            else:
                                error_details.append(f"SUGGESTION: Missing module '{missing_module}'. Install with: pip install {missing_module}")
                                error_details.append("NOTE: If this is an RFdiffusion dependency, see INSTALL_RFDIFFUSION.md for complete setup instructions.")
                        else:
                            # Fallback: check for common patterns
                            if "torch" in stderr_msg.lower():
                                error_details.append("SUGGESTION: RFdiffusion requires PyTorch. Install with: pip install torch")
                            elif "omegaconf" in stderr_msg.lower() or "hydra" in stderr_msg.lower():
                                error_details.append("SUGGESTION: Install with: pip install omegaconf hydra-core")
                            elif "rfdiffusion" in stderr_msg.lower():
                                error_details.append("SUGGESTION: RFdiffusion module not found. Install with: cd RFdiffusion && pip install -e . --no-deps")
                            else:
                                error_details.append("SUGGESTION: RFdiffusion requires dependencies in the ExoGPT environment.")
                                error_details.append("Setup: pip install omegaconf hydra-core e3nn wandb dgl<2.0 pyrsistent")
                                error_details.append("Then: cd RFdiffusion/env/SE3Transformer && pip install -e . && cd ../.. && pip install -e . --no-deps")
                                error_details.append("See INSTALL_RFDIFFUSION.md for complete installation instructions.")
                if not error_details:
                    error_details.append(f"Return code: {result.returncode}")
                
                error_msg = f"RFdiffusion: Failed to generate backbone {i}/{config.num_backbones} for {config.target_name} / {config.epitope_id}: {'; '.join(error_details)}"
                errors.append(error_msg)
                if progress_callback:
                    progress_callback(error_msg, (i / config.num_backbones) * 100)
        except subprocess.TimeoutExpired:
            error_msg = f"RFdiffusion: Timed out generating backbone {i}/{config.num_backbones} for {config.target_name} / {config.epitope_id}"
            errors.append(error_msg)
            if progress_callback:
                progress_callback(error_msg, (i / config.num_backbones) * 100)
        except Exception as e:
            error_msg = f"RFdiffusion: Error generating backbone {i}/{config.num_backbones} for {config.target_name} / {config.epitope_id}: {str(e)}"
            errors.append(error_msg)
            if progress_callback:
                progress_callback(error_msg, (i / config.num_backbones) * 100)
    
    if progress_callback:
        progress_callback(f"RFdiffusion: Completed {len(generated_files)}/{config.num_backbones} backbones for {config.target_name} / {config.epitope_id}", 100.0)
    
    result = {
        "success": len(generated_files) > 0,
        "generated_files": generated_files,
        "errors": errors,
        "output_dir": output_dir
    }
    
    # If failed, add helpful message about dependencies
    if not result["success"] and errors:
        # Check if it's a dependency issue
        dependency_errors = [e for e in errors if "ModuleNotFoundError" in e or "No module named" in e]
        if dependency_errors:
            result["error_summary"] = (
                f"RFdiffusion failed due to missing dependencies. "
                f"Install required packages: pip install torch omegaconf hydra-core. "
                f"Or use the generated scripts after setting up RFdiffusion environment. "
                f"See RFdiffusion/README.md for full installation instructions."
            )
        else:
            result["error_summary"] = f"RFdiffusion failed: {len(errors)} error(s). First error: {errors[0][:200]}"
    
    return result


def execute_proteinmpnn(
    config: DesignConfig,
    proteinmpnn_cmd: str,
    rfdiffusion_output_dir: str,
    progress_callback: Optional[Callable[[str, float], None]] = None
) -> Dict[str, Any]:
    """
    Execute ProteinMPNN to design sequences.
    
    Args:
        config: Design configuration
        proteinmpnn_cmd: Command to run ProteinMPNN
        rfdiffusion_output_dir: Directory with RFdiffusion output PDBs
        progress_callback: Optional callback for progress updates
    
    Returns:
        Dictionary with execution results
    """
    output_dir = os.path.join(
        config.output_base_dir,
        "designs",
        config.target_name.replace(" ", "_"),
        config.epitope_id.replace(" ", "_"),
        "proteinmpnn_output"
    )
    os.makedirs(output_dir, exist_ok=True)
    
    # Find RFdiffusion output files
    # RFdiffusion generates files with pattern: backbone_*.pdb (e.g., backbone_1_0.pdb, backbone_2_0.pdb)
    import glob
    complex_pdbs = glob.glob(os.path.join(rfdiffusion_output_dir, "backbone_*.pdb"))
    
    if not complex_pdbs:
        return {
            "success": False,
            "error": f"No RFdiffusion outputs found in {rfdiffusion_output_dir}",
            "generated_files": []
        }
    
    # Extract target name from config if available (passed through complex_pdbs path structure)
    target_info = "target"
    if complex_pdbs:
        # Try to extract target/epitope from path: .../designs/TARGET/EPITOPE/...
        path_parts = complex_pdbs[0].split(os.sep)
        if "designs" in path_parts:
            designs_idx = path_parts.index("designs")
            if designs_idx + 1 < len(path_parts):
                target_info = path_parts[designs_idx + 1]
            if designs_idx + 2 < len(path_parts):
                target_info += f" / {path_parts[designs_idx + 2]}"
    
    if progress_callback:
        progress_callback(f"ProteinMPNN: Now processing target {target_info} ({len(complex_pdbs)} backbones)", 0.0)
    
    generated_files = []
    errors = []
    omitted_aa_str = ",".join(config.omitted_aa)
    
    for idx, complex_pdb in enumerate(complex_pdbs):
        base_name = os.path.splitext(os.path.basename(complex_pdb))[0]
        output_fasta = os.path.join(output_dir, f"{base_name}_sequences.fasta")
        
        if progress_callback:
            progress = ((idx + 1) / len(complex_pdbs)) * 100
            progress_callback(f"ProteinMPNN: Designing sequences for backbone {idx+1}/{len(complex_pdbs)} of {target_info}", progress)
        
        # Build command
        cmd_parts = proteinmpnn_cmd.split()
        if len(cmd_parts) > 1 and cmd_parts[0] == "python":
            script_path = cmd_parts[1] if len(cmd_parts) > 1 else "protein_mpnn_run.py"
            cmd = [
                "python", script_path,
                "--pdb_path", complex_pdb,
                "--chain_id", config.binder_chain_id,
                "--num_seq_per_target", str(config.num_sequences_per_backbone),
                "--omit_AAs", omitted_aa_str,
                "--out_path", output_fasta
            ]
        else:
            cmd = [
                proteinmpnn_cmd,
                "--pdb_path", complex_pdb,
                "--chain_id", config.binder_chain_id,
                "--num_seq_per_target", str(config.num_sequences_per_backbone),
                "--omit_AAs", omitted_aa_str,
                "--out_path", output_fasta
            ]
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=600,  # 10 minute timeout
                check=False
            )
            
            if result.returncode == 0 and os.path.exists(output_fasta):
                generated_files.append(output_fasta)
            else:
                error_msg = f"ProteinMPNN failed for {base_name}: {result.stderr}"
                errors.append(error_msg)
        except subprocess.TimeoutExpired:
            errors.append(f"ProteinMPNN timed out for {base_name}")
        except Exception as e:
            errors.append(f"ProteinMPNN error for {base_name}: {str(e)}")
    
    if progress_callback:
        progress_callback(f"ProteinMPNN: Completed {len(generated_files)}/{len(complex_pdbs)} sequence files for {target_info}", 100.0)
    
    return {
        "success": len(generated_files) > 0,
        "generated_files": generated_files,
        "errors": errors,
        "output_dir": output_dir
    }


def execute_colabfold(
    config: DesignConfig,
    colabfold_cmd: str,
    proteinmpnn_output_dir: str,
    progress_callback: Optional[Callable[[str, float], None]] = None
) -> Dict[str, Any]:
    """
    Execute ColabFold for complex structure prediction.
    
    Args:
        config: Design configuration
        colabfold_cmd: Command to run ColabFold
        proteinmpnn_output_dir: Directory with ProteinMPNN output sequences
        progress_callback: Optional callback for progress updates
    
    Returns:
        Dictionary with execution results
    """
    output_dir = os.path.join(
        config.output_base_dir,
        "designs",
        config.target_name.replace(" ", "_"),
        config.epitope_id.replace(" ", "_"),
        "colabfold_output"
    )
    os.makedirs(output_dir, exist_ok=True)
    
    # Extract target sequence from PDB (using BioPython if available)
    target_seq = None
    try:
        from Bio.PDB import PDBParser
        from Bio.SeqUtils import seq1
        
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("target", config.model_path)
        chain = structure[0][config.target_chain_id]
        target_seq = seq1("".join([residue.resname for residue in chain]))
    except Exception:
        # Fallback: try to extract sequence another way or use placeholder
        if progress_callback:
            progress_callback("Warning: Could not extract target sequence from PDB. Using placeholder.", 0.0)
        target_seq = "PLACEHOLDER_SEQUENCE"  # User needs to replace this
    
    # Find ProteinMPNN output files
    import glob
    seq_files = glob.glob(os.path.join(proteinmpnn_output_dir, "*_sequences.fasta"))
    
    if not seq_files:
        return {
            "success": False,
            "error": f"No ProteinMPNN outputs found in {proteinmpnn_output_dir}",
            "generated_files": []
        }
    
    # Extract target name from config if available (passed through seq_files path structure)
    target_info = "target"
    if seq_files:
        # Try to extract target/epitope from path: .../designs/TARGET/EPITOPE/...
        path_parts = seq_files[0].split(os.sep)
        if "designs" in path_parts:
            designs_idx = path_parts.index("designs")
            if designs_idx + 1 < len(path_parts):
                target_info = path_parts[designs_idx + 1]
            if designs_idx + 2 < len(path_parts):
                target_info += f" / {path_parts[designs_idx + 2]}"
    
    if progress_callback:
        progress_callback(f"ColabFold: Now processing target {target_info} ({len(seq_files)} sequences)", 0.0)
    
    generated_files = []
    errors = []
    
    for idx, seq_file in enumerate(seq_files):
        base_name = os.path.basename(seq_file).replace("_sequences.fasta", "")
        design_output_dir = os.path.join(output_dir, base_name)
        os.makedirs(design_output_dir, exist_ok=True)
        
        if progress_callback:
            progress = ((idx + 1) / len(seq_files)) * 100
            progress_callback(f"ColabFold: Predicting structure for sequence {idx+1}/{len(seq_files)} of {target_info}", progress)
        
        # Read sequences and create batch FASTA
        batch_fasta = os.path.join(design_output_dir, f"{base_name}_batch.fasta")
        try:
            with open(seq_file, "r") as f:
                sequences = f.read()
            
            # Parse sequences and create target:binder pairs
            with open(batch_fasta, "w") as f:
                seq_count = 0
                current_seq = ""
                for line in sequences.split("\n"):
                    if line.startswith(">"):
                        if current_seq:
                            # Write previous sequence pair
                            f.write(f">{base_name}_seq{seq_count}_target\n")
                            f.write(f"{target_seq}\n")
                            f.write(f">{base_name}_seq{seq_count}_binder\n")
                            f.write(f"{current_seq}\n")
                            seq_count += 1
                            current_seq = ""
                    elif line.strip():
                        current_seq += line.strip()
                
                # Write last sequence
                if current_seq:
                    f.write(f">{base_name}_seq{seq_count}_target\n")
                    f.write(f"{target_seq}\n")
                    f.write(f">{base_name}_seq{seq_count}_binder\n")
                    f.write(f"{current_seq}\n")
        except Exception as e:
            errors.append(f"Failed to create batch FASTA for {base_name}: {str(e)}")
            continue
        
        # Run ColabFold
        cmd = [
            colabfold_cmd,
            batch_fasta,
            design_output_dir,
            "--num-models", str(config.num_af_models),
            "--num-recycles", str(config.num_af_recycles),
            "--model-type", "alphafold2_multimer_v3"
        ]
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=7200,  # 2 hour timeout
                check=False
            )
            
            if result.returncode == 0:
                # Check for output files
                output_files = glob.glob(os.path.join(design_output_dir, "*.pdb"))
                if output_files:
                    generated_files.extend(output_files)
                else:
                    errors.append(f"No output files generated for {base_name}")
            else:
                errors.append(f"ColabFold failed for {base_name}: {result.stderr}")
        except subprocess.TimeoutExpired:
            errors.append(f"ColabFold timed out for {base_name}")
        except Exception as e:
            errors.append(f"ColabFold error for {base_name}: {str(e)}")
    
    if progress_callback:
        progress_callback(f"ColabFold: Completed {len(generated_files)} structure predictions for {target_info}", 100.0)
    
    return {
        "success": len(generated_files) > 0,
        "generated_files": generated_files,
        "errors": errors,
        "output_dir": output_dir
    }


def generate_design_plan(
    step2_json_path: str,
    output_base_dir: str = "./workflows",
    binder_length: int = 90,
    num_backbones: int = 10,
    num_sequences_per_backbone: int = 8,
    num_af_models: int = 5,
    num_af_recycles: int = 3,
    use_colabfold: bool = True,
    rfdiffusion_path: Optional[str] = None,
    proteinmpnn_path: Optional[str] = None,
    colabfold_path: Optional[str] = None,
    alphafold_path: Optional[str] = None,
    execute: bool = False,
    progress_callback: Optional[Callable[[str, float], None]] = None,
) -> Dict[str, Any]:
    """
    Generate complete design plan with all scripts and configs.
    Optionally execute tools if available.
    
    Args:
        execute: If True, check for tools and execute if available
        progress_callback: Optional callback for progress updates (message, percent)
    
    Returns:
        Dictionary with design plan summary and execution results
    """
    # Load Step 2 data
    step2_data = load_step2_json(step2_json_path)
    
    # Create design configs
    configs = create_design_configs(
        step2_data,
        output_base_dir=output_base_dir,
        binder_length=binder_length,
        num_backbones=num_backbones,
        num_sequences_per_backbone=num_sequences_per_backbone,
        num_af_models=num_af_models,
        num_af_recycles=num_af_recycles,
        use_colabfold=use_colabfold,
        rfdiffusion_path=rfdiffusion_path,
        proteinmpnn_path=proteinmpnn_path,
        colabfold_path=colabfold_path,
        alphafold_path=alphafold_path,
    )
    
    if not configs:
        return {
            "success": False,
            "error": "No design configurations generated. Check Step 2 output for curated targets with epitopes.",
        }
    
    # Check tool availability if execution is requested
    tool_status = {}
    execution_results = {}
    
    if execute:
        if progress_callback:
            progress_callback("Checking tool availability...", 0.0)
        
        # Check RFdiffusion
        rf_available, rf_cmd, rf_error = check_rfdiffusion_available(rfdiffusion_path)
        tool_status["rfdiffusion"] = {
            "available": rf_available,
            "command": rf_cmd,
            "error": rf_error
        }
        
        # Check ProteinMPNN
        mpnn_available, mpnn_cmd, mpnn_error = check_proteinmpnn_available(proteinmpnn_path)
        tool_status["proteinmpnn"] = {
            "available": mpnn_available,
            "command": mpnn_cmd,
            "error": mpnn_error
        }
        
        # Check ColabFold or AlphaFold
        if use_colabfold:
            cf_available, cf_cmd, cf_error = check_colabfold_available(colabfold_path)
            tool_status["colabfold"] = {
                "available": cf_available,
                "command": cf_cmd,
                "error": cf_error
            }
        else:
            af_available, af_cmd, af_error = check_alphafold_available(alphafold_path)
            tool_status["alphafold"] = {
                "available": af_available,
                "command": af_cmd,
                "error": af_error
            }
    
    # Generate scripts and configs for each design
    generated_files = []
    
    for config in configs:
        # Create output directory structure
        design_dir = os.path.join(
            output_base_dir,
            "designs",
            config.target_name.replace(" ", "_"),
            config.epitope_id.replace(" ", "_")
        )
        os.makedirs(design_dir, exist_ok=True)
        
        # Save config files
        config_yaml_path = os.path.join(design_dir, "config.yaml")
        config_json_path = os.path.join(design_dir, "config.json")
        save_config_yaml(config, config_yaml_path)
        save_config_json(config, config_json_path)
        generated_files.append(config_yaml_path)
        generated_files.append(config_json_path)
        
        # Generate scripts
        rfdiffusion_script_path = os.path.join(design_dir, "run_rfdiffusion.sh")
        proteinmpnn_script_path = os.path.join(design_dir, "run_proteinmpnn.sh")
        
        with open(rfdiffusion_script_path, "w") as f:
            f.write(generate_rfdiffusion_script(config))
        os.chmod(rfdiffusion_script_path, 0o755)
        generated_files.append(rfdiffusion_script_path)
        
        with open(proteinmpnn_script_path, "w") as f:
            f.write(generate_proteinmpnn_script(config))
        os.chmod(proteinmpnn_script_path, 0o755)
        generated_files.append(proteinmpnn_script_path)
        
        if config.use_colabfold:
            colabfold_script_path = os.path.join(design_dir, "run_colabfold.sh")
            with open(colabfold_script_path, "w") as f:
                f.write(generate_colabfold_script(config))
            os.chmod(colabfold_script_path, 0o755)
            generated_files.append(colabfold_script_path)
        else:
            alphafold_script_path = os.path.join(design_dir, "run_alphafold.sh")
            with open(alphafold_script_path, "w") as f:
                f.write(generate_alphafold_script(config))
            os.chmod(alphafold_script_path, 0o755)
            generated_files.append(alphafold_script_path)
        
        # Generate master script
        master_script_path = os.path.join(design_dir, "run_pipeline.sh")
        with open(master_script_path, "w") as f:
            f.write(generate_master_script(config))
        os.chmod(master_script_path, 0o755)
        generated_files.append(master_script_path)
        
        # Execute tools if requested and available
        if execute:
            design_key = f"{config.target_name}_{config.epitope_id}"
            execution_results[design_key] = {}
            
            # Step 1: Execute RFdiffusion
            if tool_status.get("rfdiffusion", {}).get("available"):
                if progress_callback:
                    progress_callback(f"RFdiffusion: Starting execution for {config.target_name} / {config.epitope_id}...", 10.0)
                rf_result = execute_rfdiffusion(
                    config,
                    tool_status["rfdiffusion"]["command"],
                    progress_callback
                )
                execution_results[design_key]["rfdiffusion"] = rf_result
            else:
                execution_results[design_key]["rfdiffusion"] = {
                    "success": False,
                    "error": tool_status.get("rfdiffusion", {}).get("error", "RFdiffusion not available")
                }
            
            # Step 2: Execute ProteinMPNN (if RFdiffusion succeeded)
            if tool_status.get("proteinmpnn", {}).get("available"):
                rf_output_dir = execution_results[design_key].get("rfdiffusion", {}).get("output_dir")
                if rf_output_dir and execution_results[design_key].get("rfdiffusion", {}).get("success"):
                    if progress_callback:
                        progress_callback(f"ProteinMPNN: Starting execution for {config.target_name} / {config.epitope_id}...", 50.0)
                    mpnn_result = execute_proteinmpnn(
                        config,
                        tool_status["proteinmpnn"]["command"],
                        rf_output_dir,
                        progress_callback
                    )
                    execution_results[design_key]["proteinmpnn"] = mpnn_result
                else:
                    execution_results[design_key]["proteinmpnn"] = {
                        "success": False,
                        "error": "RFdiffusion must succeed before ProteinMPNN"
                    }
            else:
                execution_results[design_key]["proteinmpnn"] = {
                    "success": False,
                    "error": tool_status.get("proteinmpnn", {}).get("error", "ProteinMPNN not available")
                }
            
            # Step 3: Execute ColabFold/AlphaFold (if ProteinMPNN succeeded)
            if use_colabfold:
                cf_available = tool_status.get("colabfold", {}).get("available")
                if cf_available:
                    mpnn_output_dir = execution_results[design_key].get("proteinmpnn", {}).get("output_dir")
                    if mpnn_output_dir and execution_results[design_key].get("proteinmpnn", {}).get("success"):
                        if progress_callback:
                            progress_callback(f"ColabFold: Starting execution for {config.target_name} / {config.epitope_id}...", 80.0)
                        cf_result = execute_colabfold(
                            config,
                            tool_status["colabfold"]["command"],
                            mpnn_output_dir,
                            progress_callback
                        )
                        execution_results[design_key]["colabfold"] = cf_result
                    else:
                        execution_results[design_key]["colabfold"] = {
                            "success": False,
                            "error": "ProteinMPNN must succeed before ColabFold"
                        }
                else:
                    execution_results[design_key]["colabfold"] = {
                        "success": False,
                        "error": tool_status.get("colabfold", {}).get("error", "ColabFold not available")
                    }
            # Similar for AlphaFold if not using ColabFold
    
    # Create summary
    design_summary = []
    for config in configs:
        design_summary.append({
            "target_name": config.target_name,
            "gene_symbol": config.gene_symbol,
            "epitope_id": config.epitope_id,
            "epitope_residues": config.epitope_residues,
            "output_directory": os.path.join(
                output_base_dir,
                "designs",
                config.target_name.replace(" ", "_"),
                config.epitope_id.replace(" ", "_")
            ),
            "binder_length": config.binder_length,
            "num_backbones": config.num_backbones,
            "num_sequences_per_backbone": config.num_sequences_per_backbone,
            "total_sequences": config.num_backbones * config.num_sequences_per_backbone,
            "num_af_models": config.num_af_models,
            "total_af_predictions": config.num_backbones * config.num_sequences_per_backbone * config.num_af_models,
        })
    
    result = {
        "success": True,
        "disease": step2_data.get("disease", "Unknown"),
        "biofluid": step2_data.get("biofluid", "Unknown"),
        "num_designs": len(configs),
        "output_base_dir": output_base_dir,
        "designs": design_summary,
        "generated_files": generated_files,
        "pipeline_steps": [
            "1. RFdiffusion: Generate binder backbones",
            "2. ProteinMPNN: Design sequences for backbones",
            f"3. {'ColabFold' if use_colabfold else 'AlphaFold-Multimer'}: Predict complex structures"
        ],
    }
    
    # Add execution information if execution was attempted
    if execute:
        result["execution_mode"] = True
        result["tool_status"] = tool_status
        result["execution_results"] = execution_results
    else:
        result["execution_mode"] = False
    
    return result


def run_cli() -> None:
    """CLI entrypoint for Step 3."""
    parser = argparse.ArgumentParser(description="Step 3 — Nanobinder Design Orchestrator")
    parser.add_argument(
        "--step2_json",
        required=True,
        help="Path to Step 2 JSON output (curated targets with epitopes)"
    )
    parser.add_argument(
        "--out_plan",
        required=True,
        help="Path to write design plan JSON"
    )
    parser.add_argument(
        "--output_base_dir",
        default="./workflows",
        help="Base directory for all outputs (default: ./workflows)"
    )
    parser.add_argument(
        "--binder_length",
        type=int,
        default=90,
        help="Binder length in amino acids (default: 90)"
    )
    parser.add_argument(
        "--num_backbones",
        type=int,
        default=10,
        help="Number of RFdiffusion backbone designs (default: 10)"
    )
    parser.add_argument(
        "--num_sequences_per_backbone",
        type=int,
        default=8,
        help="Number of ProteinMPNN sequences per backbone (default: 8)"
    )
    parser.add_argument(
        "--num_af_models",
        type=int,
        default=5,
        help="Number of AlphaFold models per sequence (default: 5)"
    )
    parser.add_argument(
        "--num_af_recycles",
        type=int,
        default=3,
        help="Number of AlphaFold recycles (default: 3)"
    )
    parser.add_argument(
        "--use_colabfold",
        action="store_true",
        default=True,
        help="Use ColabFold (default: True). Use --no-use_colabfold for AlphaFold-Multimer"
    )
    parser.add_argument(
        "--no-use_colabfold",
        dest="use_colabfold",
        action="store_false",
        help="Use AlphaFold-Multimer instead of ColabFold"
    )
    parser.add_argument(
        "--rfdiffusion_path",
        default=None,
        help="Path to RFdiffusion executable/script (optional)"
    )
    parser.add_argument(
        "--proteinmpnn_path",
        default=None,
        help="Path to ProteinMPNN executable/script (optional)"
    )
    parser.add_argument(
        "--colabfold_path",
        default=None,
        help="Path to ColabFold executable/script (optional)"
    )
    parser.add_argument(
        "--alphafold_path",
        default=None,
        help="Path to AlphaFold-Multimer executable/script (optional)"
    )
    parser.add_argument(
        "--execute",
        action="store_true",
        default=False,
        help="Execute tools if available (default: False, only generate scripts)"
    )
    
    args = parser.parse_args()
    
    # Progress callback for CLI
    def progress_callback(message: str, percent: float):
        print(f"[{percent:.1f}%] {message}")
    
    # Generate design plan
    plan = generate_design_plan(
        step2_json_path=args.step2_json,
        output_base_dir=args.output_base_dir,
        binder_length=args.binder_length,
        num_backbones=args.num_backbones,
        num_sequences_per_backbone=args.num_sequences_per_backbone,
        num_af_models=args.num_af_models,
        num_af_recycles=args.num_af_recycles,
        use_colabfold=args.use_colabfold,
        rfdiffusion_path=args.rfdiffusion_path,
        proteinmpnn_path=args.proteinmpnn_path,
        colabfold_path=args.colabfold_path,
        alphafold_path=args.alphafold_path,
        execute=args.execute,
        progress_callback=progress_callback if args.execute else None,
    )
    
    # Ensure output directory exists
    os.makedirs(os.path.dirname(os.path.abspath(args.out_plan)), exist_ok=True)
    
    # Write design plan
    with open(args.out_plan, "w") as f:
        json.dump(plan, f, indent=2)
    
    if plan.get("success"):
        print(f"✓ Generated design plan for {plan['num_designs']} design(s)")
        print(f"✓ Wrote design plan to: {args.out_plan}")
        print(f"✓ Output base directory: {args.output_base_dir}")
        print(f"✓ Generated {len(plan['generated_files'])} files")
        
        # Show execution results if execution was attempted
        if plan.get("execution_mode"):
            print("\n" + "="*60)
            print("EXECUTION RESULTS")
            print("="*60)
            
            tool_status = plan.get("tool_status", {})
            print("\nTool Availability:")
            for tool_name, status in tool_status.items():
                if status.get("available"):
                    print(f"  ✓ {tool_name}: Available ({status.get('command')})")
                else:
                    print(f"  ✗ {tool_name}: {status.get('error', 'Not available')}")
            
            execution_results = plan.get("execution_results", {})
            if execution_results:
                print("\nExecution Results:")
                for design_key, results in execution_results.items():
                    print(f"\n  {design_key}:")
                    for step, result in results.items():
                        if result.get("success"):
                            files_count = len(result.get("generated_files", []))
                            print(f"    ✓ {step}: {files_count} files generated")
                        else:
                            error = result.get("error", "Unknown error")
                            print(f"    ✗ {step}: {error}")
        
        print("\nDesign summary:")
        for design in plan["designs"]:
            print(f"  - {design['target_name']} / {design['epitope_id']}")
            print(f"    Output: {design['output_directory']}")
            print(f"    Total sequences: {design['total_sequences']}")
            print(f"    Total AF predictions: {design['total_af_predictions']}")
    else:
        print(f"✗ Error: {plan.get('error', 'Unknown error')}")
        sys.exit(1)


if __name__ == "__main__":
    import sys
    run_cli()


