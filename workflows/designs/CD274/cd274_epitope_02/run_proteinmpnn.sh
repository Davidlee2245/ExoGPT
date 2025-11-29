#!/bin/bash
# ProteinMPNN Sequence Design Script
# Target: CD274 (CD274)
# Epitope: cd274_epitope_02

set -e  # Exit on error

# Configuration
RFDIFFUSION_OUTPUT="/home/david/.cursor-tutor/ExoGPT/./workflows/designs/CD274/cd274_epitope_02/rfdiffusion_output"
PROTEINMPNN_OUTPUT="/home/david/.cursor-tutor/ExoGPT/./workflows/designs/CD274/cd274_epitope_02/proteinmpnn_output"
BINDER_CHAIN="B"
OMITTED_AA="C,X"
NUM_SEQUENCES=8

# Create output directory
mkdir -p "$PROTEINMPNN_OUTPUT"

# ProteinMPNN command
# NOTE: Adjust the command below based on your ProteinMPNN installation
if [ -z "$PROTEINMPNN_PATH" ]; then
    PROTEINMPNN_CMD="python protein_mpnn_run.py"
else
    PROTEINMPNN_CMD="$PROTEINMPNN_PATH"
fi

# Process each RFdiffusion output
# RFdiffusion generates files with pattern: backbone_*.pdb (e.g., backbone_1_0.pdb, backbone_2_0.pdb)
for complex_pdb in "$RFDIFFUSION_OUTPUT"/backbone_*.pdb; do
    if [ ! -f "$complex_pdb" ]; then
        echo "No RFdiffusion outputs found in $RFDIFFUSION_OUTPUT"
        exit 1
    fi

    # Extract base name (e.g., backbone_1_0.pdb -> backbone_1_0)
    BASE_NAME=$(basename "$complex_pdb" .pdb)
    BACKBONE_OUTPUT_DIR="$PROTEINMPNN_OUTPUT/${BASE_NAME}"
    mkdir -p "$BACKBONE_OUTPUT_DIR"

    # Run ProteinMPNN for this backbone
    # Format: ProteinMPNN with specified chain and omitted amino acids
    # Note: ProteinMPNN uses --pdb_path_chains (not --chain_id) and --out_folder (not --out_path)
    # Output will be in $BACKBONE_OUTPUT_DIR/seqs/${BASE_NAME}.fa
    cd "$(dirname "$PROTEINMPNN_CMD")" || cd ProteinMPNN || cd .
    $PROTEINMPNN_CMD \
        --pdb_path "$complex_pdb" \
        --pdb_path_chains "$BINDER_CHAIN" \
        --num_seq_per_target $NUM_SEQUENCES \
        --omit_AAs "$OMITTED_AA" \
        --out_folder "$BACKBONE_OUTPUT_DIR"

    if [ $? -eq 0 ]; then
        echo "✓ Designed sequences for $BASE_NAME"
    else
        echo "✗ Failed to design sequences for $BASE_NAME"
        exit 1
    fi
done

echo "
ProteinMPNN completed: sequences designed for all backbones
Output directory: $PROTEINMPNN_OUTPUT
"