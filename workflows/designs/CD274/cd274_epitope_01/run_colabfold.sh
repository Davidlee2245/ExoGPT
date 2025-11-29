#!/bin/bash
# ColabFold Complex Prediction Script
# Target: CD274 (CD274)
# Epitope: cd274_epitope_01

set -e  # Exit on error

# Configuration
PROTEINMPNN_OUTPUT="/home/david/.cursor-tutor/ExoGPT/./workflows/designs/CD274/cd274_epitope_01/proteinmpnn_output"
COLABFOLD_OUTPUT="/home/david/.cursor-tutor/ExoGPT/./workflows/designs/CD274/cd274_epitope_01/colabfold_output"
TARGET_PDB="./structures/pd_l1_ecd.pdb"
TARGET_CHAIN="A"
BINDER_CHAIN="B"
NUM_MODELS=1
NUM_RECYCLES=1

# Create output directory
mkdir -p "$COLABFOLD_OUTPUT"

# Extract target sequence from PDB
# NOTE: In practice, use a tool like BioPython or PyMOL to extract sequence
# For now, this is a placeholder - you'll need to implement sequence extraction
TARGET_SEQ="TARGET_SEQUENCE_FROM_PDB"  # TODO: Extract from $TARGET_PDB chain $TARGET_CHAIN

# ColabFold command
# NOTE: Adjust the command below based on your ColabFold installation
if [ -z "$COLABFOLD_PATH" ]; then
    COLABFOLD_CMD="colabfold_batch"
else
    COLABFOLD_CMD="$COLABFOLD_PATH"
fi

# Process each ProteinMPNN output
# ProteinMPNN outputs to {output_dir}/{backbone_name}/seqs/{backbone_name}.fa
for backbone_dir in "$PROTEINMPNN_OUTPUT"/*/; do
    if [ ! -d "$backbone_dir" ]; then
        echo "No ProteinMPNN outputs found in $PROTEINMPNN_OUTPUT"
        exit 1
    fi

    # Find .fa file in seqs subdirectory
    seq_file=$(find "$backbone_dir/seqs" -name "*.fa" | head -1)
    if [ ! -f "$seq_file" ]; then
        echo "Warning: No .fa file found in $backbone_dir/seqs"
        continue
    fi

    # Extract base name from directory (e.g., backbone_1_0/seqs/backbone_1_0.fa -> backbone_1_0)
    BASE_NAME=$(basename "$backbone_dir")

    # Create batch FASTA file with target + binder pairs
    # Format: For each binder sequence, create a FASTA entry with target and binder
    BATCH_FASTA="$COLABFOLD_OUTPUT/${BASE_NAME}_batch.fasta"
    > "$BATCH_FASTA"  # Clear/create file

    # Read binder sequences from ProteinMPNN output
    # ProteinMPNN typically outputs sequences in FASTA format
    # Format: >seq_0\nSEQUENCE\n>seq_1\nSEQUENCE...
    seq_count=0
    while IFS= read -r line; do
        if [[ $line =~ ^\> ]]; then
            # This is a header line, skip or use as identifier
            continue
        elif [[ -n "$line" && ! $line =~ ^# ]]; then
            # This is a sequence line
            BINDER_SEQ="$line"
            seq_count=$((seq_count + 1))

            # Write target:binder pair to batch FASTA
            echo ">${BASE_NAME}_seq${seq_count}_target" >> "$BATCH_FASTA"
            echo "$TARGET_SEQ" >> "$BATCH_FASTA"
            echo ">${BASE_NAME}_seq${seq_count}_binder" >> "$BATCH_FASTA"
            echo "$BINDER_SEQ" >> "$BATCH_FASTA"
        fi
    done < "$seq_file"

    if [ $seq_count -eq 0 ]; then
        echo "Warning: No sequences found in $seq_file"
        continue
    fi

    # Run ColabFold batch prediction
    OUTPUT_DIR_DESIGN="$COLABFOLD_OUTPUT/${BASE_NAME}"
    mkdir -p "$OUTPUT_DIR_DESIGN"

    $COLABFOLD_CMD \
        "$BATCH_FASTA" \
        "$OUTPUT_DIR_DESIGN" \
        --num-models $NUM_MODELS \
        --num-recycle $NUM_RECYCLES \
        --model-type alphafold2_multimer_v3

    if [ $? -eq 0 ]; then
        echo "✓ Predicted structures for $BASE_NAME ($seq_count sequences)"
    else
        echo "✗ Failed to predict structures for $BASE_NAME"
        exit 1
    fi
done

echo "
ColabFold completed: structures predicted for all designs
Output directory: $COLABFOLD_OUTPUT
"