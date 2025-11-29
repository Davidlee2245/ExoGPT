#!/bin/bash
# Master Pipeline Script
# Target: CD274 (CD274)
# Epitope: cd274_epitope_02

set -e  # Exit on error

echo "========================================"
echo "Nanobinder Design Pipeline"
echo "Target: CD274"
echo "Epitope: cd274_epitope_02"
echo "========================================"

# Step 1: RFdiffusion backbone generation
echo "
Step 1/3: RFdiffusion Backbone Generation
"
bash "/home/david/.cursor-tutor/ExoGPT/./workflows/designs/CD274/cd274_epitope_02/run_rfdiffusion.sh"

# Step 2: ProteinMPNN sequence design
echo "
Step 2/3: ProteinMPNN Sequence Design
"
bash "/home/david/.cursor-tutor/ExoGPT/./workflows/designs/CD274/cd274_epitope_02/run_proteinmpnn.sh"

# Step 3: Structure prediction
echo "
Step 3/3: Complex Structure Prediction
"
bash "/home/david/.cursor-tutor/ExoGPT/./workflows/designs/CD274/cd274_epitope_02/run_colabfold.sh"

echo "
==========================================
Pipeline completed successfully!
==========================================
"