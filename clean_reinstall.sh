#!/bin/bash
# Clean Reinstall Script for ExoGPT
# This script removes the ExoGPT conda environment and all packages, then reinstalls everything

set -e  # Exit on error

echo "=========================================="
echo "ExoGPT Clean Reinstall Script"
echo "=========================================="
echo ""
echo "This will:"
echo "  1. Remove the ExoGPT conda environment (if it exists)"
echo "  2. Create a fresh ExoGPT environment with Python 3.11"
echo "  3. Install all dependencies from requirements.txt"
echo "  4. Install PyTorch with CUDA 11.8 support"
echo "  5. Install CUDA 11.6 toolkit (for DGL)"
echo "  6. Install DGL with CUDA 11.6 support"
echo "  7. Install RFdiffusion and SE3Transformer"
echo ""
read -p "Continue? (y/n) " -n 1 -r
echo ""
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Aborted."
    exit 1
fi

echo ""
echo "Step 1: Removing existing ExoGPT environment..."
conda deactivate 2>/dev/null || true
conda env remove -n ExoGPT -y 2>/dev/null || echo "  (Environment doesn't exist, skipping)"

echo ""
echo "Step 2: Creating fresh ExoGPT environment with Python 3.10..."
conda create -n ExoGPT python=3.10 -y

echo ""
echo "Step 3: Activating ExoGPT environment..."
source $(conda info --base)/etc/profile.d/conda.sh
conda activate ExoGPT

echo ""
echo "Step 4: Installing core dependencies from requirements.txt..."
pip install -r requirements.txt

echo ""
echo "Step 5: Installing PyTorch 2.1.0 with CUDA 11.8 support..."
pip install torch==2.1.0 torchvision==0.16.0 torchaudio==2.1.0 --index-url https://download.pytorch.org/whl/cu118

echo ""
echo "Step 6: Installing CUDA 11.6 toolkit (for DGL compatibility)..."
conda install cudatoolkit=11.6 -c conda-forge -y

echo ""
echo "Step 7: Fixing numpy compatibility..."
pip install "numpy<2.0"

echo ""
echo "Step 8: Installing DGL with CUDA 11.6 support from wheel file..."
pip uninstall -y dgl dgl-cu11* 2>/dev/null || true

# Check if wheel file exists in current directory
if [ -f "dgl_cu116-0.9.1.post1-cp310-cp310-manylinux1_x86_64.whl" ]; then
    echo "  Found wheel file in current directory, installing..."
    pip install dgl_cu116-0.9.1.post1-cp310-cp310-manylinux1_x86_64.whl
elif [ -f "wheels/dgl_cu116-0.9.1.post1-cp310-cp310-manylinux1_x86_64.whl" ]; then
    echo "  Found wheel file in wheels/ directory, installing..."
    pip install wheels/dgl_cu116-0.9.1.post1-cp310-cp310-manylinux1_x86_64.whl
else
    echo "  WARNING: DGL wheel file not found!"
    echo "  Please download dgl_cu116-0.9.1.post1-cp310-cp310-manylinux1_x86_64.whl"
    echo "  From: https://data.dgl.ai/wheels/repo.html or https://github.com/dmlc/dgl/releases"
    echo "  Then run: pip install dgl_cu116-0.9.1.post1-cp310-cp310-manylinux1_x86_64.whl"
    read -p "  Press Enter to continue without DGL, or Ctrl+C to abort and download the wheel file first..."
fi

echo ""
echo "Step 9: Installing RFdiffusion dependencies..."
pip install omegaconf>=2.0.0 hydra-core>=1.0.0 e3nn>=0.5.0 wandb>=0.15.0 pyrsistent>=0.18.0

echo ""
echo "Step 10: Installing SE3Transformer..."
if [ -d "RFdiffusion/env/SE3Transformer" ]; then
    cd RFdiffusion/env/SE3Transformer
    pip install -e .
    cd ../../..
else
    echo "  WARNING: RFdiffusion/env/SE3Transformer not found. Skipping SE3Transformer installation."
    echo "  You may need to clone RFdiffusion first."
fi

echo ""
echo "Step 11: Installing RFdiffusion..."
if [ -d "RFdiffusion" ]; then
    cd RFdiffusion
    pip install -e . --no-deps
    cd ..
else
    echo "  WARNING: RFdiffusion directory not found. Skipping RFdiffusion installation."
    echo "  You may need to clone RFdiffusion first."
fi

echo ""
echo "=========================================="
echo "Installation Complete!"
echo "=========================================="
echo ""
echo "Verifying installation..."
echo ""

echo "Python version:"
python --version

echo ""
echo "PyTorch version and CUDA:"
python -c "import torch; print(f'PyTorch: {torch.__version__}'); print(f'CUDA available: {torch.cuda.is_available()}'); print(f'CUDA version: {torch.version.cuda if torch.cuda.is_available() else \"N/A\"}')"

echo ""
echo "DGL version:"
python -c "import dgl; print(f'DGL: {dgl.__version__}')" || echo "  ERROR: DGL not installed correctly"

echo ""
echo "NumPy version:"
python -c "import numpy; print(f'NumPy: {numpy.__version__}')"

echo ""
echo "=========================================="
echo "Next steps:"
echo "=========================================="
echo "1. Download RFdiffusion model weights:"
echo "   cd RFdiffusion/models && wget http://files.ipd.uw.edu/pub/RFdiffusion/20230627_172345_base_ckpt_epoch=489_val_loss=0.636.ckpt -O Base_ckpt.pt"
echo ""
echo "2. Download structure files:"
echo "   python download_structures.py"
echo ""
echo "3. Activate environment:"
echo "   conda activate ExoGPT"
echo ""

