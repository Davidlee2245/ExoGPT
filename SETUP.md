# ExoGPT Setup Guide

Complete setup instructions for the ExoGPT environment with RFdiffusion, ProteinMPNN, and ColabFold support.

## 📋 Key Package Versions

| Package | Version | Installation Method |
|---------|---------|---------------------|
| **Python** | 3.10 | `conda create -n ExoGPT python=3.10` |
| **PyTorch** | 2.1.0+cu118 | `pip install torch==2.1.0 --index-url https://download.pytorch.org/whl/cu118` |
| **CUDA Toolkit** | 11.6 | `conda install cudatoolkit=11.6 -c conda-forge` |
| **DGL** | 0.9.1 (CUDA 11.6) | `pip install /path/to/dgl_cu116-0.9.1.post1-cp310-cp310-manylinux1_x86_64.whl` |
| **NumPy** | <2.0 | `pip install "numpy<2.0"` (required for PyTorch 2.1.0) |
| **omegaconf** | >=2.0.0 | `pip install omegaconf>=2.0.0` |
| **hydra-core** | >=1.0.0 | `pip install hydra-core>=1.0.0` |
| **e3nn** | >=0.5.0 | `pip install e3nn>=0.5.0` |
| **wandb** | >=0.15.0 | `pip install wandb>=0.15.0` |
| **pyrsistent** | >=0.18.0 | `pip install pyrsistent>=0.18.0` |

**Note**: PyTorch and CUDA toolkit must be installed separately (see instructions below). Other packages can be installed via `pip install -r requirements.txt`.

## System Requirements

- **OS**: Linux (tested on Ubuntu)
- **GPU**: NVIDIA GPU with CUDA support (RTX 3060 Ti or better recommended)
- **CUDA**: System CUDA 12.2+ (we use CUDA 11.6 in conda environment for DGL compatibility)
- **Python**: 3.10 (via conda)
- **Conda**: Anaconda or Miniconda

## Quick Start

### 1. Create and Activate Conda Environment

```bash
conda create -n ExoGPT python=3.10 -y
conda activate ExoGPT
```

### 2. Install Core Dependencies

```bash
# Install from requirements.txt
pip install -r requirements.txt
```

### 3. Install PyTorch with CUDA Support

**Important**: We use PyTorch 2.1.0 with CUDA 11.8 (backward compatible with CUDA 11.6 toolkit for DGL).

```bash
conda activate ExoGPT

# Uninstall any existing PyTorch
conda uninstall pytorch torchvision torchaudio pytorch-cuda -y

# Install PyTorch 2.1.0 with CUDA 11.8 (works with CUDA 11.6 toolkit)
pip install torch==2.1.0 torchvision==0.16.0 torchaudio==2.1.0 --index-url https://download.pytorch.org/whl/cu118

# Install CUDA 11.6 toolkit (needed for DGL cu116)
conda install cudatoolkit=11.6 -c conda-forge -y

# Fix numpy compatibility (PyTorch 2.1.0 requires numpy < 2.0)
pip install "numpy<2.0"
```

### 4. Install RFdiffusion

```bash
conda activate ExoGPT

# Install RFdiffusion dependencies
pip install omegaconf>=2.0.0 hydra-core>=1.0.0 e3nn>=0.5.0 wandb>=0.15.0
pip install pyrsistent>=0.18.0

# Install DGL with CUDA 11.6 support from wheel file (fixes CUDA device mismatch errors)
# First, download the wheel file:
# wget https://data.dgl.ai/wheels/repo.html/dgl_cu116-0.9.1.post1-cp310-cp310-manylinux1_x86_64.whl
# Or download from: https://github.com/dmlc/dgl/releases
pip uninstall -y dgl dgl-cu11* 2>/dev/null || true
pip install dgl_cu116-0.9.1.post1-cp310-cp310-manylinux1_x86_64.whl
# Or if you downloaded to a different location:
# pip install /path/to/dgl_cu116-0.9.1.post1-cp310-cp310-manylinux1_x86_64.whl

# Install SE3Transformer
cd RFdiffusion/env/SE3Transformer
pip install -e .
cd ../..

# Install RFdiffusion
cd RFdiffusion
pip install -e . --no-deps
cd ..
```

### 5. Download RFdiffusion Model Weights

```bash
cd RFdiffusion
mkdir -p models
cd models

# Download base checkpoint (required)
wget http://files.ipd.uw.edu/pub/RFdiffusion/20230627_172345_base_ckpt_epoch=489_val_loss=0.636.ckpt -O Base_ckpt.pt

# Download complex checkpoint (optional, for complex design)
wget http://files.ipd.uw.edu/pub/RFdiffusion/20230211_105516_complex_base_ckpt_epoch=489_val_loss=0.636.ckpt -O Complex_base_ckpt.pt

cd ../..
```

### 6. Download Structure Files

```bash
# Run the structure downloader script
python download_structures.py
```

This will download all PDB and AlphaFold structure files referenced in `data/structure_mapping.csv`.

### 7. Verify Installation

```bash
conda activate ExoGPT

# Test PyTorch
python -c "import torch; print(f'PyTorch: {torch.__version__}'); print(f'CUDA: {torch.cuda.is_available()}'); print(f'GPU: {torch.cuda.get_device_name(0) if torch.cuda.is_available() else \"N/A\"}')"

# Test DGL
python -c "import dgl; print(f'DGL: {dgl.__version__}')"

# Test RFdiffusion
python -c "import rfdiffusion; print('RFdiffusion: OK')"
```

Expected output:
```
PyTorch: 2.1.0+cu118
CUDA: True
GPU: NVIDIA GeForce RTX 3060 Ti
DGL: 1.1.3
RFdiffusion: OK
```


## Why CUDA 11.8 in Conda?

Even though your system has CUDA 12.2, we use CUDA 11.6 in the conda environment because:
- DGL 1.1.3 requires CUDA 11.6 for the cu116 build
- RFdiffusion and SE3Transformer work with CUDA 11.6
- PyTorch 2.1.0 with cu118 is backward compatible with CUDA 11.6 toolkit
- RTX 3060 Ti is fully compatible with CUDA 11.6
- Conda's cudatoolkit=11.6 provides the necessary libraries without conflicting with system CUDA

## Troubleshooting

### GPU Not Detected

If `torch.cuda.is_available()` returns `False`:

1. **Check NVIDIA driver**:
   ```bash
   nvidia-smi
   ```

2. **Verify PyTorch CUDA installation**:
   ```bash
   python -c "import torch; print(torch.version.cuda)"
   ```
   Should show `11.8` (PyTorch version, but toolkit is 11.6)

3. **Reinstall PyTorch**:
   ```bash
   pip uninstall torch torchvision torchaudio
   pip install torch==2.1.0 torchvision==0.16.0 torchaudio==2.1.0 --index-url https://download.pytorch.org/whl/cu118
   ```

### DGL Import Errors or CUDA Device Mismatch

If DGL fails to import or you get "Cannot assign edge feature on device cuda:0 to a graph on device cpu" errors:

1. **Check DGL version**:
   ```bash
   pip show dgl
   ```
   Should be `dgl-cu116==0.9.1.post1` (not the CPU-only version)

2. **Reinstall DGL with CUDA 11.6 support from wheel file**:
   ```bash
   pip uninstall -y dgl dgl-cu11*
   # Download wheel file first, then:
   pip install dgl_cu116-0.9.1.post1-cp310-cp310-manylinux1_x86_64.whl
   # Or use full path: pip install /path/to/dgl_cu116-0.9.1.post1-cp310-cp310-manylinux1_x86_64.whl
   ```

3. **Verify CUDA support**:
   ```bash
   python -c "import dgl; print(f'DGL: {dgl.__version__}'); print(f'CUDA backend: {dgl.backend.get_backend()}')"
   ```

4. **Check LD_LIBRARY_PATH** (automatically set by code):
   The code automatically sets `LD_LIBRARY_PATH` to include conda's lib directory.

### RFdiffusion Module Not Found

If `import rfdiffusion` fails:

1. **Verify installation**:
   ```bash
   cd RFdiffusion
   pip install -e . --no-deps
   ```

2. **Check PYTHONPATH** (automatically set by code):
   The code automatically adds RFdiffusion directory to PYTHONPATH.

### PDB File Not Found Errors

If you get "PDB file not found" errors:

1. **Download structure files**:
   ```bash
   python download_structures.py
   ```

2. **Check structure_mapping.csv**:
   Verify that all PDB IDs and UniProt IDs are correct in `data/structure_mapping.csv`.

### NumPy Compatibility Issues

If you get numpy-related errors:

```bash
pip install "numpy<2.0"
```

PyTorch 2.1.0 requires numpy < 2.0.

## Environment Configuration

The ExoGPT environment uses:
- **Single environment**: All tools (GUI, RFdiffusion, ProteinMPNN, ColabFold) run in the same conda environment
- **GPU support**: PyTorch uses GPU, DGL is forced to CPU mode (via `DGL_BACKEND=cpu`) to avoid CUDA compatibility issues
- **Automatic path configuration**: The code automatically sets PYTHONPATH, LD_LIBRARY_PATH, and DGL_BACKEND as needed

## Additional Tools (Optional)

### ProteinMPNN

ProteinMPNN is optional. If you want to install it:

```bash
cd ProteinMPNN
# Follow ProteinMPNN installation instructions
```

### ColabFold

ColabFold can be installed via:

```bash
pip install colabfold
# Or use the colabfold_batch command if installed system-wide
```

## Notes

- **SE3nv environment**: Not needed! Everything runs in ExoGPT environment.
- **CPU mode**: If GPU setup is problematic, the code will automatically fall back to CPU mode (slower but functional).
- **Model weights**: RFdiffusion model weights are large (~2GB each). Ensure sufficient disk space.

## Verification Checklist

Before running Step 3, verify:

- [ ] PyTorch 2.1.0+cu118 installed and CUDA available
- [ ] DGL 0.9.1 with CUDA 11.6 support installed from wheel file (dgl_cu116-0.9.1.post1)
- [ ] RFdiffusion imports successfully
- [ ] SE3Transformer installed
- [ ] RFdiffusion model weights downloaded
- [ ] Structure files downloaded (run `download_structures.py`)
- [ ] GPU detected (optional but recommended)

## Getting Help

If you encounter issues:

1. Check the error messages in the GUI progress log
2. Verify all package versions match the requirements
3. Ensure model weights are downloaded
4. Check that structure files exist in `structures/` directory
5. Review the troubleshooting section above

