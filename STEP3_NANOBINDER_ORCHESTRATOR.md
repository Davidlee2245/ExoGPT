# Step 3: Nanobinder Design Orchestrator

## Overview

Step 3 orchestrates the nanobinder design pipeline by generating reusable scripts and configuration files for:

1. **RFdiffusion** - Backbone generation
2. **ProteinMPNN** - Sequence design
3. **AlphaFold-Multimer / ColabFold** - Complex structure prediction

**Important**: Step 3 does NOT execute these tools. Instead, it generates scripts and configs that can be run on external machines/HPC clusters where these tools are installed.

## Input

Step 3 takes the JSON output from Step 2, which contains curated targets with epitope definitions. Example:

```json
{
  "disease": "Metastatic melanoma",
  "biofluid": "Plasma",
  "curated_targets": [
    {
      "target_name": "PD-L1",
      "gene_symbol": "CD274",
      "uniprot": "Q9NZQ7",
      "structure_type": "pdb",
      "model_path": "./structures/pd_l1_ecd.pdb",
      "pdb_id": "5JDR",
      "epitopes": [
        {
          "epitope_id": "cd274_epitope_01",
          "epitope_residues": [52, 54, 55, 76, 77, 120, 121],
          "epitope_description": "Surface-exposed loop patch on IgV domain"
        }
      ]
    }
  ]
}
```

## Usage

### Command Line

```bash
python -m exo_gpt.step3_nanobinder_orchestrator \
  --step2_json "./epitopes/mm_plasma_epitopes.json" \
  --output_base_dir "./workflows" \
  --binder_length 90 \
  --num_backbones 10 \
  --num_sequences_per_backbone 8 \
  --num_af_models 5 \
  --num_af_recycles 3 \
  --use_colabfold \
  --out_plan "./workflows/design_plan.json"
```

### Parameters

- `--step2_json`: Path to Step 2 JSON output (required)
- `--output_base_dir`: Base directory for all outputs (default: `./workflows`)
- `--binder_length`: Length of binder in amino acids (default: 90)
- `--num_backbones`: Number of RFdiffusion backbone designs (default: 10)
- `--num_sequences_per_backbone`: Number of ProteinMPNN sequences per backbone (default: 8)
- `--num_af_models`: Number of AlphaFold models per sequence (default: 5)
- `--num_af_recycles`: Number of AlphaFold recycles (default: 3)
- `--use_colabfold`: Use ColabFold (default: True). Use `--no-use_colabfold` for AlphaFold-Multimer
- `--rfdiffusion_path`: Optional path to RFdiffusion executable
- `--proteinmpnn_path`: Optional path to ProteinMPNN executable
- `--colabfold_path`: Optional path to ColabFold executable
- `--alphafold_path`: Optional path to AlphaFold-Multimer executable
- `--out_plan`: Path to write design plan JSON (required)

## Output Structure

For each target/epitope combination, Step 3 generates:

```
workflows/
└── designs/
    └── PD-L1/
        └── cd274_epitope_01/
            ├── config.yaml          # YAML configuration
            ├── config.json          # JSON configuration
            ├── run_rfdiffusion.sh   # RFdiffusion script
            ├── run_proteinmpnn.sh   # ProteinMPNN script
            ├── run_colabfold.sh     # ColabFold script (or run_alphafold.sh)
            ├── run_pipeline.sh      # Master script (runs all steps)
            └── [output directories created when scripts run]
                ├── rfdiffusion_output/
                ├── proteinmpnn_output/
                └── colabfold_output/
```

## Generated Scripts

### 1. RFdiffusion Script (`run_rfdiffusion.sh`)

Generates binder backbones using RFdiffusion with hotspot-guided design.

**Inputs:**
- Target PDB file
- Target chain ID
- Epitope residues (hotspot constraints)
- Binder length

**Outputs:**
- Complex PDBs (target + binder backbone) in `rfdiffusion_output/`

**Usage:**
```bash
cd workflows/designs/PD-L1/cd274_epitope_01/
bash run_rfdiffusion.sh
```

**Note**: Adjust the RFdiffusion command in the script based on your installation and version.

### 2. ProteinMPNN Script (`run_proteinmpnn.sh`)

Designs sequences for each RFdiffusion backbone.

**Inputs:**
- RFdiffusion output PDBs
- Binder chain ID
- Omitted amino acids (default: C, X)
- Number of sequences per backbone

**Outputs:**
- Designed sequences in FASTA format in `proteinmpnn_output/`

**Usage:**
```bash
bash run_proteinmpnn.sh
```

**Note**: Adjust the ProteinMPNN command based on your installation.

### 3. ColabFold Script (`run_colabfold.sh`)

Predicts complex structures for designed sequences.

**Inputs:**
- ProteinMPNN output sequences
- Target sequence (extracted from PDB)
- Number of models and recycles

**Outputs:**
- Predicted complex structures in `colabfold_output/`

**Usage:**
```bash
bash run_colabfold.sh
```

**Note**: 
- The script includes a placeholder for target sequence extraction. You'll need to implement sequence extraction from the PDB file.
- Adjust the ColabFold command based on your installation.

### 4. Master Pipeline Script (`run_pipeline.sh`)

Runs all three steps sequentially.

**Usage:**
```bash
bash run_pipeline.sh
```

## Configuration Files

Each design directory contains `config.yaml` and `config.json` with all parameters:

```yaml
target_name: PD-L1
gene_symbol: CD274
uniprot: Q9NZQ7
structure_type: pdb
model_path: ./structures/pd_l1_ecd.pdb
pdb_id: 5JDR
target_chain_id: A
epitope_id: cd274_epitope_01
epitope_residues: [52, 54, 55, 76, 77, 120, 121]
binder_length: 90
binder_style: helical_mini_domain
num_backbones: 10
num_sequences_per_backbone: 8
num_af_models: 5
num_af_recycles: 3
output_base_dir: ./workflows
use_colabfold: true
omitted_aa: [C, X]
binder_chain_id: B
```

## Example: PD-L1 Workflow

1. **Run Step 2** to get curated targets:
   ```bash
   python -m exo_gpt.step2_target_epitope_curator \
     --step1_json "./scores/mm_plasma_biomarkers.json" \
     --target "CD274" \
     --out_json "./epitopes/pd_l1_epitopes.json"
   ```

2. **Run Step 3** to generate design plan:
   ```bash
   python -m exo_gpt.step3_nanobinder_orchestrator \
     --step2_json "./epitopes/pd_l1_epitopes.json" \
     --output_base_dir "./workflows" \
     --out_plan "./workflows/design_plan.json"
   ```

3. **Review generated scripts** in `./workflows/designs/PD-L1/cd274_epitope_01/`

4. **Adjust scripts** for your RFdiffusion/ProteinMPNN/ColabFold installation

5. **Run on HPC**:
   ```bash
   cd ./workflows/designs/PD-L1/cd274_epitope_01/
   bash run_pipeline.sh
   ```

## Important Notes

1. **Tool Installation**: The generated scripts are templates. You must:
   - Install RFdiffusion, ProteinMPNN, and ColabFold/AlphaFold-Multimer
   - Adjust script commands to match your installation paths and API
   - Test each script individually before running the full pipeline

2. **Sequence Extraction**: The ColabFold script includes a placeholder for target sequence extraction. Implement this using:
   - BioPython: `Bio.PDB` and `Bio.SeqIO`
   - PyMOL: `pymol.cmd.get_fastastr()`
   - Or any other PDB parsing tool

3. **HPC Execution**: These scripts are designed for HPC/cluster execution:
   - Use job schedulers (SLURM, PBS, etc.) to submit jobs
   - Adjust resource requirements (CPU, GPU, memory) based on your needs
   - Consider parallelizing across multiple designs

4. **Output Organization**: The scripts organize outputs in structured directories. Adjust paths in the scripts if needed.

5. **Reproducibility**: All parameters are saved in config files. Use these to reproduce designs or adjust parameters.

## Design Plan Output

Step 3 generates a design plan JSON summarizing all designs:

```json
{
  "success": true,
  "disease": "Metastatic melanoma",
  "biofluid": "Plasma",
  "num_designs": 1,
  "output_base_dir": "./workflows",
  "designs": [
    {
      "target_name": "PD-L1",
      "gene_symbol": "CD274",
      "epitope_id": "cd274_epitope_01",
      "epitope_residues": [52, 54, 55, 76, 77, 120, 121],
      "output_directory": "./workflows/designs/PD-L1/cd274_epitope_01",
      "binder_length": 90,
      "num_backbones": 10,
      "num_sequences_per_backbone": 8,
      "total_sequences": 80,
      "num_af_models": 5,
      "total_af_predictions": 400
    }
  ],
  "generated_files": [...],
  "pipeline_steps": [...]
}
```

## Troubleshooting

1. **No designs generated**: Check that Step 2 output contains curated targets with epitopes.

2. **Scripts fail**: Verify tool installations and adjust script commands.

3. **Path issues**: Use absolute paths or adjust relative paths in scripts.

4. **Missing sequences**: Implement target sequence extraction in ColabFold script.

## Next Steps

After running Step 3 and executing the generated scripts:

- **Step 4**: Analyze predicted structures (binding affinity, stability, etc.)
- **Step 5**: Rank and select best designs for experimental validation


