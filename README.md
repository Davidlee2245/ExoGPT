# ExoGPT: Exosome-GPT Designer

An expert multi-agent system for end-to-end exosome biomarker discovery, nanobinder design, in-silico scoring, and wet-lab validation pipelines.

## Features

- **Step 1: EV Biomarker Finder** - Automatically searches across publications, databases (Vesiclepedia, ExoCarta, EVpedia, EVmiRNA, exoRBase), and user experiments. Normalizes disease and biofluid, aggregates EV datasets, and ranks exosomal biomarkers with a bias toward surface proteins.
- **Step 2: Target & Epitope Curator** - Curates target epitopes for nanobinder design
- **Step 3: Nanobinder Orchestrator** - Designs nanobinders using RFdiffusion, ProteinMPNN, and AlphaFold-Multimer
- **Step 4: Scoring & Off-target** - In-silico scoring using ipTM, ipAE, pLDDT, and off-target filtering
- **Step 5: Wet-lab Protocol** - Generates wet-lab protocols for expression, purification, binding assays, and exosome capture

## Quick Start

### Prerequisites

- Python 3.8+
- Node.js 18+
- Conda (recommended)

### Installation

1. Clone the repository:
```bash
git clone https://github.com/Davidlee2245/ExoGPT.git
cd ExoGPT
```

2. Set up Python environment:
```bash
conda create -n ExoGPT python=3.10
conda activate ExoGPT
pip install -r requirements.txt
```

3. Set up frontend:
```bash
cd gui
npm install
```

### Running the Application

1. Start the backend (in one terminal):
```bash
conda activate ExoGPT
cd backend
python app.py
```

2. Start the frontend (in another terminal):
```bash
conda activate ExoGPT
cd gui
npm run dev
```

3. Open your browser to `http://localhost:5173`

## Project Structure

```
ExoGPT/
├── exo_gpt/              # Python package with pipeline steps
│   ├── step1_ev_biomarker_finder.py
│   ├── data_downloader.py
│   └── ...
├── backend/              # Flask API server
│   ├── app.py
│   └── requirements.txt
├── gui/                  # React frontend (Gemini-style UI)
│   ├── src/
│   └── ...
├── data/                 # EV biomarker datasets
│   ├── publications/
│   ├── databases/
│   └── experiments/
└── scores/              # Output files
```

## Data System

The system automatically searches for EV biomarker data in:
- `./data/publications/` - Extracted from published articles
- `./data/databases/` - Pre-packaged from Vesiclepedia, ExoCarta, EVpedia, EVmiRNA, exoRBase
- `./data/experiments/` - User-generated experimental data

If data is not found locally, the system automatically downloads datasets from EV databases with progress tracking.

## Documentation

- [Quick Start Guide](QUICK_START_GUI.md)
- [Data System Documentation](DATA_SYSTEM.md)
- [Auto-Download Feature](AUTO_DOWNLOAD_FEATURE.md)

## License

[Add your license here]

## Citation

[Add citation information here]

