# Exosome-GPT GUI

Gemini-style React GUI for the Exosome-GPT nanobinder design pipeline.

## Quick Start

### Prerequisites
- Node.js >= 18 (installed in conda environment)
- Backend server running (see below)

### Setup

```bash
# Activate conda environment
conda activate ExoGPT

# Install dependencies (if not already installed)
cd gui
npm install
```

### Run

**Terminal 1 - Backend:**
```bash
conda activate ExoGPT
cd /home/david/.cursor-tutor/ExoGPT/backend
pip install flask flask-cors  # if not already installed
python app.py
```

**Terminal 2 - Frontend:**
```bash
conda activate ExoGPT
cd /home/david/.cursor-tutor/ExoGPT/gui
npm run dev
```

Then open **http://localhost:5173** in your browser.

## Features

- **Gemini-style sidebar** with 5 pipeline steps
- **Step 1**: EV Biomarker Finder with Run button and results table
- **Step 2-5**: Input forms and command generation
- **Dark theme** with modern UI
- **Responsive design**

## Build for Production

```bash
npm run build
```

Output will be in `dist/` directory.

