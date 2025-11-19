# Quick Start - Exosome-GPT GUI

## 🚀 Start the GUI (2 terminals)

### Terminal 1 - Backend Server:
```bash
conda activate ExoGPT
cd /home/david/.cursor-tutor/ExoGPT/backend

# Install Flask if needed
pip install flask flask-cors

# Start backend
python app.py
```

Backend will run on **http://127.0.0.1:5000**

### Terminal 2 - Frontend (React):
```bash
conda activate ExoGPT
cd /home/david/.cursor-tutor/ExoGPT/gui

# Start dev server
npm run dev
```

Frontend will run on **http://localhost:5173**

## 🎨 Open in Browser

Open **http://localhost:5173** to see the Gemini-style interface!

## ✨ Features

- **Sidebar Navigation** - Click any step to switch
- **Step 1** - EV Biomarker Finder with:
  - Input fields for disease and biofluid
  - **Run button** to execute directly
  - Results table showing biomarkers from all data sources
  - Evidence levels, scores, and supporting studies
- **Steps 2-5** - Input forms and command generation
- **Dark Theme** - Modern Gemini-style UI

## 📊 Step 1 Example

1. Enter: **Disease**: "Metastatic melanoma"
2. Enter: **Biofluid**: "Plasma"
3. Click **"Run Step 1"**
4. View results:
   - 55+ biomarkers found
   - From 7 data sources (publications, databases)
   - Ranked by evidence level
   - CD274 (PD-L1) typically #1 with high evidence

## 🔧 Troubleshooting

**If backend not found:**
```bash
pip install flask flask-cors
```

**If frontend errors:**
```bash
cd gui
rm -rf node_modules package-lock.json
npm install
```

**If Node.js version error:**
```bash
conda activate ExoGPT
conda install -y -c conda-forge nodejs=18
cd gui
rm -rf node_modules package-lock.json
npm install
```

## 📝 Notes

- Backend must be running for Step 1 "Run" button to work
- Steps 2-5 show commands that can be run manually
- All data is automatically loaded from `./data/` subdirectories

