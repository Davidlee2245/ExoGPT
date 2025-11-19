"""
Flask backend API for Exosome-GPT GUI
Executes Python modules and returns results
"""

import os
import sys
import json
import subprocess
from pathlib import Path
from flask import Flask, request, jsonify
from flask_cors import CORS

app = Flask(__name__)
CORS(app)  # Enable CORS for React frontend

# Add the project root to Python path
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


@app.route("/api/health", methods=["GET"])
def health():
    """Health check endpoint"""
    return jsonify({"status": "ok"})


@app.route("/api/step1/check", methods=["POST"])
def check_step1_data():
    """Check if local data exists for disease/biofluid"""
    try:
        data = request.json
        disease = data.get("disease", "")
        biofluid = data.get("biofluid", "")
        
        if not disease or not biofluid:
            return jsonify({"error": "disease and biofluid are required"}), 400
        
        from exo_gpt.data_downloader import DataDownloader
        
        downloader = DataDownloader()
        has_data = downloader.check_local_data(disease, biofluid)
        
        return jsonify({
            "success": True,
            "has_local_data": has_data
        })
        
    except Exception as e:
        import traceback
        return jsonify({
            "success": False,
            "error": str(e),
            "traceback": traceback.format_exc()
        }), 500


@app.route("/api/step1/download", methods=["POST"])
def download_step1_data():
    """Download datasets for disease/biofluid"""
    try:
        data = request.json
        disease = data.get("disease", "")
        biofluid = data.get("biofluid", "")
        
        if not disease or not biofluid:
            return jsonify({"error": "disease and biofluid are required"}), 400
        
        from exo_gpt.data_downloader import DataDownloader
        
        progress_messages = []
        
        def progress_callback(message: str, percent: float):
            progress_messages.append({"message": message, "percent": percent})
        
        downloader = DataDownloader(progress_callback=progress_callback)
        downloaded_files = downloader.download_all_sources(disease, biofluid)
        
        return jsonify({
            "success": True,
            "downloaded_files": downloaded_files,
            "progress": progress_messages
        })
        
    except Exception as e:
        import traceback
        return jsonify({
            "success": False,
            "error": str(e),
            "traceback": traceback.format_exc()
        }), 500


@app.route("/api/step1/run", methods=["POST"])
def run_step1():
    """Execute Step 1: EV Biomarker Finder"""
    try:
        data = request.json
        disease = data.get("disease", "")
        biofluid = data.get("biofluid", "")
        auto_download = data.get("auto_download", True)  # Default to True
        
        if not disease or not biofluid:
            return jsonify({"error": "disease and biofluid are required"}), 400
        
        # Check if data exists, download if needed
        from exo_gpt.data_downloader import download_data_if_needed
        
        progress_log = []
        
        def progress_callback(message: str, percent: float):
            progress_log.append({"message": message, "percent": percent})
        
        if auto_download:
            has_data, downloaded = download_data_if_needed(
                disease, biofluid, progress_callback
            )
        else:
            has_data = True
            downloaded = []
        
        # Import and run the module directly (faster than subprocess)
        from exo_gpt.step1_ev_biomarker_finder import aggregate_biomarkers, biomarkers_to_markdown
        
        # Run the aggregation
        result = aggregate_biomarkers(disease, biofluid)
        
        # Generate markdown table
        markdown_table = biomarkers_to_markdown(result["ev_biomarkers"], top_n=50)
        
        # Clean NaN values for JSON serialization
        import pandas as pd
        import numpy as np
        
        def clean_for_json(obj):
            if isinstance(obj, dict):
                return {k: clean_for_json(v) for k, v in obj.items()}
            elif isinstance(obj, list):
                return [clean_for_json(item) for item in obj]
            elif isinstance(obj, float):
                if pd.isna(obj) or np.isnan(obj):
                    return None
            elif isinstance(obj, type(pd.NA)) if hasattr(pd, 'NA') else False:
                return None
            return obj
        
        cleaned_result = clean_for_json(result)
        
        return jsonify({
            "success": True,
            "result": cleaned_result,
            "markdown_table": markdown_table,
            "download_progress": progress_log,
            "downloaded_files": downloaded if auto_download else []
        })
        
    except Exception as e:
        import traceback
        return jsonify({
            "success": False,
            "error": str(e),
            "traceback": traceback.format_exc()
        }), 500


if __name__ == "__main__":
    app.run(host="127.0.0.1", port=5000, debug=True)

