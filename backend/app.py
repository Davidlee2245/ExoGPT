"""
Flask backend API for Exosome-GPT GUI
Executes Python modules and returns results
"""

import os
import sys
import json
import subprocess
import threading
import queue
from pathlib import Path
from flask import Flask, request, jsonify, Response, stream_with_context
from flask_cors import CORS

app = Flask(__name__)
# Enable CORS for React frontend, including SSE endpoints
CORS(app, resources={
    r"/api/*": {
        "origins": "*",
        "methods": ["GET", "POST", "OPTIONS"],
        "allow_headers": ["Content-Type", "Authorization"]
    }
})

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
        
        # Optionally save results to file for Step 2 to use
        scores_dir = os.path.join(PROJECT_ROOT, "scores")
        os.makedirs(scores_dir, exist_ok=True)
        
        # Generate filename from disease and biofluid
        safe_disease = disease.lower().replace(" ", "_").replace("/", "_")
        safe_biofluid = biofluid.lower().replace(" ", "_")
        json_filename = f"{safe_disease}_{safe_biofluid}_biomarkers.json"
        json_path = os.path.join(scores_dir, json_filename)
        
        try:
            with open(json_path, "w") as f:
                json.dump(cleaned_result, f, indent=2)
        except Exception as e:
            # Don't fail if file write fails, just log it
            print(f"Warning: Could not save Step 1 results to {json_path}: {e}")
        
        return jsonify({
            "success": True,
            "result": cleaned_result,
            "markdown_table": markdown_table,
            "download_progress": progress_log,
            "downloaded_files": downloaded if auto_download else [],
            "saved_json_path": json_path if os.path.exists(json_path) else None
        })
        
    except Exception as e:
        import traceback
        return jsonify({
            "success": False,
            "error": str(e),
            "traceback": traceback.format_exc()
        }), 500


@app.route("/api/step2/run", methods=["POST"])
def run_step2():
    """Execute Step 2: Target & Epitope Curator"""
    try:
        data = request.json
        step1_json_path = data.get("step1_json_path")
        step1_json_data = data.get("step1_json_data")  # Alternative: pass JSON directly
        target = data.get("target")  # Optional specific target
        max_targets = data.get("max_targets", 5)
        structure_mapping_path = data.get(
            "structure_mapping_path",
            os.path.join(PROJECT_ROOT, "data", "structure_mapping.csv")
        )
        
        # Either step1_json_path or step1_json_data must be provided
        if not step1_json_path and not step1_json_data:
            return jsonify({"error": "step1_json_path or step1_json_data is required"}), 400
        
        # Import Step 2 module
        from exo_gpt.step2_target_epitope_curator import (
            curate_targets,
            curated_targets_to_markdown
        )
        
        # If JSON data is provided directly, write it to a temp file
        if step1_json_data:
            import tempfile
            with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
                json.dump(step1_json_data, f)
                step1_json_path = f.name
        else:
            # Ensure step1_json_path is absolute and normalize path
            if step1_json_path:
                # Remove leading ./ if present
                if step1_json_path.startswith('./'):
                    step1_json_path = step1_json_path[2:]
                # Make absolute if relative
                if not os.path.isabs(step1_json_path):
                    step1_json_path = os.path.join(PROJECT_ROOT, step1_json_path)
                # Normalize path (remove redundant separators, resolve .., etc.)
                step1_json_path = os.path.normpath(step1_json_path)
        
        if not os.path.exists(step1_json_path):
            return jsonify({"error": f"Step 1 JSON file not found: {step1_json_path}"}), 400
        
        # Ensure structure_mapping_path is absolute
        if not os.path.isabs(structure_mapping_path):
            structure_mapping_path = os.path.join(PROJECT_ROOT, structure_mapping_path)
        
        # Run curation
        result = curate_targets(
            step1_json_path=step1_json_path,
            structure_mapping_path=structure_mapping_path,
            user_specified_target=target,
            max_targets=max_targets
        )
        
        # Generate markdown table
        markdown_table = curated_targets_to_markdown(result["curated_targets"])
        
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
        
        # Clean up temp file if created
        if step1_json_data and os.path.exists(step1_json_path):
            try:
                os.unlink(step1_json_path)
            except:
                pass
        
        return jsonify({
            "success": True,
            "result": cleaned_result,
            "markdown_table": markdown_table
        })
        
    except Exception as e:
        import traceback
        return jsonify({
            "success": False,
            "error": str(e),
            "traceback": traceback.format_exc()
        }), 500


@app.route("/api/step3/stream", methods=["POST", "OPTIONS"])
def run_step3_stream():
    """Execute Step 3 with Server-Sent Events (SSE) for real-time progress streaming"""
    # Handle CORS preflight
    if request.method == "OPTIONS":
        response = jsonify({})
        response.headers.add("Access-Control-Allow-Origin", "*")
        response.headers.add("Access-Control-Allow-Headers", "Content-Type")
        response.headers.add("Access-Control-Allow-Methods", "POST, OPTIONS")
        return response
    
    try:
        data = request.json
        step2_json_path = data.get("step2_json_path")
        step2_json_data = data.get("step2_json_data")
        output_base_dir = data.get("output_base_dir", "./workflows")
        binder_length = data.get("binder_length", 90)
        num_backbones = data.get("num_backbones", 10)
        num_sequences_per_backbone = data.get("num_sequences_per_backbone", 8)
        num_af_models = data.get("num_af_models", 5)
        num_af_recycles = data.get("num_af_recycles", 3)
        use_colabfold = data.get("use_colabfold", True)
        rfdiffusion_path = data.get("rfdiffusion_path")
        proteinmpnn_path = data.get("proteinmpnn_path")
        colabfold_path = data.get("colabfold_path")
        alphafold_path = data.get("alphafold_path")
        execute = data.get("execute", False)
        
        # Either step2_json_path or step2_json_data must be provided
        if not step2_json_path and not step2_json_data:
            return jsonify({"error": "step2_json_path or step2_json_data is required"}), 400
        
        # Import Step 3 module
        from exo_gpt.step3_nanobinder_orchestrator import generate_design_plan
        
        # If JSON data is provided directly, write it to a temp file
        temp_file_created = False
        if step2_json_data:
            import tempfile
            with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
                json.dump(step2_json_data, f)
                step2_json_path = f.name
                temp_file_created = True
        else:
            # Ensure step2_json_path is absolute and normalize path
            if step2_json_path:
                # Remove leading ./ if present
                if step2_json_path.startswith('./'):
                    step2_json_path = step2_json_path[2:]
                # Make absolute if relative
                if not os.path.isabs(step2_json_path):
                    step2_json_path = os.path.join(PROJECT_ROOT, step2_json_path)
                # Normalize path
                step2_json_path = os.path.normpath(step2_json_path)
        
        if not os.path.exists(step2_json_path):
            return jsonify({"error": f"Step 2 JSON file not found: {step2_json_path}"}), 400
        
        # Ensure output_base_dir is absolute
        if not os.path.isabs(output_base_dir):
            output_base_dir = os.path.join(PROJECT_ROOT, output_base_dir)
        
        # Queue for progress messages
        progress_queue = queue.Queue()
        plan_result = {"plan": None, "error": None}
        
        def progress_callback(message: str, percent: float):
            """Progress callback that sends messages via queue"""
            progress_queue.put({"message": message, "percent": percent})
        
        def run_generation():
            """Run the design plan generation in a separate thread"""
            try:
                plan = generate_design_plan(
                    step2_json_path=step2_json_path,
                    output_base_dir=output_base_dir,
                    binder_length=binder_length,
                    num_backbones=num_backbones,
                    num_sequences_per_backbone=num_sequences_per_backbone,
                    num_af_models=num_af_models,
                    num_af_recycles=num_af_recycles,
                    use_colabfold=use_colabfold,
                    rfdiffusion_path=rfdiffusion_path,
                    proteinmpnn_path=proteinmpnn_path,
                    colabfold_path=colabfold_path,
                    alphafold_path=alphafold_path,
                    execute=execute,
                    progress_callback=progress_callback if execute else None,
                )
                plan_result["plan"] = plan
            except Exception as e:
                import traceback
                plan_result["error"] = {
                    "message": str(e),
                    "traceback": traceback.format_exc()
                }
            finally:
                # Signal completion
                progress_queue.put(None)
                # Clean up temp file if created
                if temp_file_created and os.path.exists(step2_json_path):
                    try:
                        os.unlink(step2_json_path)
                    except:
                        pass
        
        # Start generation in background thread
        thread = threading.Thread(target=run_generation, daemon=True)
        thread.start()
        
        def generate():
            """Generator function for SSE streaming"""
            try:
                completed = False
                while not completed:
                    try:
                        # Get message from queue with timeout
                        item = progress_queue.get(timeout=1.0)
                        
                        # None signals completion
                        if item is None:
                            completed = True
                            # Wait a bit for plan_result to be set
                            import time
                            for _ in range(10):  # Wait up to 1 second
                                if plan_result["plan"] is not None or plan_result["error"] is not None:
                                    break
                                time.sleep(0.1)
                            
                            # Send final result
                            if plan_result["plan"]:
                                yield f"data: {json.dumps({'type': 'result', 'plan': plan_result['plan']})}\n\n"
                            elif plan_result["error"]:
                                yield f"data: {json.dumps({'type': 'error', 'error': plan_result['error']})}\n\n"
                            else:
                                yield f"data: {json.dumps({'type': 'error', 'error': {'message': 'Generation completed but no result available'}})}\n\n"
                            break
                        
                        # Send progress message
                        yield f"data: {json.dumps({'type': 'progress', 'message': item['message'], 'percent': item['percent']})}\n\n"
                    except queue.Empty:
                        # Send keepalive to keep connection alive
                        yield ": keepalive\n\n"
                        # Check if thread is still alive
                        if not thread.is_alive() and not completed:
                            # Thread finished, wait for result
                            import time
                            time.sleep(0.2)
                            if plan_result["plan"] is not None:
                                yield f"data: {json.dumps({'type': 'result', 'plan': plan_result['plan']})}\n\n"
                                completed = True
                            elif plan_result["error"] is not None:
                                yield f"data: {json.dumps({'type': 'error', 'error': plan_result['error']})}\n\n"
                                completed = True
                        continue
            except Exception as e:
                import traceback
                yield f"data: {json.dumps({'type': 'error', 'error': {'message': str(e), 'traceback': traceback.format_exc()}})}\n\n"
        
        response = Response(
            stream_with_context(generate()),
            mimetype='text/event-stream',
            headers={
                'Cache-Control': 'no-cache',
                'X-Accel-Buffering': 'no',
                'Connection': 'keep-alive',
                'Access-Control-Allow-Origin': '*',
                'Access-Control-Allow-Headers': 'Content-Type',
                'Access-Control-Allow-Methods': 'POST, OPTIONS'
            }
        )
        return response
        
    except Exception as e:
        import traceback
        # For streaming endpoint, return JSON error with CORS headers
        error_response = jsonify({
            "success": False,
            "error": str(e),
            "traceback": traceback.format_exc()
        })
        error_response.headers.add("Access-Control-Allow-Origin", "*")
        return error_response, 500


@app.route("/api/step3/run", methods=["POST"])
def run_step3():
    """Execute Step 3: Nanobinder Design Orchestrator (non-streaming, for non-execute mode)"""
    try:
        data = request.json
        step2_json_path = data.get("step2_json_path")
        step2_json_data = data.get("step2_json_data")  # Alternative: pass JSON directly
        output_base_dir = data.get("output_base_dir", "./workflows")
        binder_length = data.get("binder_length", 90)
        num_backbones = data.get("num_backbones", 10)
        num_sequences_per_backbone = data.get("num_sequences_per_backbone", 8)
        num_af_models = data.get("num_af_models", 5)
        num_af_recycles = data.get("num_af_recycles", 3)
        use_colabfold = data.get("use_colabfold", True)
        rfdiffusion_path = data.get("rfdiffusion_path")
        proteinmpnn_path = data.get("proteinmpnn_path")
        colabfold_path = data.get("colabfold_path")
        alphafold_path = data.get("alphafold_path")
        execute = data.get("execute", False)  # New: execute tools if available
        
        # Either step2_json_path or step2_json_data must be provided
        if not step2_json_path and not step2_json_data:
            return jsonify({"error": "step2_json_path or step2_json_data is required"}), 400
        
        # Progress callback for API (store messages in list)
        progress_messages = []
        
        def progress_callback(message: str, percent: float):
            progress_messages.append({"message": message, "percent": percent})
        
        # Import Step 3 module
        from exo_gpt.step3_nanobinder_orchestrator import generate_design_plan
        
        # If JSON data is provided directly, write it to a temp file
        if step2_json_data:
            import tempfile
            with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
                json.dump(step2_json_data, f)
                step2_json_path = f.name
        else:
            # Ensure step2_json_path is absolute and normalize path
            if step2_json_path:
                # Remove leading ./ if present
                if step2_json_path.startswith('./'):
                    step2_json_path = step2_json_path[2:]
                # Make absolute if relative
                if not os.path.isabs(step2_json_path):
                    step2_json_path = os.path.join(PROJECT_ROOT, step2_json_path)
                # Normalize path
                step2_json_path = os.path.normpath(step2_json_path)
        
        if not os.path.exists(step2_json_path):
            return jsonify({"error": f"Step 2 JSON file not found: {step2_json_path}"}), 400
        
        # Ensure output_base_dir is absolute
        if not os.path.isabs(output_base_dir):
            output_base_dir = os.path.join(PROJECT_ROOT, output_base_dir)
        
        # Generate design plan
        plan = generate_design_plan(
            step2_json_path=step2_json_path,
            output_base_dir=output_base_dir,
            binder_length=binder_length,
            num_backbones=num_backbones,
            num_sequences_per_backbone=num_sequences_per_backbone,
            num_af_models=num_af_models,
            num_af_recycles=num_af_recycles,
            use_colabfold=use_colabfold,
            rfdiffusion_path=rfdiffusion_path,
            proteinmpnn_path=proteinmpnn_path,
            colabfold_path=colabfold_path,
            alphafold_path=alphafold_path,
            execute=execute,
            progress_callback=progress_callback if execute else None,
        )
        
        # Clean up temp file if created
        if step2_json_data and os.path.exists(step2_json_path):
            try:
                os.unlink(step2_json_path)
            except:
                pass
        
        return jsonify({
            "success": plan.get("success", False),
            "plan": plan,
            "error": plan.get("error") if not plan.get("success") else None,
            "progress": progress_messages if execute else []
        })
        
    except Exception as e:
        import traceback
        return jsonify({
            "success": False,
            "error": str(e),
            "traceback": traceback.format_exc()
        }), 500


@app.route("/api/file/read", methods=["POST"])
def read_file():
    """Read file contents for display in GUI"""
    try:
        data = request.json
        file_path = data.get("file_path")
        
        if not file_path:
            return jsonify({"error": "file_path is required"}), 400
        
        # Ensure file_path is absolute and normalize path
        if file_path.startswith('./'):
            file_path = file_path[2:]
        if not os.path.isabs(file_path):
            file_path = os.path.join(PROJECT_ROOT, file_path)
        file_path = os.path.normpath(file_path)
        
        # Security: ensure file is within project root
        if not file_path.startswith(str(PROJECT_ROOT)):
            return jsonify({"error": "File path outside project root"}), 403
        
        if not os.path.exists(file_path):
            return jsonify({"error": f"File not found: {file_path}"}), 404
        
        if not os.path.isfile(file_path):
            return jsonify({"error": "Path is not a file"}), 400
        
        # Read file contents
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                contents = f.read()
        except UnicodeDecodeError:
            # Try binary mode for non-text files
            with open(file_path, 'rb') as f:
                contents = f.read().decode('utf-8', errors='replace')
        
        # Determine file type for syntax highlighting
        file_ext = os.path.splitext(file_path)[1].lower()
        file_type = "text"
        if file_ext in ['.sh', '.bash']:
            file_type = "bash"
        elif file_ext in ['.py']:
            file_type = "python"
        elif file_ext in ['.yaml', '.yml']:
            file_type = "yaml"
        elif file_ext in ['.json']:
            file_type = "json"
        elif file_ext in ['.md']:
            file_type = "markdown"
        
        return jsonify({
            "success": True,
            "contents": contents,
            "file_type": file_type,
            "file_name": os.path.basename(file_path),
            "file_path": file_path
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

