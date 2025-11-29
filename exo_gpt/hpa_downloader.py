"""
HPA (Human Protein Atlas) Pathology Data Downloader
-----------------------------------------------------

Downloads pathology data from HPA API for cancer types.
API endpoint: https://www.proteinatlas.org/api/search_download.php
"""

from __future__ import annotations
import os
import json
import re
from typing import Optional, Tuple
import requests
import pandas as pd


def slugify(text: str) -> str:
    """Convert text to filesystem-safe slug."""
    # Convert to lowercase and replace spaces/special chars with underscores
    slug = text.lower().strip()
    slug = re.sub(r'[^\w\s-]', '', slug)  # Remove special characters
    slug = re.sub(r'[-\s]+', '_', slug)   # Replace spaces and hyphens with underscore
    slug = slug.strip('_')                 # Remove leading/trailing underscores
    return slug


def download_hpa_pathology(
    cancer_name: str,
    output_base_dir: str = "data/hpa_raw",
    api_url: str = "https://www.proteinatlas.org/api/search_download.php"
) -> Tuple[Optional[str], Optional[str]]:
    """
    Download HPA pathology data for a given cancer type.
    
    Args:
        cancer_name: Human-readable cancer name (e.g., "colorectal cancer", "breast cancer")
        output_base_dir: Base directory for output files (default: "data/hpa_raw")
        api_url: HPA API endpoint URL
    
    Returns:
        Tuple of (json_path, csv_path) if successful, (None, None) if failed
    
    Example:
        >>> json_path, csv_path = download_hpa_pathology("colorectal cancer")
        >>> print(f"Saved JSON to: {json_path}")
        >>> print(f"Saved CSV to: {csv_path}")
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_base_dir, exist_ok=True)
    
    # Generate filesystem-safe slug
    slug = slugify(cancer_name)
    
    # Define output paths
    json_path = os.path.join(output_base_dir, f"{slug}_hpa_raw.json")
    csv_path = os.path.join(output_base_dir, f"{slug}_hpa_proteins.csv")
    
    # Prepare API request
    params = {
        "search": cancer_name,
        "format": "json"
    }
    
    try:
        # Call HPA API
        print(f"Downloading HPA pathology data for: {cancer_name}")
        print(f"API URL: {api_url}")
        print(f"Parameters: {params}")
        
        response = requests.get(api_url, params=params, timeout=30)
        
        # Check if response is HTML (error page) instead of JSON
        content_type = response.headers.get('content-type', '').lower()
        if 'text/html' in content_type or response.text.strip().startswith('<!DOCTYPE') or response.text.strip().startswith('<html'):
            error_msg = f"API returned HTML error page instead of JSON. The search term '{cancer_name}' may not be valid for HPA API."
            print(f"❌ {error_msg}")
            print(f"   Response preview: {response.text[:300] if response.text else 'No response text'}")
            print(f"   Try using HPA's exact terminology (e.g., 'Skin Cutaneous Melanoma' for melanoma)")
            return None, None
        
        # Check for 400 Bad Request - might indicate invalid search term
        if response.status_code == 400:
            error_msg = f"400 Bad Request - The search term '{cancer_name}' may not be valid for HPA API."
            print(f"❌ {error_msg}")
            print(f"   Response text: {response.text[:200] if response.text else 'No response text'}")
            print(f"   Try using HPA's exact terminology (e.g., 'Skin Cutaneous Melanoma' for melanoma, 'Colorectal cancer' for CRC)")
            return None, None
        
        response.raise_for_status()
        
        # Check if response is empty
        if not response.text or len(response.text.strip()) == 0:
            print(f"⚠️  Warning: Empty response from HPA API for '{cancer_name}'")
            return None, None
        
        # Try to parse as JSON
        try:
            data = response.json()
        except json.JSONDecodeError as e:
            print(f"❌ Error: Invalid JSON response from HPA API for '{cancer_name}': {e}")
            print(f"   Response preview: {response.text[:200] if response.text else 'No response text'}")
            return None, None
        
        # Check if data is empty or invalid
        if not data:
            print(f"⚠️  Warning: Empty data returned from HPA API for '{cancer_name}'")
            return None, None
        
        # Save raw JSON
        with open(json_path, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=2, ensure_ascii=False)
        print(f"✓ Saved raw JSON to: {json_path}")
        
        # Convert to DataFrame - handle different JSON formats
        try:
            # Try reading as JSON first
            if isinstance(data, list):
                # If it's a list of records, convert directly
                df = pd.DataFrame(data)
            elif isinstance(data, dict):
                # If it's a dict, try to find the data array
                # HPA API might return {"data": [...]} or similar
                if "data" in data and isinstance(data["data"], list):
                    df = pd.DataFrame(data["data"])
                elif "results" in data and isinstance(data["results"], list):
                    df = pd.DataFrame(data["results"])
                else:
                    # Try to convert dict to DataFrame (might be a single record)
                    df = pd.DataFrame([data])
            else:
                # Fallback: try pd.read_json
                df = pd.read_json(json_path)
        except Exception as e:
            print(f"❌ Error: Failed to parse JSON as DataFrame: {e}")
            print(f"   JSON data type: {type(data)}")
            if isinstance(data, dict):
                print(f"   JSON keys: {list(data.keys())[:10]}")
            return json_path, None
        
        # Check if DataFrame is empty
        if df.empty:
            print(f"⚠️  Warning: Empty DataFrame for '{cancer_name}' - no proteins found")
            return json_path, None
        
        print(f"✓ Parsed JSON: {len(df)} rows, {len(df.columns)} columns")
        
        # Preserve all columns as returned by HPA
        # Save as CSV
        df.to_csv(csv_path, index=False, encoding='utf-8')
        print(f"✓ Saved CSV ({len(df)} rows) to: {csv_path}")
        
        return json_path, csv_path
        
    except requests.exceptions.HTTPError as e:
        if e.response and e.response.status_code == 400:
            print(f"❌ Error: 400 Bad Request - The search term '{cancer_name}' is not valid for HPA API.")
            print(f"   Common valid terms: 'colorectal cancer', 'breast cancer', 'lung cancer', 'skin cancer', etc.")
            print(f"   Check HPA website for valid search terms: https://www.proteinatlas.org/")
        else:
            print(f"❌ Error: HTTP {e.response.status_code if e.response else 'Unknown'} - Failed to download HPA data for '{cancer_name}': {e}")
        return None, None
    except requests.exceptions.RequestException as e:
        print(f"❌ Error: Failed to download HPA data for '{cancer_name}': {e}")
        return None, None
    except Exception as e:
        print(f"❌ Error: Unexpected error while processing '{cancer_name}': {e}")
        return None, None


def download_multiple_cancers(
    cancer_names: list[str],
    output_base_dir: str = "data/hpa_raw"
) -> dict[str, Tuple[Optional[str], Optional[str]]]:
    """
    Download HPA pathology data for multiple cancer types.
    
    Args:
        cancer_names: List of cancer names to download
        output_base_dir: Base directory for output files
    
    Returns:
        Dictionary mapping cancer names to (json_path, csv_path) tuples
    """
    results = {}
    for cancer_name in cancer_names:
        json_path, csv_path = download_hpa_pathology(cancer_name, output_base_dir)
        results[cancer_name] = (json_path, csv_path)
        print()  # Empty line for readability
    
    return results


if __name__ == "__main__":
    # Example usage for multiple cancer types
    cancers = [
        "colorectal cancer",
        "breast cancer",
        "lung cancer",
        "prostate cancer",
        "melanoma"
    ]
    
    print("=" * 60)
    print("HPA Pathology Data Downloader")
    print("=" * 60)
    print()
    
    results = download_multiple_cancers(cancers)
    
    print()
    print("=" * 60)
    print("Download Summary")
    print("=" * 60)
    for cancer_name, (json_path, csv_path) in results.items():
        status = "✓ Success" if csv_path else "✗ Failed"
        print(f"{cancer_name:30s} {status}")
        if csv_path:
            print(f"  CSV: {csv_path}")

