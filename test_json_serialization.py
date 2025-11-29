#!/usr/bin/env python3
"""Test JSON serialization with NaN values"""

import pandas as pd
import numpy as np
import json
import os

excel_path = 'data/publications/publications_data.xlsx'

if not os.path.exists(excel_path):
    print("Excel file not found")
    exit(1)

# Read biomarkers
df_bio = pd.read_excel(excel_path, sheet_name='Biomarkers')
print("Biomarkers with NaN values:")
print(df_bio[['gene_symbol', 'logfc', 'p_value', 'fdr']].head().to_string())

# Test JSON serialization
def clean_for_json(obj):
    """Recursively clean NaN/None values for JSON serialization"""
    if isinstance(obj, dict):
        return {k: clean_for_json(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [clean_for_json(item) for item in obj]
    elif isinstance(obj, float):
        if pd.isna(obj) or np.isnan(obj):
            return None
    elif isinstance(obj, type(pd.NA)) if hasattr(pd, 'NA') else False:
        return None
    elif pd.isna(obj) if hasattr(pd, 'isna') else False:
        return None
    return obj

test_data = df_bio.head().to_dict('records')
print("\nBefore cleaning (has NaN):")
for item in test_data[:2]:
    print(f"  {item}")

cleaned = clean_for_json(test_data)
print("\nAfter cleaning:")
for item in cleaned[:2]:
    print(f"  {item}")

# Test JSON serialization
try:
    json_str = json.dumps(cleaned, indent=2)
    print("\n✓ JSON serialization successful!")
    print(f"JSON length: {len(json_str)} characters")
    print("\nSample JSON (first 300 chars):")
    print(json_str[:300])
except Exception as e:
    print(f"\n❌ JSON serialization failed: {e}")

