#!/usr/bin/env python3
"""
Helper script to process a publication file using OpenAI biomarker extractor.

Usage:
    python process_publication.py /path/to/publication.pdf
    python process_publication.py data/publications/my_paper.pdf
"""

import os
import sys
from pathlib import Path

# Add project root to path
PROJECT_ROOT = Path(__file__).parent
sys.path.insert(0, str(PROJECT_ROOT))

def main():
    if len(sys.argv) < 2:
        print("Usage: python process_publication.py <path_to_publication_file>")
        print("\nExample:")
        print("  python process_publication.py data/publications/my_paper.pdf")
        sys.exit(1)
    
    file_path = sys.argv[1]
    
    if not os.path.exists(file_path):
        print(f"Error: File not found: {file_path}")
        sys.exit(1)
    
    # Check for API key
    if not os.getenv("OPENAI_API_KEY"):
        print("Error: OPENAI_API_KEY environment variable not set!")
        print("\nTo set it:")
        print("  export OPENAI_API_KEY='sk-your-api-key-here'")
        print("\nOr add it to your .env file or ~/.bashrc")
        sys.exit(1)
    
    # Get publications directory
    publications_dir = os.path.join(PROJECT_ROOT, "data", "publications")
    os.makedirs(publications_dir, exist_ok=True)
    
    # If file is not in publications directory, copy it there
    if not file_path.startswith(publications_dir):
        import shutil
        filename = os.path.basename(file_path)
        dest_path = os.path.join(publications_dir, filename)
        
        # Avoid overwriting
        if os.path.exists(dest_path):
            import time
            name, ext = os.path.splitext(filename)
            filename = f"{name}_{int(time.time())}{ext}"
            dest_path = os.path.join(publications_dir, filename)
        
        print(f"Copying {file_path} to {dest_path}...")
        shutil.copy2(file_path, dest_path)
        file_path = dest_path
    
    # Process the file
    try:
        from exo_gpt.publication_extractor import extract_publication_data
        
        print(f"\nProcessing: {os.path.basename(file_path)}")
        print("This may take a few minutes depending on file size...\n")
        
        result = extract_publication_data(
            file_path,
            publications_dir,
            excel_filename="publications_data.xlsx"
        )
        
        # Display results
        print("\n" + "="*60)
        print("✓ Extraction Complete!")
        print("="*60)
        
        pub_info = result.get('publication', {})
        biomarkers = result.get('biomarkers', [])
        
        print(f"\nPublication Information:")
        if pub_info.get('title'):
            print(f"  Title: {pub_info['title']}")
        if pub_info.get('journal'):
            print(f"  Journal: {pub_info['journal']}")
        if pub_info.get('authors'):
            print(f"  Authors: {pub_info['authors']}")
        if pub_info.get('year'):
            print(f"  Year: {pub_info['year']}")
        
        print(f"\nBiomarkers Found: {len(biomarkers)}")
        if biomarkers:
            print("\nTop 10 Biomarkers:")
            for i, bm in enumerate(biomarkers[:10], 1):
                gene = bm.get('gene_symbol', 'N/A')
                uniprot = bm.get('uniprot', 'N/A')
                disease = bm.get('disease', 'N/A')
                print(f"  {i}. {gene} ({uniprot}) - {disease}")
            if len(biomarkers) > 10:
                print(f"  ... and {len(biomarkers) - 10} more")
        
        excel_path = os.path.join(publications_dir, "publications_data.xlsx")
        print(f"\n✓ Data saved to: {excel_path}")
        print("\nYou can now view this publication in the GUI!")
        
    except ImportError as e:
        print(f"Error: Failed to import publication extractor: {e}")
        print("Make sure you've installed all dependencies:")
        print("  pip install openai openpyxl PyPDF2 python-docx")
        sys.exit(1)
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()

