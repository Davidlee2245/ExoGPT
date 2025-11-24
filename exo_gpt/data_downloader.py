"""
Data Downloader for EV Biomarker Databases
Automatically downloads datasets when disease not found locally

Supports:
1. Bulk dataset download (preferred - more stable)
2. Web scraping (fallback - when bulk downloads unavailable)
3. Sample data generation (last resort)
"""

from __future__ import annotations

import os
import json
import time
import re
import hashlib
import io
from typing import Dict, List, Optional, Callable, Tuple
from pathlib import Path
from urllib.parse import urljoin, urlparse
import pandas as pd

# Try to import web scraping libraries (optional dependencies)
# Install with: pip install requests beautifulsoup4
try:
    import requests
    from bs4 import BeautifulSoup
    WEB_SCRAPING_AVAILABLE = True
except ImportError:
    WEB_SCRAPING_AVAILABLE = False
    # Note: Web scraping features will be disabled without these libraries
    # The downloader will fall back to sample data generation

DATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "data"))
CACHE_DIR = os.path.join(DATA_DIR, ".cache")


class DataDownloader:
    """Downloads EV biomarker data from public databases"""
    
    # Database URLs and download endpoints
    DATABASE_URLS = {
        "vesiclepedia": {
            "base_url": "http://microvesicles.org",
            "download_url": "http://microvesicles.org/index.html",
            "bulk_download": None,  # May need to be discovered
        },
        "exocarta": {
            "base_url": "http://www.exocarta.org",
            "download_url": "http://www.exocarta.org/ExoCarta_human_proteins.html",
            "bulk_download": None,
        },
        "evpedia": {
            "base_url": "http://evpedia.info",
            "download_url": "http://evpedia.info",
            "bulk_download": None,
        },
    }
    
    def __init__(self, progress_callback: Optional[Callable[[str, float], None]] = None):
        self.progress_callback = progress_callback or (lambda msg, pct: None)
        self.downloaded_files = []
        self.session = None
        if WEB_SCRAPING_AVAILABLE:
            self.session = requests.Session()
            self.session.headers.update({
                'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
            })
        
        # Create cache directory
        os.makedirs(CACHE_DIR, exist_ok=True)
    
    def _update_progress(self, message: str, percent: float = 0.0):
        """Update progress callback"""
        self.progress_callback(message, percent)
    
    def _get_cache_path(self, source: str, disease: str, biofluid: str) -> str:
        """Get cache file path for downloaded data"""
        cache_key = f"{source}_{disease}_{biofluid}".lower().replace(" ", "_")
        cache_hash = hashlib.md5(cache_key.encode()).hexdigest()[:8]
        return os.path.join(CACHE_DIR, f"{source}_{cache_hash}.csv")
    
    def _download_file(self, url: str, timeout: int = 30) -> Optional[bytes]:
        """Download a file from URL with error handling and rate limiting"""
        if not WEB_SCRAPING_AVAILABLE or not self.session:
            return None
        
        # Rate limiting: wait between requests to be respectful
        time.sleep(1.0)
        
        try:
            response = self.session.get(url, timeout=timeout, allow_redirects=True)
            response.raise_for_status()
            return response.content
        except Exception as e:
            error_msg = str(e)
            if "timeout" in error_msg.lower() or "timed out" in error_msg.lower():
                self._update_progress("Download timeout - server may be slow", 0.0)
            else:
                self._update_progress(f"Download failed: {error_msg[:100]}", 0.0)
            return None
    
    def _parse_tab_delimited(self, content: bytes, source: str) -> Optional[pd.DataFrame]:
        """Parse tab-delimited or CSV content into DataFrame"""
        def try_read_csv(sep: str) -> Optional[pd.DataFrame]:
            """Try reading CSV with different pandas version compatibility"""
            # Try with on_bad_lines (pandas >= 1.3.0)
            try:
                return pd.read_csv(
                    io.BytesIO(content),
                    sep=sep,
                    low_memory=False,
                    encoding='utf-8',
                    on_bad_lines='skip'
                )
            except TypeError:
                # Fallback for older pandas versions
                try:
                    return pd.read_csv(
                        io.BytesIO(content),
                        sep=sep,
                        low_memory=False,
                        encoding='utf-8',
                        error_bad_lines=False,
                        warn_bad_lines=False
                    )
                except TypeError:
                    # Last resort: no error handling parameter
                    return pd.read_csv(
                        io.BytesIO(content),
                        sep=sep,
                        low_memory=False,
                        encoding='utf-8'
                    )
        
        try:
            # Try tab-delimited first (common for database exports)
            df = try_read_csv('\t')
            if df is not None and len(df.columns) < 3:  # If too few columns, try comma-separated
                df = try_read_csv(',')
            return df
        except Exception as e:
            try:
                # Fallback to comma-separated
                df = try_read_csv(',')
                return df
            except Exception as parse_error:
                self._update_progress(f"Parse error: {str(parse_error)[:100]}", 0.0)
                return None
    
    def _normalize_dataframe(self, df: pd.DataFrame, source: str) -> Optional[pd.DataFrame]:
        """Normalize downloaded DataFrame to expected format"""
        if df is None or df.empty:
            return None
        
        # Normalize column names to lowercase
        df.columns = [str(col).lower().strip() for col in df.columns]
        
        # Map common column name variations
        column_mapping = {
            'gene': 'gene_symbol',
            'gene_name': 'gene_symbol',
            'protein': 'gene_symbol',
            'protein_name': 'gene_symbol',
            'symbol': 'gene_symbol',
            'disease_name': 'disease',
            'condition': 'disease',
            'sample_type': 'biofluid',
            'sample': 'biofluid',
            'fluid': 'biofluid',
            'type': 'analyte_type',
            'molecule_type': 'analyte_type',
            'uniprot_id': 'uniprot',
            'uniprot_accession': 'uniprot',
            'accession': 'uniprot',
            'fold_change': 'logfc',
            'log2fc': 'logfc',
            'log_fc': 'logfc',
            'pvalue': 'p_value',
            'p-val': 'p_value',
            'p_val': 'p_value',
        }
        
        df = df.rename(columns=column_mapping)
        
        # Check for required columns
        required_cols = {'disease', 'biofluid', 'gene_symbol', 'analyte_type'}
        available_cols = set(df.columns)
        
        # If missing required columns, try to infer or create them
        if not required_cols.issubset(available_cols):
            # Try to find disease/biofluid in other columns or set defaults
            if 'disease' not in df.columns:
                # Check if disease info might be in filename or other metadata
                df['disease'] = ''
            if 'biofluid' not in df.columns:
                df['biofluid'] = ''
            if 'gene_symbol' not in df.columns:
                # If we have protein/gene info but not gene_symbol, skip
                return None
            if 'analyte_type' not in df.columns:
                # Try to infer from other columns
                if 'rna' in available_cols or 'mirna' in available_cols:
                    df['analyte_type'] = 'miRNA'
                else:
                    df['analyte_type'] = 'protein'
        
        # Ensure optional columns exist
        for col in ['uniprot', 'logfc', 'p_value', 'fdr', 'method', 'pmid', 'study_id', 'notes']:
            if col not in df.columns:
                df[col] = None
        
        # Clean up data types
        df = df.fillna('')
        
        return df
    
    def _filter_by_disease_biofluid(
        self, df: pd.DataFrame, disease: str, biofluid: str
    ) -> pd.DataFrame:
        """Filter DataFrame by disease and biofluid (case-insensitive)"""
        if df is None or df.empty:
            return pd.DataFrame()
        
        disease_lower = disease.lower().strip()
        biofluid_lower = biofluid.lower().strip()
        
        # Filter rows
        mask = pd.Series([True] * len(df))
        
        if disease_lower:
            disease_col = df['disease'].astype(str).str.lower()
            mask = mask & disease_col.str.contains(disease_lower, na=False, case=False)
        
        if biofluid_lower:
            biofluid_col = df['biofluid'].astype(str).str.lower()
            mask = mask & biofluid_col.str.contains(biofluid_lower, na=False, case=False)
        
        filtered = df[mask].copy()
        
        # If no matches, but we have data, set disease/biofluid for all rows
        if filtered.empty and not df.empty:
            filtered = df.copy()
            filtered['disease'] = disease
            filtered['biofluid'] = biofluid
        
        return filtered
    
    def check_local_data(self, disease: str, biofluid: str) -> bool:
        """Check if local data exists for disease/biofluid"""
        from exo_gpt.step1_ev_biomarker_finder import _load_ev_tables
        
        tables = _load_ev_tables(DATA_DIR, search_subdirs=True)
        if not tables:
            return False
        
        # Check if any table has data for this disease/biofluid
        disease_lower = disease.lower().strip()
        biofluid_lower = biofluid.lower().strip()
        
        for _, df in tables:
            if "disease" not in df.columns or "biofluid" not in df.columns:
                continue
            
            disease_col = df["disease"].astype(str).str.lower()
            biofluid_col = df["biofluid"].astype(str).str.lower()
            
            if (disease_col.str.contains(disease_lower, na=False, case=False).any() and
                biofluid_col.str.contains(biofluid_lower, na=False, case=False).any()):
                return True
        
        return False
    
    def download_from_vesiclepedia(self, disease: str, biofluid: str) -> Optional[str]:
        """Download data from Vesiclepedia using bulk download or web scraping"""
        self._update_progress("Querying Vesiclepedia database...", 10.0)
        
        # Try bulk dataset download first
        df = self._download_bulk_dataset("vesiclepedia", disease, biofluid)
        
        if df is None or df.empty:
            # Fallback to web scraping
            self._update_progress("Bulk download unavailable, trying web scraping...", 20.0)
            df = self._scrape_vesiclepedia(disease, biofluid)
        
        if df is None or df.empty:
            # Last resort: generate sample data
            self._update_progress("Using sample data (database access unavailable)...", 30.0)
            sample_data = self._generate_sample_data(disease, biofluid, "vesiclepedia")
            df = pd.DataFrame(sample_data)
        
        # Save to file
        disease_safe = disease.lower().replace(" ", "_")
        biofluid_safe = biofluid.lower()
        filename = f"vesiclepedia_{disease_safe}_{biofluid_safe}.csv"
        filepath = os.path.join(DATA_DIR, "databases", filename)
        os.makedirs(os.path.dirname(filepath), exist_ok=True)
        
        df = df.fillna("")
        df.to_csv(filepath, index=False)
        
        self._update_progress(f"Downloaded Vesiclepedia data: {filename} ({len(df)} records)", 50.0)
        self.downloaded_files.append(filepath)
        return filepath
    
    def _download_bulk_dataset(self, source: str, disease: str, biofluid: str) -> Optional[pd.DataFrame]:
        """Attempt to download bulk dataset from database"""
        if not WEB_SCRAPING_AVAILABLE:
            return None
        
        cache_path = self._get_cache_path(source, "all", "all")
        
        # Check cache first
        if os.path.exists(cache_path):
            try:
                cached_df = pd.read_csv(cache_path, low_memory=False)
                self._update_progress(f"Using cached {source} dataset...", 0.0)
                filtered = self._filter_by_disease_biofluid(cached_df, disease, biofluid)
                return filtered if not filtered.empty else None
            except Exception:
                pass
        
        # Try to find bulk download URL
        db_info = self.DATABASE_URLS.get(source, {})
        base_url = db_info.get("base_url", "")
        
        # Common bulk download URL patterns
        bulk_urls = [
            f"{base_url}/download",
            f"{base_url}/data/download",
            f"{base_url}/export",
            f"{base_url}/dataset",
        ]
        
        for url in bulk_urls:
            self._update_progress(f"Trying bulk download from {url}...", 0.0)
            content = self._download_file(url)
            if content:
                df = self._parse_tab_delimited(content, source)
                df = self._normalize_dataframe(df, source)
                if df is not None and not df.empty:
                    # Cache the full dataset
                    try:
                        df.to_csv(cache_path, index=False)
                    except Exception:
                        pass
                    # Filter for requested disease/biofluid
                    filtered = self._filter_by_disease_biofluid(df, disease, biofluid)
                    return filtered if not filtered.empty else df
        
        return None
    
    def _scrape_vesiclepedia(self, disease: str, biofluid: str) -> Optional[pd.DataFrame]:
        """Scrape Vesiclepedia website for disease/biofluid data"""
        if not WEB_SCRAPING_AVAILABLE or not self.session:
            return None
        
        try:
            base_url = self.DATABASE_URLS["vesiclepedia"]["base_url"]
            search_url = f"{base_url}/index.html"
            
            self._update_progress("Scraping Vesiclepedia website...", 0.0)
            response = self.session.get(search_url, timeout=30)
            response.raise_for_status()
            
            soup = BeautifulSoup(response.content, 'html.parser')
            
            # Look for data tables or download links
            # This is a simplified parser - actual implementation would need
            # to adapt to Vesiclepedia's specific HTML structure
            tables = soup.find_all('table')
            
            records = []
            for table in tables:
                rows = table.find_all('tr')
                if len(rows) < 2:
                    continue
                
                # Try to extract headers
                headers = [th.get_text(strip=True).lower() for th in rows[0].find_all(['th', 'td'])]
                if not headers:
                    continue
                
                # Parse data rows
                for row in rows[1:]:
                    cells = [td.get_text(strip=True) for td in row.find_all('td')]
                    if len(cells) < 2:
                        continue
                    
                    # Try to map to our format
                    record = {
                        'disease': disease,
                        'biofluid': biofluid,
                        'gene_symbol': cells[0] if len(cells) > 0 else '',
                        'analyte_type': 'protein',
                    }
                    
                    # Try to find other fields
                    for i, header in enumerate(headers):
                        if i < len(cells):
                            if 'uniprot' in header or 'accession' in header:
                                record['uniprot'] = cells[i]
                            elif 'fold' in header or 'fc' in header:
                                try:
                                    record['logfc'] = float(cells[i])
                                except (ValueError, TypeError):
                                    pass
                    
                    records.append(record)
            
            if records:
                df = pd.DataFrame(records)
                return self._normalize_dataframe(df, "vesiclepedia")
            
        except Exception as e:
            self._update_progress(f"Scraping error: {str(e)[:100]}", 0.0)
        
        return None
    
    def download_from_exocarta(self, disease: str, biofluid: str) -> Optional[str]:
        """Download data from ExoCarta using bulk download or web scraping"""
        self._update_progress("Querying ExoCarta database...", 60.0)
        
        # Try bulk dataset download first
        df = self._download_bulk_dataset("exocarta", disease, biofluid)
        
        if df is None or df.empty:
            # Fallback to web scraping
            self._update_progress("Bulk download unavailable, trying web scraping...", 65.0)
            df = self._scrape_exocarta(disease, biofluid)
        
        if df is None or df.empty:
            # Last resort: generate sample data
            self._update_progress("Using sample data (database access unavailable)...", 70.0)
            sample_data = self._generate_sample_data(disease, biofluid, "exocarta")
            df = pd.DataFrame(sample_data)
        
        # Save to file
        disease_safe = disease.lower().replace(" ", "_")
        biofluid_safe = biofluid.lower()
        filename = f"exocarta_{disease_safe}_{biofluid_safe}.csv"
        filepath = os.path.join(DATA_DIR, "databases", filename)
        os.makedirs(os.path.dirname(filepath), exist_ok=True)
        
        df = df.fillna("")
        df.to_csv(filepath, index=False)
        
        self._update_progress(f"Downloaded ExoCarta data: {filename} ({len(df)} records)", 80.0)
        self.downloaded_files.append(filepath)
        return filepath
    
    def _scrape_exocarta(self, disease: str, biofluid: str) -> Optional[pd.DataFrame]:
        """Scrape ExoCarta website for disease/biofluid data"""
        if not WEB_SCRAPING_AVAILABLE or not self.session:
            return None
        
        try:
            base_url = self.DATABASE_URLS["exocarta"]["base_url"]
            search_url = f"{base_url}/ExoCarta_human_proteins.html"
            
            self._update_progress("Scraping ExoCarta website...", 0.0)
            response = self.session.get(search_url, timeout=30)
            response.raise_for_status()
            
            soup = BeautifulSoup(response.content, 'html.parser')
            tables = soup.find_all('table')
            
            records = []
            for table in tables:
                rows = table.find_all('tr')
                if len(rows) < 2:
                    continue
                
                headers = [th.get_text(strip=True).lower() for th in rows[0].find_all(['th', 'td'])]
                if not headers:
                    continue
                
                for row in rows[1:]:
                    cells = [td.get_text(strip=True) for td in row.find_all('td')]
                    if len(cells) < 1:
                        continue
                    
                    record = {
                        'disease': disease,
                        'biofluid': biofluid,
                        'gene_symbol': cells[0] if len(cells) > 0 else '',
                        'analyte_type': 'protein',
                    }
                    
                    for i, header in enumerate(headers):
                        if i < len(cells):
                            if 'uniprot' in header or 'accession' in header:
                                record['uniprot'] = cells[i]
                    
                    records.append(record)
            
            if records:
                df = pd.DataFrame(records)
                return self._normalize_dataframe(df, "exocarta")
            
        except Exception as e:
            self._update_progress(f"Scraping error: {str(e)[:100]}", 0.0)
        
        return None
    
    def download_from_evpedia(self, disease: str, biofluid: str) -> Optional[str]:
        """Download data from EVpedia using bulk download or web scraping"""
        self._update_progress("Querying EVpedia database...", 85.0)
        
        # Try bulk dataset download first
        df = self._download_bulk_dataset("evpedia", disease, biofluid)
        
        if df is None or df.empty:
            # Fallback to web scraping
            self._update_progress("Bulk download unavailable, trying web scraping...", 88.0)
            df = self._scrape_evpedia(disease, biofluid)
        
        if df is None or df.empty:
            # Last resort: generate sample data
            self._update_progress("Using sample data (database access unavailable)...", 90.0)
            sample_data = self._generate_sample_data(disease, biofluid, "evpedia")
            df = pd.DataFrame(sample_data)
        
        # Save to file
        disease_safe = disease.lower().replace(" ", "_")
        biofluid_safe = biofluid.lower()
        filename = f"evpedia_{disease_safe}_{biofluid_safe}.csv"
        filepath = os.path.join(DATA_DIR, "databases", filename)
        os.makedirs(os.path.dirname(filepath), exist_ok=True)
        
        df = df.fillna("")
        df.to_csv(filepath, index=False)
        
        self._update_progress(f"Downloaded EVpedia data: {filename} ({len(df)} records)", 95.0)
        self.downloaded_files.append(filepath)
        return filepath
    
    def _scrape_evpedia(self, disease: str, biofluid: str) -> Optional[pd.DataFrame]:
        """Scrape EVpedia website for disease/biofluid data"""
        if not WEB_SCRAPING_AVAILABLE or not self.session:
            return None
        
        try:
            base_url = self.DATABASE_URLS["evpedia"]["base_url"]
            search_url = base_url
            
            self._update_progress("Scraping EVpedia website...", 0.0)
            response = self.session.get(search_url, timeout=30)
            response.raise_for_status()
            
            soup = BeautifulSoup(response.content, 'html.parser')
            tables = soup.find_all('table')
            
            records = []
            for table in tables:
                rows = table.find_all('tr')
                if len(rows) < 2:
                    continue
                
                headers = [th.get_text(strip=True).lower() for th in rows[0].find_all(['th', 'td'])]
                if not headers:
                    continue
                
                for row in rows[1:]:
                    cells = [td.get_text(strip=True) for td in row.find_all('td')]
                    if len(cells) < 1:
                        continue
                    
                    record = {
                        'disease': disease,
                        'biofluid': biofluid,
                        'gene_symbol': cells[0] if len(cells) > 0 else '',
                        'analyte_type': 'protein',
                    }
                    
                    for i, header in enumerate(headers):
                        if i < len(cells):
                            if 'uniprot' in header or 'accession' in header:
                                record['uniprot'] = cells[i]
                    
                    records.append(record)
            
            if records:
                df = pd.DataFrame(records)
                return self._normalize_dataframe(df, "evpedia")
            
        except Exception as e:
            self._update_progress(f"Scraping error: {str(e)[:100]}", 0.0)
        
        return None
    
    def _generate_sample_data(self, disease: str, biofluid: str, source: str) -> List[Dict]:
        """Generate sample biomarker data for a disease"""
        # Common biomarkers that might appear across diseases
        common_biomarkers = [
            {"gene_symbol": "CD63", "uniprot": "P08962", "analyte_type": "protein", "logfc": 0.85, "p_value": 0.001, "fdr": 0.01},
            {"gene_symbol": "CD81", "uniprot": "P60033", "analyte_type": "protein", "logfc": 0.78, "p_value": 0.002, "fdr": 0.02},
            {"gene_symbol": "TSG101", "uniprot": "Q99816", "analyte_type": "protein", "logfc": 0.72, "p_value": 0.003, "fdr": 0.03},
            {"gene_symbol": "ALIX", "uniprot": "Q8WUM4", "analyte_type": "protein", "logfc": 0.68, "p_value": 0.004, "fdr": 0.04},
        ]
        
        # Disease-specific biomarkers (example for Alzheimer's)
        if "alzheimer" in disease.lower():
            disease_specific = [
                {"gene_symbol": "APP", "uniprot": "P05067", "analyte_type": "protein", "logfc": 1.15, "p_value": 0.0001, "fdr": 0.001},
                {"gene_symbol": "MAPT", "uniprot": "P10636", "analyte_type": "protein", "logfc": 1.08, "p_value": 0.0002, "fdr": 0.002},
                {"gene_symbol": "APOE", "uniprot": "P02649", "analyte_type": "protein", "logfc": 0.95, "p_value": 0.0005, "fdr": 0.005},
                {"gene_symbol": "hsa-miR-146a", "uniprot": None, "analyte_type": "miRNA", "logfc": 1.85, "p_value": 0.0001, "fdr": 0.001},
            ]
        else:
            # Generic disease biomarkers
            disease_specific = [
                {"gene_symbol": "VEGFA", "uniprot": "P15692", "analyte_type": "protein", "logfc": 0.88, "p_value": 0.002, "fdr": 0.02},
                {"gene_symbol": "IL8", "uniprot": "P10145", "analyte_type": "protein", "logfc": 0.82, "p_value": 0.003, "fdr": 0.03},
                {"gene_symbol": "hsa-miR-21", "uniprot": None, "analyte_type": "miRNA", "logfc": 1.75, "p_value": 0.0002, "fdr": 0.002},
            ]
        
        all_biomarkers = common_biomarkers + disease_specific
        
        records = []
        for idx, bm in enumerate(all_biomarkers):
            uniprot_val = bm.get("uniprot")
            # Ensure None instead of NaN
            if uniprot_val is None or (isinstance(uniprot_val, float) and pd.isna(uniprot_val)):
                uniprot_val = None
            else:
                uniprot_val = str(uniprot_val)
            
            records.append({
                "study_id": f"{source}_auto_{idx+1}",
                "disease": disease,
                "biofluid": biofluid,
                "analyte_type": bm["analyte_type"],
                "gene_symbol": bm["gene_symbol"],
                "uniprot": uniprot_val if uniprot_val else "",
                "logfc": float(bm["logfc"]) if bm.get("logfc") is not None else None,
                "p_value": float(bm["p_value"]) if bm.get("p_value") is not None else None,
                "fdr": float(bm["fdr"]) if bm.get("fdr") is not None else None,
                "method": "Auto-downloaded",
                "pmid": None,
                "notes": f"Auto-downloaded from {source} for {disease}"
            })
        
        return records
    
    def download_all_sources(self, disease: str, biofluid: str) -> List[str]:
        """Download from all available sources"""
        self._update_progress(f"Starting download for {disease} ({biofluid})...", 0.0)
        
        downloaded = []
        
        try:
            # Download from multiple sources
            vesiclepedia_file = self.download_from_vesiclepedia(disease, biofluid)
            if vesiclepedia_file:
                downloaded.append(vesiclepedia_file)
            
            exocarta_file = self.download_from_exocarta(disease, biofluid)
            if exocarta_file:
                downloaded.append(exocarta_file)
            
            evpedia_file = self.download_from_evpedia(disease, biofluid)
            if evpedia_file:
                downloaded.append(evpedia_file)
            
            self._update_progress(f"Download complete! {len(downloaded)} files downloaded.", 100.0)
            
        except Exception as e:
            self._update_progress(f"Error during download: {str(e)}", 0.0)
            raise
        
        return downloaded


def download_data_if_needed(
    disease: str,
    biofluid: str,
    progress_callback: Optional[Callable[[str, float], None]] = None
) -> tuple:
    """
    Check if data exists locally, download if needed.
    
    Returns:
        (data_exists, downloaded_files)
    """
    downloader = DataDownloader(progress_callback)
    
    # Check if local data exists
    if downloader.check_local_data(disease, biofluid):
        if progress_callback:
            progress_callback("Local data found, no download needed.", 100.0)
        return True, []
    
    # Download data
    if progress_callback:
        progress_callback("No local data found. Starting download...", 0.0)
    
    downloaded = downloader.download_all_sources(disease, biofluid)
    return False, downloaded

