"""
Step 3-1 — Existing Antibody Finder (Commercial Lookup)
--------------------------------------------------------

Finds existing commercial antibodies that match targets/epitopes from Step 2.

Input: Step 2 JSON output (curated targets with epitopes)
Output: Commercial antibody matches per epitope
"""

from __future__ import annotations

import json
import os
import time
from dataclasses import dataclass, asdict
from typing import List, Dict, Any, Optional, Callable

# Try to import web scraping libraries
try:
    import requests
    from bs4 import BeautifulSoup
    WEB_SCRAPING_AVAILABLE = True
except ImportError:
    WEB_SCRAPING_AVAILABLE = False
    requests = None
    BeautifulSoup = None


@dataclass
class CommercialAntibody:
    """Commercial antibody information."""
    vendor: str
    product_id: str
    clone_name: Optional[str] = None
    target_gene: str = ""
    target_uniprot: Optional[str] = None
    host_species: str = ""
    isotype: Optional[str] = None
    applications: List[str] = None
    reported_epitope: Optional[str] = None
    epitope_residues: Optional[List[int]] = None
    url: Optional[str] = None
    validation_status: Optional[str] = None
    references: Optional[int] = None  # Number of references/validations
    
    def __post_init__(self):
        if self.applications is None:
            self.applications = []


@dataclass
class AntibodyMatch:
    """Match result for a single epitope."""
    epitope_id: str
    match_type: str
    commercial_antibodies: List[CommercialAntibody]
    structural_antibodies: List[Dict[str, Any]]  # Empty for now
    match_confidence: str
    notes: Optional[str] = None


@dataclass
class AntibodyFinderResult:
    """Complete antibody finder output."""
    target_name: str
    gene_symbol: str
    uniprot: Optional[str]
    epitope_matches: List[AntibodyMatch]
    search_metadata: Dict[str, Any]


class CommercialAntibodyFinder:
    """Finds commercial antibodies for given targets."""
    
    def __init__(self, progress_callback: Optional[Callable[[str, float], None]] = None):
        self.progress_callback = progress_callback or (lambda msg, pct: None)
        self.session = None
        
        if WEB_SCRAPING_AVAILABLE:
            self.session = requests.Session()
            self.session.headers.update({
                'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36',
                'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8',
                'Accept-Language': 'en-US,en;q=0.5',
                'Accept-Encoding': 'gzip, deflate',
                'Connection': 'keep-alive',
                'Upgrade-Insecure-Requests': '1'
            })
    
    def _update_progress(self, message: str, percent: float = 0.0):
        """Update progress callback"""
        self.progress_callback(message, percent)
    
    def find_commercial_antibodies(
        self,
        gene_symbol: str,
        uniprot: Optional[str] = None,
        max_results: int = 50
    ) -> List[CommercialAntibody]:
        """
        Find commercial antibodies for a target gene.
        
        Args:
            gene_symbol: Gene symbol (e.g., "CD63", "CD274")
            uniprot: Optional UniProt ID
            max_results: Maximum number of results to return
        
        Returns:
            List of CommercialAntibody objects
        """
        if not WEB_SCRAPING_AVAILABLE:
            self._update_progress("Web scraping not available. Install beautifulsoup4: pip install beautifulsoup4", 0.0)
            return []
        
        antibodies = []
        
        # Always try gene symbol first (more reliable for Antibodypedia)
        # Only try UniProt if gene symbol search fails
        self._update_progress(f"Searching for antibodies targeting {gene_symbol}...", 0.0)
        
        # Strategy 1: Query Antibodypedia with gene symbol
        try:
            antibodypedia_results = self._query_antibodypedia(gene_symbol, max_results)
            antibodies.extend(antibodypedia_results)
            time.sleep(1)  # Rate limiting
        except Exception as e:
            self._update_progress(f"Antibodypedia query failed: {str(e)[:100]}", 0.0)
        
        # Strategy 2: Query vendor websites (BioLegend, CST, Abcam) with gene symbol
        if len(antibodies) < max_results:
            try:
                vendor_results = self._query_vendor_websites(gene_symbol, max_results - len(antibodies))
                antibodies.extend(vendor_results)
                time.sleep(1)  # Rate limiting
            except Exception as e:
                self._update_progress(f"Vendor query failed: {str(e)[:100]}", 0.0)
        
        # Deduplicate by vendor + product_id
        seen = set()
        unique_antibodies = []
        for ab in antibodies:
            key = (ab.vendor.lower(), ab.product_id.lower())
            if key not in seen:
                seen.add(key)
                unique_antibodies.append(ab)
        
        return unique_antibodies[:max_results]
    
    def _query_antibodypedia(self, query: str, max_results: int = 50) -> List[CommercialAntibody]:
        """Query Antibodypedia database."""
        antibodies = []
        
        try:
            # Step 1: Search for the gene
            # Antibodypedia search URL format - use gene symbol, not UniProt ID
            # Skip if query looks like a UniProt ID (starts with letter + numbers, e.g., P08962, Q9NZQ7)
            if query and len(query) >= 6 and query[0].isalpha() and query[1:].replace('-', '').isdigit():
                self._update_progress(f"Skipping UniProt ID {query} - using gene symbol instead", 0.0)
                print(f"DEBUG: Skipping UniProt ID {query}")
                return antibodies
            
            # Try direct explore URL first (based on user's manual testing)
            # The search redirects to /explore/{gene_lowercase}
            explore_url = f"https://www.antibodypedia.com/explore/{query.lower()}"
            print(f"DEBUG: Step 1 - Trying direct explore URL: {explore_url}")
            self._update_progress(f"Searching Antibodypedia for {query}...", 0.0)
            
            response = self.session.get(explore_url, timeout=15, allow_redirects=True)
            print(f"DEBUG: Explore URL response status: {response.status_code}")
            print(f"DEBUG: Final URL after redirect: {response.url}")
            
            # If explore URL doesn't work, try search URL as fallback
            if response.status_code == 404:
                print(f"DEBUG: Explore URL returned 404, trying search URL as fallback...")
                search_url = f"https://www.antibodypedia.com/search?q={query}"
                print(f"DEBUG: Trying search URL: {search_url}")
                response = self.session.get(search_url, timeout=15, allow_redirects=True)
                print(f"DEBUG: Search URL response status: {response.status_code}")
                print(f"DEBUG: Final URL after redirect: {response.url}")
            
            # Handle 404 gracefully - might mean no results or wrong URL format
            if response.status_code == 404:
                error_msg = f"Antibodypedia returned 404 for {query} - may not have this gene"
                self._update_progress(error_msg, 0.0)
                print(f"DEBUG: {error_msg}")
                return antibodies
            
            response.raise_for_status()
            
            # Check final URL
            final_url = response.url
            print(f"DEBUG: Parsing HTML from {final_url}")
            soup = BeautifulSoup(response.content, 'html.parser')
            print(f"DEBUG: HTML content length: {len(response.content)} bytes")
            
            # Step 2: Find the gene page link
            # Antibodypedia gene pages are at: /gene/{gene_id}/{gene_symbol}
            gene_page_link = None
            
            print(f"DEBUG: Step 2 - Looking for gene page link...")
            
            # Strategy 1: Look for links with /gene/ pattern
            all_links = soup.find_all('a', href=True)
            print(f"DEBUG: Found {len(all_links)} total links on page")
            
            gene_links = [link for link in all_links if '/gene/' in link.get('href', '')]
            print(f"DEBUG: Found {len(gene_links)} links containing '/gene/'")
            
            for link in gene_links:
                href = link.get('href', '')
                text = link.text.strip()
                print(f"DEBUG: Gene link found: {href} (text: {text})")
                if '/gene/' in href:
                    gene_page_link = href
                    print(f"DEBUG: Selected gene page link: {gene_page_link}")
                    break
            
            # Strategy 2: If on explore page, look for the gene name link in results table
            if not gene_page_link and '/explore/' in final_url:
                print(f"DEBUG: On explore page, searching in table rows...")
                # Look in table rows for gene links
                table_rows = soup.find_all('tr')
                print(f"DEBUG: Found {len(table_rows)} table rows")
                
                # Display table contents for debugging
                for row_idx, row in enumerate(table_rows):
                    print(f"DEBUG: --- Row {row_idx} ---")
                    cells = row.find_all(['td', 'th'])
                    for cell_idx, cell in enumerate(cells):
                        cell_text = cell.get_text(strip=True)
                        links = cell.find_all('a', href=True)
                        print(f"DEBUG:   Cell {cell_idx}: '{cell_text}'")
                        for link in links:
                            href = link.get('href', '')
                            link_text = link.get_text(strip=True)
                            print(f"DEBUG:     Link: href='{href}', text='{link_text}'")
                    
                    # Look for links that match the query
                    links = row.find_all('a', href=True)
                    for link in links:
                        href = link.get('href', '')
                        text = link.text.strip().upper()
                        print(f"DEBUG: Checking link: href='{href}', text='{text}', query='{query.upper()}'")
                        
                        # Look for links that match the query (might not have /gene/ in explore page)
                        if query.upper() in text:
                            # Check if it's a gene link or construct it
                            if '/gene/' in href:
                                gene_page_link = href
                                print(f"DEBUG: Found gene link in table: {gene_page_link}")
                                break
                            # If it's just the gene name, we might need to construct the URL
                            elif text == query.upper() and href:
                                # Try to see if this link leads to gene page
                                print(f"DEBUG: Found matching link, checking if it's a gene page link...")
                                # Follow the link to see where it goes
                                test_url = href if href.startswith('http') else f"https://www.antibodypedia.com{href}"
                                print(f"DEBUG: Would follow: {test_url}")
                                # For now, if the text matches exactly, assume it's the gene link
                                if not gene_page_link:
                                    gene_page_link = href
                                    print(f"DEBUG: Using as gene page link: {gene_page_link}")
                    
                    if gene_page_link:
                        break
            
            # Strategy 3: Try to find gene ID from page source or construct URL
            if not gene_page_link:
                print(f"DEBUG: Trying to extract gene ID from page HTML...")
                import re
                page_text = str(soup)
                
                # Look for patterns like "gene/2770" in the HTML
                gene_id_pattern = r'/gene/(\d+)/'
                matches = re.findall(gene_id_pattern, page_text)
                print(f"DEBUG: Found {len(matches)} gene ID matches in HTML: {matches}")
                
                # Also check JavaScript/JSON data
                script_tags = soup.find_all('script')
                for script in script_tags:
                    script_text = script.string or ""
                    # Look for patterns like geneId: 2770 or "gene_id": 2770
                    id_patterns = [
                        r'gene[Ii]d["\']?\s*[:=]\s*(\d+)',
                        r'gene_id["\']?\s*[:=]\s*(\d+)',
                        r'"id"\s*:\s*(\d+)',
                    ]
                    for pattern in id_patterns:
                        id_matches = re.findall(pattern, script_text)
                        if id_matches:
                            print(f"DEBUG: Found gene ID in script: {id_matches[0]}")
                            matches.extend(id_matches)
                
                # Strategy 4: Try common gene IDs for known EV biomarkers
                # This is a fallback - we can maintain a small mapping
                known_gene_ids = {
                    'CD63': '2770',
                    'CD81': '975',
                    'CD9': '928',
                    'CD274': '29126',  # PD-L1
                    'CD44': '960',
                }
                
                gene_id_match = None
                if matches:
                    gene_id_match = matches[0]
                    print(f"DEBUG: Using extracted gene ID: {gene_id_match}")
                elif query.upper() in known_gene_ids:
                    gene_id_match = known_gene_ids[query.upper()]
                    print(f"DEBUG: Using known gene ID mapping: {gene_id_match}")
                
                if gene_id_match:
                    gene_page_link = f"/gene/{gene_id_match}/{query}"
                    print(f"DEBUG: Constructed gene page link: {gene_page_link}")
                else:
                    print(f"DEBUG: Could not determine gene ID - will try direct gene page access")
                    # Last resort: try the direct gene page URL pattern (might work for some genes)
                    # We'll try this in the next step
            
            # If we still don't have a gene page link, try known gene IDs
            if not gene_page_link:
                known_gene_ids = {
                    'CD63': '2770',
                    'CD81': '975',
                    'CD9': '928',
                    'CD274': '29126',  # PD-L1
                    'CD44': '960',
                }
                if query.upper() in known_gene_ids:
                    gene_id = known_gene_ids[query.upper()]
                    gene_page_link = f"/gene/{gene_id}/{query}"
                    print(f"DEBUG: Using known gene ID {gene_id} for {query}")
            
            if not gene_page_link:
                error_msg = f"Could not find gene page link for {query} on Antibodypedia"
                self._update_progress(error_msg, 0.0)
                print(f"DEBUG: {error_msg}")
                print(f"DEBUG: Page title: {soup.title.string if soup.title else 'No title'}")
                # Save HTML for inspection
                with open(f"/tmp/antibodypedia_search_{query}.html", 'w') as f:
                    f.write(str(soup))
                print(f"DEBUG: Saved HTML to /tmp/antibodypedia_search_{query}.html for inspection")
                return antibodies
            
            # Step 3: Navigate to the gene page
            if not gene_page_link.startswith('http'):
                gene_page_url = f"https://www.antibodypedia.com{gene_page_link}"
            else:
                gene_page_url = gene_page_link
            
            time.sleep(1)  # Rate limiting
            self._update_progress(f"Navigating to gene page: {gene_page_url}", 0.0)
            print(f"DEBUG: Step 3 - Navigating to gene page: {gene_page_url}")
            response = self.session.get(gene_page_url, timeout=15)
            print(f"DEBUG: Gene page response status: {response.status_code}")
            response.raise_for_status()
            soup = BeautifulSoup(response.content, 'html.parser')
            print(f"DEBUG: Gene page HTML length: {len(response.content)} bytes")
            
            # Step 4: Parse ONLY "TOP VALIDATED ANTIBODIES" table on the gene page
            # DO NOT parse "All Antibodies" or other sections
            print(f"DEBUG: Step 4 - Looking for 'TOP VALIDATED ANTIBODIES' table ONLY")
            
            # Save HTML for inspection
            with open(f"/tmp/antibodypedia_gene_{query}.html", 'w', encoding='utf-8') as f:
                f.write(str(soup))
            print(f"DEBUG: Saved gene page HTML to /tmp/antibodypedia_gene_{query}.html")
            
            top_validated_table = None
            page_text = soup.get_text()
            
            # Strategy 1: Look for table by ID "featured_antibodies" (most reliable)
            print(f"DEBUG: Strategy 1 - Looking for table with id='featured_antibodies'")
            top_validated_table = soup.find('table', id='featured_antibodies')
            if top_validated_table:
                print(f"DEBUG: ✓ Found table by ID 'featured_antibodies'")
            
            # Strategy 2: Look for table where first row contains "Top validated antibodies"
            if not top_validated_table:
                print(f"DEBUG: Strategy 2 - Looking for table with 'Top validated antibodies' in first row")
                all_tables = soup.find_all('table')
                for table_idx, table in enumerate(all_tables):
                    first_row = table.find('tr')
                    if first_row:
                        first_row_text = first_row.get_text().upper()
                        if 'TOP VALIDATED' in first_row_text or 'TOP VALIDATED ANTIBODIES' in first_row_text:
                            # Check that this is not the "All Antibodies" section
                            # The "All Antibodies" section comes after "Top validated antibodies"
                            # Check if there's an h2 with "All Antibodies" before this table
                            prev_h2 = table.find_previous('h2')
                            if prev_h2 and 'ALL ANTIBODIES' in prev_h2.get_text().upper():
                                print(f"DEBUG: Table {table_idx + 1} has TOP VALIDATED but 'All Antibodies' h2 comes before it - skipping")
                                continue
                            top_validated_table = table
                            print(f"DEBUG: ✓ Found table (Table {table_idx + 1}) with 'Top validated antibodies' in first row")
                            break
            
            # Strategy 3: Fallback - look for table with Provider column that comes before "All Antibodies" h2
            if not top_validated_table:
                print(f"DEBUG: Strategy 3 - Fallback: Looking for first table with Provider-like structure")
                all_tables = soup.find_all('table')
                all_antibodies_h2 = soup.find('h2', string=lambda text: text and 'All Antibodies' in text)
                
                for table_idx, table in enumerate(all_tables):
                    # Check if this table comes before "All Antibodies" h2
                    if all_antibodies_h2:
                        # Check if table comes before the h2
                        table_pos = str(soup).find(str(table))
                        h2_pos = str(soup).find(str(all_antibodies_h2))
                        if h2_pos != -1 and table_pos != -1 and table_pos > h2_pos:
                            print(f"DEBUG: Table {table_idx + 1} comes after 'All Antibodies' h2 - skipping")
                            continue
                    
                    # Check if first data row (after header) has provider-like structure
                    rows = table.find_all('tr')
                    if len(rows) >= 2:
                        # First row might be header, second row should have data
                        data_row = rows[1] if rows[0].find('h2') or 'TOP VALIDATED' in rows[0].get_text().upper() else rows[0]
                        cells = data_row.find_all('td')
                        if len(cells) >= 3:
                            # First cell should be provider name (text, not empty)
                            first_cell_text = cells[0].get_text(strip=True)
                            if first_cell_text and len(first_cell_text) > 2 and not first_cell_text.isdigit():
                                # Second cell might be logo (empty or has class)
                                # Third cell should have antibody link
                                third_cell = cells[2] if len(cells) > 2 else None
                                if third_cell and third_cell.find('a', href=True):
                                    top_validated_table = table
                                    print(f"DEBUG: ✓ Found table (Table {table_idx + 1}) with provider-like structure")
                                    break
            
            if not top_validated_table:
                print(f"DEBUG: ERROR - Could not find TOP VALIDATED ANTIBODIES table")
                print(f"DEBUG: Page text sample (first 2000 chars): {page_text[:2000]}")
                print(f"DEBUG: Searching for 'TOP' in page text...")
                if 'TOP' in page_text.upper():
                    top_idx = page_text.upper().find('TOP')
                    print(f"DEBUG: Found 'TOP' at position {top_idx}, context: '{page_text[max(0, top_idx-50):top_idx+100]}'")
                return antibodies
            
            print(f"DEBUG: Found TOP VALIDATED ANTIBODIES table, parsing...")
            
            # Parse ONLY this table
            rows = top_validated_table.find_all('tr')
            print(f"DEBUG: TOP VALIDATED table has {len(rows)} rows")
            
            # Skip the first row if it contains "Top validated antibodies" header
            data_rows = []
            for row in rows:
                row_text = row.get_text().upper()
                # Skip header row that contains "TOP VALIDATED" or has colspan (header row)
                if 'TOP VALIDATED' in row_text or row.find('td', colspan=True):
                    print(f"DEBUG: Skipping header row: '{row.get_text(strip=True)[:50]}'")
                    continue
                data_rows.append(row)
            
            print(f"DEBUG: Found {len(data_rows)} data rows in TOP VALIDATED table")
            
            for row_idx, row in enumerate(data_rows):
                cells = row.find_all('td')
                # Need at least: Provider (0), Logo (1, skip), Product ID (2), References (3), Type (4), Applications (5)
                if len(cells) < 5:
                    print(f"DEBUG: Skipping row {row_idx} - only {len(cells)} cells (need at least 5)")
                    continue
                
                print(f"DEBUG: Processing TOP VALIDATED row {row_idx} with {len(cells)} cells")
                
                # Parse columns according to actual TOP VALIDATED ANTIBODIES format:
                # Column 0: Provider (vendor name)
                # Column 1: Provider logo (empty, skip this)
                # Column 2: Antibody ID / Catalog Number (with link)
                # Column 3: References (number with abbr tag)
                # Column 4: Type (Monoclonal/Polyclonal)
                # Column 5: Applications (complex div structure with images and spans)
                
                # Extract Provider (Column 0)
                vendor = cells[0].get_text(strip=True)
                print(f"DEBUG:   Provider: '{vendor}'")
                
                # Skip if vendor is "Unknown Provider" or empty
                if not vendor or vendor.upper() == "UNKNOWN PROVIDER" or vendor.upper().startswith("ENST"):
                    print(f"DEBUG:   Skipping - invalid vendor: '{vendor}'")
                    continue
                
                # Extract Product ID / Catalog Number (Column 2, skip logo column 1)
                product_id = None
                clone_name = None
                full_url = None
                if len(cells) > 2:
                    antibody_cell = cells[2]
                    link = antibody_cell.find('a', href=True)
                    if link:
                        product_id = link.get_text(strip=True)
                        clone_name = product_id  # In Antibodypedia, product ID often serves as clone name
                        href = link.get('href', '')
                        if not href.startswith('http'):
                            full_url = f"https://www.antibodypedia.com{href}" if href.startswith('/') else f"https://www.antibodypedia.com/{href}"
                        else:
                            full_url = href
                    else:
                        product_id = antibody_cell.get_text(strip=True)
                        clone_name = product_id
                
                if not product_id:
                    print(f"DEBUG:   Skipping - no product_id found")
                    continue
                
                # Skip if product_id looks like a transcript ID or UniProt ID
                if product_id.startswith('ENST') or product_id.startswith('ENSP') or (len(product_id) >= 6 and product_id[0].isalpha() and product_id[1:].replace('-', '').isdigit() and product_id[0] == 'P'):
                    print(f"DEBUG:   Skipping - product_id looks like transcript/UniProt: '{product_id}'")
                    continue
                
                print(f"DEBUG:   Catalog Number: '{product_id}'")
                
                # Extract References (Column 3)
                references = None
                if len(cells) > 3:
                    ref_cell = cells[3]
                    # Look for abbr tag with title containing "references"
                    abbr = ref_cell.find('abbr')
                    if abbr:
                        ref_text = abbr.get_text(strip=True)
                        try:
                            references = int(ref_text)
                        except ValueError:
                            pass
                    if references is None:
                        # Fallback: extract number from text
                        ref_text = ref_cell.get_text(strip=True)
                        import re
                        ref_match = re.search(r'(\d+)\s*references?', ref_text, re.IGNORECASE)
                        if ref_match:
                            references = int(ref_match.group(1))
                    print(f"DEBUG:   References: {references}")
                
                # Extract Type (Column 4)
                antibody_type = None
                if len(cells) > 4:
                    type_cell = cells[4]
                    type_text = type_cell.get_text(strip=True)
                    if type_text in ['Monoclonal', 'Polyclonal']:
                        antibody_type = type_text
                    print(f"DEBUG:   Type: '{antibody_type}'")
                
                # Extract Applications from Column 5 (complex div structure)
                applications = []
                if len(cells) > 5:
                    app_cell = cells[5]
                    # Look for spans with application codes (WB, EL, ICC, IP, IHC, FC, IEM, FA, OA)
                    # The structure is: <div><img ...><span>WB</span></div>
                    app_spans = app_cell.find_all('span')
                    for span in app_spans:
                        app_text = span.get_text(strip=True).upper()
                        # Check if this span has a corresponding image (validation indicator)
                        # The image should be in the same div or previous sibling
                        parent_div = span.find_parent('div')
                        if parent_div:
                            # Check if there's an img tag in this div (validation indicator)
                            img = parent_div.find('img')
                            if img or app_text in ['WB', 'EL', 'ICC', 'IP', 'IHC', 'FC', 'IEM', 'FA', 'OA']:
                                if app_text in ['WB', 'EL', 'ICC', 'IP', 'IHC', 'FC', 'IEM', 'FA', 'OA']:
                                    if app_text not in applications:
                                        applications.append(app_text)
                                        print(f"DEBUG:   Found application: {app_text}")
                
                # Only add if we have valid vendor and product_id
                if vendor and product_id and vendor.upper() != "UNKNOWN PROVIDER":
                    ab = CommercialAntibody(
                        vendor=vendor,
                        product_id=product_id,
                        clone_name=clone_name,
                        target_gene=query,
                        applications=applications,
                        isotype=antibody_type,
                        references=references,
                        url=full_url or gene_page_url
                    )
                    antibodies.append(ab)
                    print(f"DEBUG: ✓ Added: {vendor} | {product_id} | {references} refs | {antibody_type} | {', '.join(applications) if applications else 'N/A'}")
                    
                    if len(antibodies) >= max_results:
                        break
                else:
                    print(f"DEBUG: ✗ Skipped - missing vendor or product_id")
            
            print(f"DEBUG: Total antibodies extracted: {len(antibodies)}")
            
            # If we still don't have results, check if page mentions antibodies
            if len(antibodies) == 0:
                page_text = soup.get_text()
                import re
                antibody_count_match = re.search(r'(\d+)\s+antibod', page_text, re.IGNORECASE)
                if antibody_count_match:
                    # There are antibodies, but we couldn't parse them
                    # Return a placeholder indicating antibodies exist
                    ab = CommercialAntibody(
                        vendor="Antibodypedia",
                        product_id="multiple",
                        target_gene=query,
                        url=gene_page_url,
                        notes=f"Found {antibody_count_match.group(1)} antibodies on Antibodypedia (parsing needs improvement)"
                    )
                    antibodies.append(ab)
            
        except Exception as e:
            # Log the error for debugging
            import traceback
            error_msg = f"Antibodypedia query error: {str(e)[:200]}"
            self._update_progress(error_msg, 0.0)
            # Print to console for debugging (can be removed in production)
            print(f"DEBUG: ERROR - Antibodypedia query failed for {query}")
            print(f"DEBUG: Error message: {str(e)}")
            print(f"DEBUG: Full traceback:")
            print(traceback.format_exc())
        
        self._update_progress(f"Found {len(antibodies)} antibodies from Antibodypedia for {query}", 0.0)
        print(f"DEBUG: Returning {len(antibodies)} antibodies for {query}")
        return antibodies
    
    def _query_vendor_websites(self, query: str, max_results: int = 50) -> List[CommercialAntibody]:
        """Query vendor websites (BioLegend, CST, Abcam)."""
        antibodies = []
        
        vendors = [
            {
                "name": "BioLegend",
                "search_url": f"https://www.biolegend.com/en-us/search?q={query}",
                "base_url": "https://www.biolegend.com"
            },
            {
                "name": "Cell Signaling Technology",
                "search_url": f"https://www.cellsignal.com/products/primary-antibodies?N={query}",
                "base_url": "https://www.cellsignal.com"
            },
            {
                "name": "Abcam",
                "search_url": f"https://www.abcam.com/products?query={query}",
                "base_url": "https://www.abcam.com"
            }
        ]
        
        for vendor_info in vendors:
            if len(antibodies) >= max_results:
                break
            
            try:
                self._update_progress(f"Querying {vendor_info['name']}...", 0.0)
                response = self.session.get(vendor_info['search_url'], timeout=15)
                response.raise_for_status()
                
                soup = BeautifulSoup(response.content, 'html.parser')
                
                # Generic parsing - look for product links or listings
                # This is a basic implementation; vendor-specific parsing can be enhanced
                product_links = soup.find_all('a', href=lambda x: x and ('product' in x.lower() or 'antibody' in x.lower()))
                
                for link in product_links[:10]:  # Limit per vendor
                    try:
                        href = link.get('href', '')
                        if not href.startswith('http'):
                            href = vendor_info['base_url'] + href
                        
                        product_text = link.text.strip() if link.text else ""
                        product_id = href.split('/')[-1] if href else "N/A"
                        
                        ab = CommercialAntibody(
                            vendor=vendor_info['name'],
                            product_id=product_id,
                            clone_name=None,  # May need deeper parsing
                            target_gene=query,
                            url=href,
                            applications=[]  # May need deeper parsing
                        )
                        antibodies.append(ab)
                    except Exception:
                        continue
                
                time.sleep(1)  # Rate limiting between vendors
                
            except Exception as e:
                # Continue to next vendor if one fails
                continue
        
        return antibodies


def find_existing_antibodies(
    step2_json_path: str,
    output_json_path: Optional[str] = None,
    search_mode: str = "commercial",
    commercial_sources: Optional[List[str]] = None,
    epitope_overlap_threshold: float = 0.3,
    max_results_per_epitope: int = 50,
    progress_callback: Optional[Callable[[str, float], None]] = None
) -> AntibodyFinderResult:
    """
    Find existing antibodies for targets/epitopes from Step 2.
    
    Args:
        step2_json_path: Path to Step 2 JSON output
        output_json_path: Optional path to save results JSON
        search_mode: "commercial", "structural", or "hybrid" (only commercial implemented)
        commercial_sources: List of vendors to query (None = all available)
        epitope_overlap_threshold: Minimum overlap score for structural matches (not used yet)
        max_results_per_epitope: Maximum results to return per epitope
        progress_callback: Optional callback for progress updates
    
    Returns:
        AntibodyFinderResult with matches for each epitope
    """
    # Load Step 2 JSON
    with open(step2_json_path, 'r') as f:
        step2_data = json.load(f)
    
    if progress_callback:
        progress_callback("Starting antibody search...", 0.0)
    
    # Initialize finder
    finder = CommercialAntibodyFinder(progress_callback=progress_callback)
    
    # Process each target
    all_epitope_matches = []
    
    for target_idx, target in enumerate(step2_data.get("curated_targets", [])):
        gene_symbol = target.get("gene_symbol", "")
        uniprot = target.get("uniprot")
        target_name = target.get("target_name", gene_symbol)
        
        if progress_callback:
            progress = (target_idx + 1) / len(step2_data.get("curated_targets", [1]))
            progress_callback(f"Searching antibodies for {target_name} ({gene_symbol})...", progress * 0.8)
        
        # Find commercial antibodies for this target
        commercial_antibodies = finder.find_commercial_antibodies(
            gene_symbol=gene_symbol,
            uniprot=uniprot,
            max_results=max_results_per_epitope
        )
        
        # Create matches for each epitope
        for epitope in target.get("epitopes", []):
            epitope_id = epitope.get("epitope_id", f"epitope_{target_idx}")
            
            # For now, all epitopes from the same target get the same antibodies
            # In future, we can filter by epitope-specific criteria
            match = AntibodyMatch(
                epitope_id=epitope_id,
                match_type=search_mode,
                commercial_antibodies=commercial_antibodies,
                structural_antibodies=[],  # Not implemented yet
                match_confidence="high" if len(commercial_antibodies) > 5 else "medium" if len(commercial_antibodies) > 0 else "low",
                notes=f"Found {len(commercial_antibodies)} commercial antibodies for {gene_symbol}"
            )
            all_epitope_matches.append(match)
    
    # Create result (using first target for summary)
    first_target = step2_data.get("curated_targets", [{}])[0] if step2_data.get("curated_targets") else {}
    
    result = AntibodyFinderResult(
        target_name=first_target.get("target_name", "Unknown"),
        gene_symbol=first_target.get("gene_symbol", "Unknown"),
        uniprot=first_target.get("uniprot"),
        epitope_matches=all_epitope_matches,
        search_metadata={
            "search_mode": search_mode,
            "commercial_sources": commercial_sources or ["all"],
            "epitope_overlap_threshold": epitope_overlap_threshold,
            "max_results_per_epitope": max_results_per_epitope,
            "total_epitopes_searched": len(all_epitope_matches),
            "total_antibodies_found": sum(len(m.commercial_antibodies) for m in all_epitope_matches)
        }
    )
    
    # Save to JSON if path provided
    if output_json_path:
        result_dict = {
            "target_name": result.target_name,
            "gene_symbol": result.gene_symbol,
            "uniprot": result.uniprot,
            "epitope_matches": [
                {
                    "epitope_id": m.epitope_id,
                    "match_type": m.match_type,
                    "commercial_antibodies": [asdict(ab) for ab in m.commercial_antibodies],
                    "structural_antibodies": m.structural_antibodies,
                    "match_confidence": m.match_confidence,
                    "notes": m.notes
                }
                for m in result.epitope_matches
            ],
            "search_metadata": result.search_metadata
        }
        
        with open(output_json_path, 'w') as f:
            json.dump(result_dict, f, indent=2)
        
        if progress_callback:
            progress_callback(f"Results saved to {output_json_path}", 1.0)
    
    return result


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Find existing antibodies for Step 2 targets")
    parser.add_argument("--step2_json", required=True, help="Path to Step 2 JSON output")
    parser.add_argument("--output_json", help="Path to save results JSON")
    parser.add_argument("--search_mode", default="commercial", choices=["commercial", "structural", "hybrid"])
    parser.add_argument("--max_results_per_epitope", type=int, default=50)
    
    args = parser.parse_args()
    
    def progress_cb(msg, pct):
        print(f"[{pct*100:.1f}%] {msg}")
    
    result = find_existing_antibodies(
        step2_json_path=args.step2_json,
        output_json_path=args.output_json,
        search_mode=args.search_mode,
        max_results_per_epitope=args.max_results_per_epitope,
        progress_callback=progress_cb
    )
    
    print(f"\nFound {len(result.epitope_matches)} epitope matches")
    print(f"Total antibodies: {result.search_metadata['total_antibodies_found']}")

