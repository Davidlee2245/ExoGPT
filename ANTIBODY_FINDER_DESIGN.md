# Existing Antibody Finder - Design Document

## Executive Summary

This document describes the design and implementation strategy for an **"Existing Antibody Finder"** feature that complements Step 3 (Nanobinder Design Orchestrator) in the Exosome-GPT pipeline. The feature enables users to search for existing antibodies that bind to targets/epitopes identified in Step 2, providing an alternative to *de novo* nanobinder design.

---

## 1. Strategy Selection and Justification

### Recommended Approach: **Hybrid Strategy (Option 2 + Option 3)**

**Primary: Option 2 (Target-based Commercial Antibody Lookup)**
- **Why**: Most practical for experimental planning. Provides immediately actionable information (vendor, clone, applications) that researchers need for validation experiments.
- **Implementation**: Gene/UniProt-based queries to commercial antibody databases.

**Secondary: Option 3 (PDB Antibody-Antigen Mining with Epitope Overlap)**
- **Why**: Provides structural validation and epitope-specific matches when available. Complements commercial data with structural insights.
- **Implementation**: PDB mining for antibody-antigen complexes with residue overlap analysis.

**Why NOT Option 1 alone**: SAbDab/IEDB are excellent but require more complex structural alignment and may have limited coverage for all targets. Better as a future enhancement.

### Implementation Priority

1. **Phase 1**: Implement Option 2 (commercial antibody lookup) - Fast, practical, high user value
2. **Phase 2**: Add Option 3 (PDB structural matching) - Provides structural validation
3. **Phase 3** (Future): Integrate SAbDab/IEDB for comprehensive structural epitope matching

---

## 2. Module Architecture

### 2.1 Package Structure

```
exo_gpt/
├── step3_antibody_finder.py          # Main module
├── antibody_finder/
│   ├── __init__.py
│   ├── commercial_lookup.py          # Option 2: Commercial antibody databases
│   ├── pdb_mining.py                 # Option 3: PDB antibody-antigen mining
│   ├── epitope_matcher.py            # Epitope overlap computation
│   └── data_sources.py               # Database connectors and parsers
```

### 2.2 Core Data Classes

```python
from dataclasses import dataclass
from typing import List, Optional, Dict, Any

@dataclass
class CommercialAntibody:
    """Commercial antibody information."""
    vendor: str                    # e.g., "BioLegend", "CST", "Abcam"
    product_id: str                # Catalog number
    clone_name: Optional[str]      # Clone identifier
    target_gene: str                # Gene symbol
    target_uniprot: Optional[str]  # UniProt ID
    host_species: str               # e.g., "Mouse", "Rabbit"
    isotype: Optional[str]          # e.g., "IgG1", "IgG2a"
    applications: List[str]          # ["IHC", "WB", "FC", "ELISA"]
    reported_epitope: Optional[str]  # e.g., "Extracellular domain", "aa 50-100"
    epitope_residues: Optional[List[int]]  # If specific residues known
    url: Optional[str]              # Product page URL
    validation_status: Optional[str] # e.g., "Validated", "Uncharacterized"

@dataclass
class StructuralAntibody:
    """Antibody from PDB antibody-antigen complex."""
    pdb_id: str                     # PDB entry ID
    antibody_chain_ids: List[str]   # e.g., ["H", "L"]
    antigen_chain_id: str            # Antigen chain in complex
    antibody_name: Optional[str]     # If available
    clone_name: Optional[str]
    epitope_overlap_score: float     # Jaccard index or residue overlap
    overlapping_residues: List[int] # Residues in Step 2 epitope that overlap
    interface_residues: List[int]   # Antigen residues in contact with antibody
    resolution: Optional[float]     # Structure resolution
    method: str                      # "X-ray", "Cryo-EM", "NMR"

@dataclass
class AntibodyMatch:
    """Match result for a single epitope."""
    epitope_id: str
    match_type: str                  # "commercial" or "structural"
    commercial_antibodies: List[CommercialAntibody]
    structural_antibodies: List[StructuralAntibody]
    match_confidence: str            # "high", "medium", "low"
    notes: Optional[str]

@dataclass
class AntibodyFinderResult:
    """Complete antibody finder output."""
    target_name: str
    gene_symbol: str
    uniprot: Optional[str]
    epitope_matches: List[AntibodyMatch]
    search_metadata: Dict[str, Any]  # Timestamps, databases queried, etc.
```

### 2.3 Main Function Interface

```python
def find_existing_antibodies(
    step2_json_path: str,
    output_json_path: Optional[str] = None,
    search_mode: str = "hybrid",  # "commercial", "structural", "hybrid"
    commercial_sources: Optional[List[str]] = None,  # ["all"] or specific vendors
    epitope_overlap_threshold: float = 0.3,  # Minimum Jaccard index for structural matches
    max_results_per_epitope: int = 50,
    progress_callback: Optional[Callable[[str, float], None]] = None
) -> AntibodyFinderResult:
    """
    Find existing antibodies for targets/epitopes from Step 2.
    
    Args:
        step2_json_path: Path to Step 2 JSON output
        output_json_path: Optional path to save results JSON
        search_mode: "commercial", "structural", or "hybrid"
        commercial_sources: List of vendors to query (None = all available)
        epitope_overlap_threshold: Minimum overlap score for structural matches
        max_results_per_epitope: Maximum results to return per epitope
        progress_callback: Optional callback for progress updates
    
    Returns:
        AntibodyFinderResult with matches for each epitope
    """
    pass
```

### 2.4 Integration with Pipeline

**Input**: Step 2 JSON (unchanged structure)
```json
{
  "disease": "Metastatic melanoma",
  "biofluid": "Plasma",
  "curated_targets": [
    {
      "target_name": "PD-L1",
      "gene_symbol": "CD274",
      "uniprot": "Q9NZQ7",
      "epitopes": [
        {
          "epitope_id": "cd274_epitope_01",
          "epitope_residues": [52, 54, 55, 76, 77, 120, 121],
          "epitope_description": "Surface-exposed loop patch"
        }
      ]
    }
  ]
}
```

**Output**: Extended Step 2 JSON with `existing_antibodies` field
```json
{
  "disease": "Metastatic melanoma",
  "biofluid": "Plasma",
  "curated_targets": [
    {
      "target_name": "PD-L1",
      "gene_symbol": "CD274",
      "uniprot": "Q9NZQ7",
      "epitopes": [
        {
          "epitope_id": "cd274_epitope_01",
          "epitope_residues": [52, 54, 55, 76, 77, 120, 121],
          "existing_antibodies": {
            "commercial": [
              {
                "vendor": "BioLegend",
                "product_id": "329701",
                "clone_name": "29E.2A3",
                "applications": ["FC", "IHC"],
                "url": "https://..."
              }
            ],
            "structural": [
              {
                "pdb_id": "5JDR",
                "epitope_overlap_score": 0.71,
                "overlapping_residues": [52, 54, 55, 76, 77]
              }
            ]
          }
        }
      ]
    }
  ]
}
```

---

## 3. Data Sources and Matching Logic

### 3.1 Commercial Antibody Lookup (Option 2)

#### Data Sources

1. **Antibodypedia** (Primary)
   - **Access**: Web scraping or API (if available)
   - **Coverage**: Comprehensive commercial antibody database
   - **Data**: Vendor, clone, applications, validation status

2. **Vendor APIs/Web Scraping** (Secondary)
   - BioLegend, Cell Signaling Technology (CST), Abcam, R&D Systems
   - **Strategy**: Scrape product search pages or use vendor APIs
   - **Fallback**: Manual curation of common targets

3. **Local Database** (Fallback)
   - Pre-compiled JSON/CSV of common EV biomarkers and their antibodies
   - Updated periodically

#### Matching Logic

```python
def find_commercial_antibodies(
    gene_symbol: str,
    uniprot: Optional[str] = None,
    vendors: Optional[List[str]] = None
) -> List[CommercialAntibody]:
    """
    Find commercial antibodies for a target.
    
    Pseudocode:
    1. Query Antibodypedia by gene_symbol (primary) or UniProt (secondary)
    2. For each vendor in vendors list (or all):
       a. Query vendor search API/page
       b. Parse product listings
       c. Extract: clone, applications, epitope info, URL
    3. Deduplicate by clone_name + vendor
    4. Rank by validation status and application diversity
    5. Return top N results
    """
    results = []
    
    # Primary: Antibodypedia
    antibodypedia_results = query_antibodypedia(gene_symbol, uniprot)
    results.extend(antibodypedia_results)
    
    # Secondary: Vendor-specific queries
    if vendors:
        for vendor in vendors:
            vendor_results = query_vendor_database(vendor, gene_symbol, uniprot)
            results.extend(vendor_results)
    
    # Deduplicate and rank
    deduplicated = deduplicate_antibodies(results)
    ranked = rank_antibodies(deduplicated)
    
    return ranked[:max_results]
```

#### Implementation Details

**Antibodypedia Query**:
```python
def query_antibodypedia(gene_symbol: str, uniprot: Optional[str] = None) -> List[CommercialAntibody]:
    """
    Query Antibodypedia database.
    
    Strategy:
    - Use requests + BeautifulSoup for web scraping
    - Search URL: https://www.antibodypedia.com/search?q={gene_symbol}
    - Parse product table: vendor, clone, applications, validation
    - Extract product URLs for detailed info
    """
    import requests
    from bs4 import BeautifulSoup
    
    search_url = f"https://www.antibodypedia.com/search?q={gene_symbol}"
    response = requests.get(search_url, headers={"User-Agent": "Exosome-GPT/1.0"})
    soup = BeautifulSoup(response.content, 'html.parser')
    
    antibodies = []
    # Parse product listings from HTML
    for product_row in soup.find_all('tr', class_='product-row'):
        vendor = extract_vendor(product_row)
        clone = extract_clone(product_row)
        applications = extract_applications(product_row)
        # ... parse other fields
        
        ab = CommercialAntibody(
            vendor=vendor,
            clone_name=clone,
            target_gene=gene_symbol,
            applications=applications,
            # ...
        )
        antibodies.append(ab)
    
    return antibodies
```

**Vendor-Specific Queries**:
```python
VENDOR_APIS = {
    "BioLegend": {
        "search_url": "https://www.biolegend.com/en-us/search?q={query}",
        "parser": parse_biolegend_results
    },
    "CST": {
        "search_url": "https://www.cellsignal.com/products/primary-antibodies?N={query}",
        "parser": parse_cst_results
    },
    # ... other vendors
}
```

### 3.2 PDB Structural Matching (Option 3)

#### Data Sources

1. **PDB** (Primary)
   - **Access**: RCSB PDB API and bulk downloads
   - **Query**: Antibody-antigen complexes where antigen matches target

2. **SAbDab** (Future enhancement)
   - Pre-processed antibody structures
   - Faster lookup than raw PDB

#### Matching Logic

```python
def find_structural_antibodies(
    target_uniprot: str,
    target_structure_path: str,
    epitope_residues: List[int],
    overlap_threshold: float = 0.3
) -> List[StructuralAntibody]:
    """
    Find PDB antibody-antigen complexes with epitope overlap.
    
    Pseudocode:
    1. Query PDB for complexes containing target protein:
       - Search by UniProt ID in PDB metadata
       - Filter for antibody-antigen complexes (entity type)
    
    2. For each matching PDB entry:
       a. Download structure file
       b. Identify antigen chain and antibody chains (H, L)
       c. Compute antibody-antigen interface residues
       d. Map interface residues to target structure numbering
       e. Compute overlap with Step 2 epitope residues
    
    3. Filter by overlap_threshold (Jaccard index)
    
    4. Rank by overlap score and structure quality (resolution)
    
    5. Return top N results
    """
    # Step 1: Query PDB
    pdb_entries = query_pdb_by_uniprot(target_uniprot)
    antibody_complexes = filter_antibody_antigen_complexes(pdb_entries)
    
    matches = []
    for pdb_id in antibody_complexes:
        # Step 2: Download and analyze
        structure = download_pdb_structure(pdb_id)
        antigen_chain, antibody_chains = identify_chains(structure)
        
        # Compute interface
        interface_residues = compute_interface_residues(
            structure, antigen_chain, antibody_chains
        )
        
        # Map to target structure numbering
        mapped_interface = map_residue_numbers(
            interface_residues, 
            structure, 
            target_structure_path
        )
        
        # Compute overlap
        overlap_score = compute_epitope_overlap(
            epitope_residues, 
            mapped_interface
        )
        
        if overlap_score >= overlap_threshold:
            matches.append(StructuralAntibody(
                pdb_id=pdb_id,
                epitope_overlap_score=overlap_score,
                overlapping_residues=find_overlapping_residues(
                    epitope_residues, mapped_interface
                ),
                interface_residues=mapped_interface,
                # ...
            ))
    
    # Rank and return
    ranked = sorted(matches, key=lambda x: x.epitope_overlap_score, reverse=True)
    return ranked[:max_results]
```

#### Epitope Overlap Computation

```python
def compute_epitope_overlap(
    epitope_residues: List[int],
    interface_residues: List[int],
    distance_cutoff: float = 5.0  # Angstroms for spatial proximity
) -> float:
    """
    Compute overlap between Step 2 epitope and antibody interface.
    
    Methods:
    1. Residue-based Jaccard index (exact match)
    2. Spatial proximity (if structures available)
    
    Returns:
        Overlap score [0, 1]
    """
    # Method 1: Exact residue overlap (Jaccard)
    epitope_set = set(epitope_residues)
    interface_set = set(interface_residues)
    
    intersection = epitope_set & interface_set
    union = epitope_set | interface_set
    
    jaccard = len(intersection) / len(union) if union else 0.0
    
    # Method 2: Spatial proximity (if structures available)
    # Compute Cα distances between epitope and interface residues
    # Count residues within distance_cutoff
    # spatial_overlap = count_close_residues / len(epitope_residues)
    
    return jaccard  # or weighted combination
```

#### Residue Numbering Mapping

```python
def map_residue_numbers(
    source_residues: List[int],
    source_structure: PDB.Structure,
    target_structure_path: str
) -> List[int]:
    """
    Map residue numbers from PDB complex to target structure.
    
    Strategy:
    1. Align sequences (if available) or structures
    2. Use alignment to map residue indices
    3. Handle domain extraction (ectodomain vs full length)
    """
    from Bio import pairwise2
    from Bio.SeqUtils import seq1
    
    # Extract sequences
    source_seq = extract_sequence(source_structure, antigen_chain)
    target_seq = extract_sequence_from_pdb(target_structure_path)
    
    # Align
    alignment = pairwise2.align.globalxx(source_seq, target_seq)[0]
    
    # Map residues
    mapped = []
    for res in source_residues:
        mapped_res = map_via_alignment(res, alignment)
        if mapped_res:
            mapped.append(mapped_res)
    
    return mapped
```

---

## 4. GUI Design

### 4.1 Integration Point: Step 3 Panel Enhancement

**Location**: `gui/src/components/Step3Panel.tsx`

**Design**: Add a **mode selector** at the top of Step 3 panel:

```
┌─────────────────────────────────────────────────────────┐
│ Step 3 · Binder Optimization                            │
├─────────────────────────────────────────────────────────┤
│ Mode: ○ Design New Nanobinders  ● Find Existing Antibodies│
├─────────────────────────────────────────────────────────┤
│ [When "Find Existing Antibodies" selected]               │
│                                                          │
│ Step 2 JSON: [./epitopes/..._epitopes.json] [Browse...]  │
│                                                          │
│ Search Options:                                          │
│ ☑ Commercial Antibodies (BioLegend, CST, Abcam, etc.)   │
│ ☑ Structural Matches (PDB antibody-antigen complexes)    │
│                                                          │
│ Epitope Overlap Threshold: [0.3] (for structural)       │
│ Max Results per Epitope: [50]                            │
│                                                          │
│ [Run Antibody Search]                                    │
└─────────────────────────────────────────────────────────┘
```

### 4.2 Results Display Component

**New Component**: `gui/src/components/AntibodyFinderResults.tsx`

```tsx
interface AntibodyFinderResultsProps {
  results: AntibodyFinderResult;
  onAntibodySelect?: (antibody: CommercialAntibody | StructuralAntibody) => void;
}

export const AntibodyFinderResults: React.FC<AntibodyFinderResultsProps> = ({
  results,
  onAntibodySelect
}) => {
  return (
    <div className="antibody-results">
      {/* Target Summary */}
      <div className="target-summary">
        <h3>{results.target_name} ({results.gene_symbol})</h3>
        <p>Found {results.epitope_matches.length} epitope(s) with matches</p>
      </div>

      {/* Epitope Tabs */}
      <div className="epitope-tabs">
        {results.epitope_matches.map((match) => (
          <button key={match.epitope_id} className="epitope-tab">
            {match.epitope_id}
          </button>
        ))}
      </div>

      {/* Results Table */}
      <div className="results-table-container">
        <table className="antibody-table">
          <thead>
            <tr>
              <th>Type</th>
              <th>Vendor/PDB</th>
              <th>Clone/Name</th>
              <th>Applications</th>
              <th>Epitope Match</th>
              <th>Actions</th>
            </tr>
          </thead>
          <tbody>
            {/* Commercial Antibodies */}
            {match.commercial_antibodies.map((ab) => (
              <tr key={`${ab.vendor}-${ab.product_id}`}>
                <td><span className="badge badge-commercial">Commercial</span></td>
                <td>{ab.vendor}</td>
                <td>{ab.clone_name || "N/A"}</td>
                <td>
                  <div className="application-tags">
                    {ab.applications.map(app => (
                      <span key={app} className="tag">{app}</span>
                    ))}
                  </div>
                </td>
                <td>
                  {ab.reported_epitope ? (
                    <span className="epitope-info">{ab.reported_epitope}</span>
                  ) : (
                    <span className="text-muted">Not specified</span>
                  )}
                </td>
                <td>
                  <a href={ab.url} target="_blank" className="link-button">
                    View Product →
                  </a>
                </td>
              </tr>
            ))}

            {/* Structural Antibodies */}
            {match.structural_antibodies.map((ab) => (
              <tr key={ab.pdb_id}>
                <td><span className="badge badge-structural">Structural</span></td>
                <td>
                  <a 
                    href={`https://www.rcsb.org/structure/${ab.pdb_id}`}
                    target="_blank"
                  >
                    {ab.pdb_id}
                  </a>
                </td>
                <td>{ab.antibody_name || ab.clone_name || "N/A"}</td>
                <td>
                  <span className="text-muted">Structure Analysis</span>
                </td>
                <td>
                  <div className="overlap-score">
                    <span className="score-value">{ab.epitope_overlap_score.toFixed(2)}</span>
                    <span className="score-label">overlap</span>
                    <div className="overlap-residues">
                      {ab.overlapping_residues.length} residues
                    </div>
                  </div>
                </td>
                <td>
                  <button 
                    className="action-button"
                    onClick={() => onAntibodySelect?.(ab)}
                  >
                    View Structure →
                  </button>
                </td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>

      {/* Filters */}
      <div className="results-filters">
        <label>
          Filter by Vendor:
          <select>
            <option>All</option>
            <option>BioLegend</option>
            <option>CST</option>
            {/* ... */}
          </select>
        </label>
        <label>
          Filter by Application:
          <select>
            <option>All</option>
            <option>IHC</option>
            <option>WB</option>
            {/* ... */}
          </select>
        </label>
        <label>
          Min Overlap Score:
          <input type="range" min="0" max="1" step="0.1" />
        </label>
      </div>
    </div>
  );
};
```

### 4.3 User Flow

1. **User completes Step 2** → Gets curated targets with epitopes
2. **User navigates to Step 3** → Sees mode selector
3. **User selects "Find Existing Antibodies"** → Panel switches to antibody finder mode
4. **User selects Step 2 JSON** (auto-filled if available from Step 2)
5. **User configures search options**:
   - Toggle commercial vs structural search
   - Set overlap threshold
   - Select vendors (optional)
6. **User clicks "Run Antibody Search"** → Progress indicator shows
7. **Results display**:
   - Tabbed view by epitope
   - Table with filters
   - Links to vendor pages / PDB entries
8. **User can export results** → JSON or CSV

### 4.4 Visual Design Elements

**Badges**:
- `badge-commercial`: Blue background for commercial antibodies
- `badge-structural`: Green background for PDB structural matches
- `badge-high-confidence`: Gold for high overlap scores (>0.7)

**Tags**:
- Application tags (IHC, WB, FC, etc.) with color coding
- Epitope match quality indicators

**Tables**:
- Sortable columns
- Expandable rows for detailed info
- Hover effects for better UX

---

## 5. Backend API Integration

### 5.1 New API Endpoint

**File**: `backend/app.py`

```python
@app.route("/api/step3/antibody-finder", methods=["POST", "OPTIONS"])
def run_antibody_finder():
    """Find existing antibodies for Step 2 targets/epitopes."""
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
        search_mode = data.get("search_mode", "hybrid")
        commercial_sources = data.get("commercial_sources")
        epitope_overlap_threshold = data.get("epitope_overlap_threshold", 0.3)
        max_results_per_epitope = data.get("max_results_per_epitope", 50)
        
        from exo_gpt.step3_antibody_finder import find_existing_antibodies
        
        # Handle JSON data or file path
        if step2_json_data:
            import tempfile
            with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
                json.dump(step2_json_data, f)
                step2_json_path = f.name
        
        # Run antibody finder
        result = find_existing_antibodies(
            step2_json_path=step2_json_path,
            search_mode=search_mode,
            commercial_sources=commercial_sources,
            epitope_overlap_threshold=epitope_overlap_threshold,
            max_results_per_epitope=max_results_per_epitope
        )
        
        # Convert to dict for JSON response
        result_dict = asdict(result) if hasattr(result, '__dict__') else result
        
        return jsonify({
            "success": True,
            "result": result_dict
        }), 200
        
    except Exception as e:
        return jsonify({
            "success": False,
            "error": str(e),
            "traceback": traceback.format_exc()
        }), 500
```

### 5.2 Streaming Version (for long searches)

```python
@app.route("/api/step3/antibody-finder/stream", methods=["POST", "OPTIONS"])
def run_antibody_finder_stream():
    """Find existing antibodies with SSE streaming."""
    # Similar to run_step3_stream, but for antibody finder
    # Stream progress updates as databases are queried
    pass
```

---

## 6. Implementation Phases

### Phase 1: Commercial Antibody Lookup (Weeks 1-2)

**Deliverables**:
- `commercial_lookup.py` module
- Antibodypedia scraper
- Basic vendor query functions (2-3 vendors)
- JSON output format
- Basic GUI integration

**Testing**:
- Test with 5-10 common EV biomarkers (CD9, CD63, CD81, PD-L1, etc.)
- Validate data extraction accuracy

### Phase 2: PDB Structural Matching (Weeks 3-4)

**Deliverables**:
- `pdb_mining.py` module
- PDB API integration
- Interface residue computation
- Epitope overlap algorithm
- Residue numbering mapping

**Testing**:
- Test with known antibody-antigen complexes
- Validate overlap scores against manual inspection

### Phase 3: GUI Enhancement (Week 5)

**Deliverables**:
- Mode selector in Step 3 panel
- `AntibodyFinderResults` component
- Filters and sorting
- Export functionality

### Phase 4: Integration & Polish (Week 6)

**Deliverables**:
- Full pipeline integration
- Error handling
- Documentation
- Example workflows

---

## 7. Compatibility Considerations

### 7.1 JSON Structure Compatibility

- **Input**: Standard Step 2 JSON (no changes required)
- **Output**: Extended Step 2 JSON with `existing_antibodies` field (backward compatible)
- **Optional output**: Separate JSON file for antibody results only

### 7.2 Backend Compatibility

- Reuses existing Flask API structure
- Follows same error handling patterns
- Compatible with existing Step 3 orchestrator (doesn't interfere)

### 7.3 Execution Model

- **Local execution**: Can run on user's machine (no HPC required)
- **Database access**: Requires internet for commercial/PDB queries
- **Caching**: Local cache for frequently queried targets

---

## 8. Future Enhancements

1. **SAbDab Integration**: Direct SAbDab API for faster structural queries
2. **IEDB Integration**: Immune epitope database for validated epitopes
3. **Machine Learning**: Predict antibody-epitope compatibility
4. **Batch Processing**: Process multiple targets/epitopes in parallel
5. **Custom Database**: Allow users to add their own antibody databases

---

## 9. Example Usage

### Command Line

```bash
python -m exo_gpt.step3_antibody_finder \
  --step2_json "./epitopes/mm_plasma_epitopes.json" \
  --output_json "./antibody_results/mm_plasma_antibodies.json" \
  --search_mode "hybrid" \
  --epitope_overlap_threshold 0.3 \
  --max_results_per_epitope 50
```

### Python API

```python
from exo_gpt.step3_antibody_finder import find_existing_antibodies

result = find_existing_antibodies(
    step2_json_path="./epitopes/mm_plasma_epitopes.json",
    search_mode="hybrid",
    commercial_sources=["BioLegend", "CST"],
    epitope_overlap_threshold=0.3
)

# Access results
for target in result.curated_targets:
    for epitope_match in target.epitope_matches:
        print(f"Epitope {epitope_match.epitope_id}:")
        print(f"  Commercial: {len(epitope_match.commercial_antibodies)}")
        print(f"  Structural: {len(epitope_match.structural_antibodies)}")
```

---

## 10. Conclusion

This design provides a **practical, implementable solution** for finding existing antibodies that complements the *de novo* nanobinder design workflow. The hybrid approach (commercial + structural) maximizes both practical utility and scientific rigor, while maintaining full compatibility with the existing Exosome-GPT pipeline.

**Key Advantages**:
- ✅ No changes to Step 2 JSON structure
- ✅ Optional feature (doesn't break existing workflows)
- ✅ Practical for experimental planning
- ✅ Structurally validated when possible
- ✅ Clean GUI integration
- ✅ Extensible for future enhancements

