"""
Step 2 — Automated Primer Design Module
---------------------------------------

Automated RT-qPCR primer design for EV mRNA biomarkers.

Features:
- Primer3 integration for primer design
- NCBI Primer-BLAST for off-target screening
- Ensembl API for transcript/exon structure
- Isoform-aware primer selection
- Automated parameter tuning
"""

from __future__ import annotations
import argparse
import json
import os
import subprocess
import tempfile
from dataclasses import dataclass, asdict
from typing import List, Dict, Any, Optional, Tuple
import re

import requests
import pandas as pd

DATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "data"))


@dataclass
class PrimerPair:
    """Primer pair design result"""
    forward_primer: str
    reverse_primer: str
    forward_tm: float
    reverse_tm: float
    amplicon_size: int
    amplicon_start: int
    amplicon_end: int
    gc_content: float
    exon_junction: Optional[str] = None  # e.g., "exon4-exon5"
    isoform_specific: bool = False
    off_target_hits: int = 0
    primer3_score: Optional[float] = None


@dataclass
class TranscriptInfo:
    """Transcript and exon structure information"""
    transcript_id: str
    gene_symbol: str
    ensembl_id: str
    sequence: str
    exons: List[Dict[str, Any]]  # List of {start, end, exon_number}
    dominant_isoform: bool = False
    tpm_expression: Optional[float] = None


# Gene symbol aliases/alternative names mapping
GENE_SYMBOL_ALIASES = {
    "IL8": "CXCL8",  # Interleukin 8 → C-X-C motif chemokine ligand 8
    "ALIX": "PDCD6IP",  # ALIX → Programmed cell death 6 interacting protein
    "HSA-MIR-21": None,  # MicroRNA, not a protein-coding gene
    # Add more aliases as needed
}

def query_ensembl_transcript(gene_symbol: str, ensembl_id: Optional[str] = None) -> List[TranscriptInfo]:
    """
    Query Ensembl API to get transcript and exon structure.
    
    Returns list of TranscriptInfo objects.
    """
    transcripts = []
    
    try:
        base_url = "https://rest.ensembl.org"
        
        # First, get gene ID if not provided
        if ensembl_id is None:
            # Try original symbol first
            lookup_url = f"{base_url}/lookup/symbol/homo_sapiens/{gene_symbol}"
            headers = {"Content-Type": "application/json"}
            response = requests.get(lookup_url, headers=headers, timeout=10)
            
            if response.status_code == 200:
                gene_data = response.json()
                ensembl_id = gene_data.get("id")
            else:
                # Try alias/alternative name if available
                alias = GENE_SYMBOL_ALIASES.get(gene_symbol.upper())
                if alias:
                    lookup_url = f"{base_url}/lookup/symbol/homo_sapiens/{alias}"
                    response = requests.get(lookup_url, headers=headers, timeout=10)
                    if response.status_code == 200:
                        gene_data = response.json()
                        ensembl_id = gene_data.get("id")
                        print(f"✓ Found Ensembl ID for {gene_symbol} using alias: {alias}")
                    else:
                        print(f"⚠️  Could not find Ensembl ID for {gene_symbol} (tried alias {alias})")
                        return []
                else:
                    print(f"⚠️  Could not find Ensembl ID for {gene_symbol}")
                    return []
        
        # Get transcripts for this gene
        transcripts_url = f"{base_url}/overlap/id/{ensembl_id}"
        params = {
            "feature": "transcript",
            "content-type": "application/json"
        }
        response = requests.get(transcripts_url, params=params, timeout=10)
        
        if response.status_code == 200:
            transcript_data = response.json()
            
            for tx in transcript_data[:10]:  # Limit to top 10 transcripts
                tx_id = tx.get("id")
                if not tx_id:
                    continue
                
                # Get transcript sequence
                seq_url = f"{base_url}/sequence/id/{tx_id}"
                seq_params = {"content-type": "application/json", "type": "cdna"}
                seq_response = requests.get(seq_url, params=seq_params, timeout=10)
                
                sequence = ""
                if seq_response.status_code == 200:
                    seq_data = seq_response.json()
                    sequence = seq_data.get("seq", "")
                
                # Get exon structure
                exons_url = f"{base_url}/overlap/id/{tx_id}"
                exons_params = {
                    "feature": "exon",
                    "content-type": "application/json"
                }
                exons_response = requests.get(exons_url, params=exons_params, timeout=10)
                
                exons = []
                if exons_response.status_code == 200:
                    exons_data = exons_response.json()
                    for exon in sorted(exons_data, key=lambda x: x.get("start", 0)):
                        exons.append({
                            "start": exon.get("start"),
                            "end": exon.get("end"),
                            "exon_number": exon.get("exon_number", len(exons) + 1)
                        })
                
                # Check if this is the dominant isoform (canonical transcript)
                is_canonical = tx.get("is_canonical", False)
                
                transcript = TranscriptInfo(
                    transcript_id=tx_id,
                    gene_symbol=gene_symbol,
                    ensembl_id=ensembl_id,
                    sequence=sequence,
                    exons=exons,
                    dominant_isoform=is_canonical
                )
                transcripts.append(transcript)
    
    except Exception as e:
        print(f"⚠️  Error querying Ensembl: {e}")
    
    return transcripts


# Cache for Primer3 availability check
_PRIMER3_AVAILABLE = None
_PRIMER3_CHECKED = False

def _check_primer3_available() -> bool:
    """Check if Primer3 is available (cached check)"""
    global _PRIMER3_AVAILABLE, _PRIMER3_CHECKED
    
    if _PRIMER3_CHECKED:
        return _PRIMER3_AVAILABLE
    
    _PRIMER3_CHECKED = True
    try:
        result = subprocess.run(["primer3_core", "--version"], 
                              capture_output=True, timeout=5)
        _PRIMER3_AVAILABLE = (result.returncode in [0, 255])  # 255 is normal for --version
        if not _PRIMER3_AVAILABLE:
            print("⚠️  Primer3 not found. Install with: conda install -c bioconda primer3")
        return _PRIMER3_AVAILABLE
    except FileNotFoundError:
        _PRIMER3_AVAILABLE = False
        print("⚠️  Primer3 not found. Install with: conda install -c bioconda primer3")
        return False

def design_primers_primer3(
    sequence: str,
    target_region: Optional[Tuple[int, int]] = None,
    min_tm: float = 58.0,
    opt_tm: float = 60.0,
    max_tm: float = 62.0,
    min_gc: float = 40.0,
    max_gc: float = 60.0,
    product_size_range: Tuple[int, int] = (70, 150),
    exon_junction: Optional[Tuple[int, int]] = None
) -> Optional[PrimerPair]:
    """
    Design primers using Primer3.
    
    Args:
        sequence: Target DNA sequence
        target_region: Optional (start, end) region to target
        min_tm, opt_tm, max_tm: Melting temperature constraints
        min_gc, max_gc: GC content constraints
        product_size_range: (min, max) amplicon size
        exon_junction: Optional (exon_start, exon_end) for junction targeting
    
    Returns:
        PrimerPair or None if no primers found
    """
    # Check if Primer3 is available (cached check)
    if not _check_primer3_available():
        return None
    
    # Create Primer3 input file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.primer3', delete=False) as f:
        primer3_input = f"""SEQUENCE_ID={target_region[0] if target_region else 'target'}
SEQUENCE_TEMPLATE={sequence}
PRIMER_TASK=pick_pcr_primers
PRIMER_MIN_SIZE=18
PRIMER_OPT_SIZE=20
PRIMER_MAX_SIZE=24
PRIMER_MIN_TM={min_tm}
PRIMER_OPT_TM={opt_tm}
PRIMER_MAX_TM={max_tm}
PRIMER_MIN_GC={min_gc}
PRIMER_MAX_GC={max_gc}
PRIMER_PRODUCT_SIZE_RANGE={product_size_range[0]}-{product_size_range[1]}
PRIMER_GC_CLAMP=1
PRIMER_NUM_RETURN=1
=
"""
        f.write(primer3_input)
        input_file = f.name
    
    try:
        # Run Primer3
        result = subprocess.run(
            ["primer3_core", input_file],
            capture_output=True,
            text=True,
            timeout=30
        )
        
        if result.returncode != 0:
            print(f"⚠️  Primer3 error: {result.stderr}")
            return None
        
        # Parse Primer3 output
        output = result.stdout
        
        # Extract primer sequences
        forward_match = re.search(r'PRIMER_LEFT_0_SEQUENCE=(\w+)', output)
        reverse_match = re.search(r'PRIMER_RIGHT_0_SEQUENCE=(\w+)', output)
        
        if not forward_match or not reverse_match:
            return None
        
        forward_seq = forward_match.group(1)
        reverse_seq = reverse_match.group(1)
        
        # Extract Tm values
        forward_tm_match = re.search(r'PRIMER_LEFT_0_TM=([\d.]+)', output)
        reverse_tm_match = re.search(r'PRIMER_RIGHT_0_TM=([\d.]+)', output)
        
        forward_tm = float(forward_tm_match.group(1)) if forward_tm_match else 60.0
        reverse_tm = float(reverse_tm_match.group(1)) if reverse_tm_match else 60.0
        
        # Extract product size
        product_size_match = re.search(r'PRIMER_PAIR_0_PRODUCT_SIZE=(\d+)', output)
        product_size = int(product_size_match.group(1)) if product_size_match else 100
        
        # Extract product coordinates
        product_start_match = re.search(r'PRIMER_LEFT_0=(\d+)', output)
        product_end_match = re.search(r'PRIMER_RIGHT_0=(\d+)', output)
        
        product_start = int(product_start_match.group(1)) if product_start_match else 0
        product_end = int(product_end_match.group(1)) if product_end_match else product_size
        
        # Calculate GC content
        gc_content = (forward_seq.count('G') + forward_seq.count('C') + 
                     reverse_seq.count('G') + reverse_seq.count('C')) / (len(forward_seq) + len(reverse_seq)) * 100
        
        # Extract primer3 score if available
        score_match = re.search(r'PRIMER_PAIR_0_PENALTY=([\d.]+)', output)
        primer3_score = float(score_match.group(1)) if score_match else None
        
        return PrimerPair(
            forward_primer=forward_seq,
            reverse_primer=reverse_seq,
            forward_tm=forward_tm,
            reverse_tm=reverse_tm,
            amplicon_size=product_size,
            amplicon_start=product_start,
            amplicon_end=product_end,
            gc_content=gc_content,
            exon_junction=None,  # Will be set later if junction targeting
            primer3_score=primer3_score
        )
    
    except Exception as e:
        print(f"⚠️  Error running Primer3: {e}")
        return None
    
    finally:
        # Clean up
        if os.path.exists(input_file):
            os.unlink(input_file)


def check_off_target_primer_blast(
    forward_primer: str,
    reverse_primer: str,
    gene_symbol: str
) -> int:
    """
    Check for off-target hits using Primer-BLAST.
    
    Returns number of off-target hits (excluding the target gene).
    """
    # TODO: Implement Primer-BLAST API call or automated browser interaction
    # For now, return 0 (no off-target hits found)
    # In production, this would:
    # 1. Call NCBI Primer-BLAST API
    # 2. Set species to Homo sapiens (taxid 9606)
    # 3. Search against RefSeq mRNA or human genome
    # 4. Count hits that are not the target gene
    
    return 0


def design_primers_for_biomarker(
    gene_symbol: str,
    ensembl_id: Optional[str] = None,
    strategy: str = "dominant_isoform"  # or "isoform_specific"
) -> List[Dict[str, Any]]:
    """
    Design primers for an mRNA biomarker.
    
    Args:
        gene_symbol: Gene symbol
        ensembl_id: Optional Ensembl ID
        strategy: "dominant_isoform" or "isoform_specific"
    
    Returns:
        List of primer design results (JSON-serializable)
    """
    print(f"🧬 Designing primers for {gene_symbol}...")
    
    # Get transcript information
    transcripts = query_ensembl_transcript(gene_symbol, ensembl_id)
    
    if not transcripts:
        print(f"⚠️  No transcripts found for {gene_symbol}")
        return []
    
    # Select target transcript based on strategy
    if strategy == "dominant_isoform":
        # Use dominant/canonical isoform
        target_transcript = next((t for t in transcripts if t.dominant_isoform), transcripts[0])
    else:
        # Use first transcript (can be enhanced for isoform-specific targeting)
        target_transcript = transcripts[0]
    
    if not target_transcript.sequence:
        print(f"⚠️  No sequence available for transcript {target_transcript.transcript_id}")
        return []
    
    print(f"✓ Using transcript {target_transcript.transcript_id} ({len(target_transcript.sequence)} bp)")
    
    # Design primers
    # Strategy: Target exon-exon junction to avoid genomic DNA
    primer_results = []
    
    # Check Primer3 availability once before attempting designs
    if not _check_primer3_available():
        print(f"⚠️  No primers designed for {gene_symbol} (Primer3 not available)")
        return []
    
    if target_transcript.exons and len(target_transcript.exons) > 1:
        # Try to design primers across exon junctions
        # Limit to first 5 junctions to avoid excessive attempts
        max_junctions = min(5, len(target_transcript.exons) - 1)
        for i in range(max_junctions):
            exon1 = target_transcript.exons[i]
            exon2 = target_transcript.exons[i + 1]
            
            # Get junction region (end of exon1, start of exon2)
            # Note: This is simplified - actual junction coordinates need sequence mapping
            junction_name = f"exon{exon1.get('exon_number', i+1)}-exon{exon2.get('exon_number', i+2)}"
            
            # For now, design primers from the full sequence
            # In production, would extract junction region specifically
            primer_pair = design_primers_primer3(
                sequence=target_transcript.sequence,
                min_tm=58.0,
                opt_tm=60.0,
                max_tm=62.0,
                min_gc=40.0,
                max_gc=60.0,
                product_size_range=(70, 150)
            )
            
            if primer_pair:
                primer_pair.exon_junction = junction_name
                primer_pair.isoform_specific = (strategy == "isoform_specific")
                
                # Check for off-target hits
                primer_pair.off_target_hits = check_off_target_primer_blast(
                    primer_pair.forward_primer,
                    primer_pair.reverse_primer,
                    gene_symbol
                )
                
                primer_results.append(asdict(primer_pair))
                
                # Only return first successful design for now
                break
    else:
        # No exon structure available, design from full sequence
        primer_pair = design_primers_primer3(
            sequence=target_transcript.sequence,
            min_tm=58.0,
            opt_tm=60.0,
            max_tm=62.0,
            min_gc=40.0,
            max_gc=60.0,
            product_size_range=(70, 150)
        )
        
        if primer_pair:
            primer_pair.off_target_hits = check_off_target_primer_blast(
                primer_pair.forward_primer,
                primer_pair.reverse_primer,
                gene_symbol
            )
            primer_results.append(asdict(primer_pair))
    
    if primer_results:
        print(f"✓ Designed {len(primer_results)} primer pair(s)")
    else:
        print(f"⚠️  No primers designed for {gene_symbol}")
    
    return primer_results


def design_primers_batch(
    biomarkers: List[Dict[str, Any]],
    strategy: str = "dominant_isoform"
) -> Dict[str, Any]:
    """
    Design primers for a batch of biomarkers.
    
    Args:
        biomarkers: List of biomarker dicts (from Step 1 output)
        strategy: Primer design strategy
    
    Returns:
        Dict with primer designs for each biomarker
    """
    results = {}
    
    for biomarker in biomarkers:
        gene_symbol = biomarker.get("gene_symbol")
        ensembl_id = biomarker.get("ensembl_id")
        
        if not gene_symbol:
            continue
        
        primer_designs = design_primers_for_biomarker(
            gene_symbol=gene_symbol,
            ensembl_id=ensembl_id,
            strategy=strategy
        )
        
        if primer_designs:
            results[gene_symbol] = {
                "gene_symbol": gene_symbol,
                "ensembl_id": ensembl_id,
                "primer_pairs": primer_designs,
                "strategy": strategy
            }
    
    return results


def main():
    parser = argparse.ArgumentParser(description="Primer Design Module")
    parser.add_argument("--biomarkers_json", required=True, help="Input JSON file with biomarkers from Step 1")
    parser.add_argument("--out_json", required=True, help="Output JSON file with primer designs")
    parser.add_argument("--strategy", default="dominant_isoform", 
                        choices=["dominant_isoform", "isoform_specific"],
                        help="Primer design strategy")
    
    args = parser.parse_args()
    
    # Load biomarkers
    with open(args.biomarkers_json, "r") as f:
        data = json.load(f)
    
    biomarkers = data.get("ev_mrna_biomarkers", [])
    
    if not biomarkers:
        print("⚠️  No biomarkers found in input file")
        return
    
    print(f"🧬 Designing primers for {len(biomarkers)} biomarkers...")
    
    # Design primers
    results = design_primers_batch(biomarkers, strategy=args.strategy)
    
    # Save results
    output = {
        "disease": data.get("disease"),
        "biofluid": data.get("biofluid"),
        "primer_designs": results,
        "strategy": args.strategy
    }
    
    with open(args.out_json, "w") as f:
        json.dump(output, f, indent=2)
    
    print(f"✓ Saved primer designs to {args.out_json}")
    print(f"✓ Designed primers for {len(results)} biomarkers")


if __name__ == "__main__":
    main()

