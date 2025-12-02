import React, { useState, useEffect } from "react";

interface PrimerPair {
  forward_primer: string;
  reverse_primer: string;
  forward_tm: number;
  reverse_tm: number;
  amplicon_size: number;
  amplicon_start: number;
  amplicon_end: number;
  gc_content: number;
  exon_junction?: string;
  isoform_specific: boolean;
  off_target_hits: number;
  primer3_score?: number;
}

interface PrimerDesign {
  gene_symbol: string;
  ensembl_id?: string;
  primer_pairs: PrimerPair[];
  strategy: string;
}

interface MrnaStep2Result {
  disease: string;
  biofluid: string;
  primer_designs: Record<string, PrimerDesign>;
  strategy: string;
}

export const MrnaStep2Panel: React.FC = () => {
  const [loading, setLoading] = useState(false);
  const [result, setResult] = useState<MrnaStep2Result | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [strategy, setStrategy] = useState<"dominant_isoform" | "isoform_specific">("dominant_isoform");

  // Restore last results
  useEffect(() => {
    try {
      const savedResult = localStorage.getItem("mrna_step2_result");
      if (savedResult) {
        const parsed: MrnaStep2Result = JSON.parse(savedResult);
        setResult(parsed);
        setStrategy(parsed.strategy as "dominant_isoform" | "isoform_specific");
      }
    } catch (e) {
      console.warn("Failed to restore mRNA Step 2 results from localStorage:", e);
    }
  }, []);

  const handleRun = async () => {
    setLoading(true);
    setError(null);
    setResult(null);

    try {
      // Get Step 1 results from localStorage
      const step1Result = localStorage.getItem("mrna_step1_result");
      if (!step1Result) {
        throw new Error("Please run mRNA Step 1 first to get biomarkers");
      }

      const step1Data = JSON.parse(step1Result);

      const response = await fetch("http://localhost:5000/api/mrna/step2/run", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
        body: JSON.stringify({
          step1_json_data: step1Data,
          strategy: strategy,
        }),
      });

      const data = await response.json();

      if (!data.success) {
        throw new Error(data.error || "Failed to design primers");
      }

      setResult(data.result);

      // Save to localStorage
      localStorage.setItem("mrna_step2_result", JSON.stringify(data.result));
    } catch (err) {
      setError(err instanceof Error ? err.message : "An error occurred");
    } finally {
      setLoading(false);
    }
  };

  return (
    <div className="panel">
      <div className="panel-header">
        <h2>mRNA Step 2: Primer Design</h2>
        <p>Automated RT-qPCR primer design using Primer3, Primer-BLAST, and Ensembl API</p>
      </div>

      <div className="panel-body">
        <div className="panel-grid">
          <div className="panel-card">
            <h3>Primer Design Parameters</h3>
            <div className="field">
              <label className="field-label">Design Strategy</label>
              <select
                value={strategy}
                onChange={(e) => setStrategy(e.target.value as "dominant_isoform" | "isoform_specific")}
                disabled={loading}
                style={{
                  borderRadius: "0.55rem",
                  border: "1px solid var(--border-subtle)",
                  background: "rgba(15, 23, 42, 0.9)",
                  color: "var(--text)",
                  padding: "6px 8px",
                  fontSize: "0.84rem",
                }}
              >
                <option value="dominant_isoform">Dominant Isoform Targeting</option>
                <option value="isoform_specific">Isoform-Specific Detection</option>
              </select>
              <div className="field-help">
                {strategy === "dominant_isoform"
                  ? "Targets the major/canonical transcript"
                  : "Targets isoform-unique exon-exon junctions"}
              </div>
            </div>
            <button
              className="run-button"
              onClick={handleRun}
              disabled={loading}
            >
              {loading ? "Designing Primers..." : "Run Primer Design"}
            </button>
            <div className="field-help" style={{ marginTop: "12px" }}>
              <strong>Note:</strong> Requires Step 1 results. Primers are designed for top biomarkers.
            </div>
          </div>

          <div className="panel-card">
            <h3>Primer Design Criteria</h3>
            <ul style={{ fontSize: "0.85rem", color: "var(--muted)", marginTop: "8px" }}>
              <li>Amplicon size: 70-150 bp</li>
              <li>Primer Tm: 58-62°C (forward-reverse diff ≤ 1°C)</li>
              <li>GC content: 40-60%</li>
              <li>3' GC clamp: 1-2 G/C</li>
              <li>Targets exon-exon junctions (avoids genomic DNA)</li>
              <li>Off-target screening via Primer-BLAST</li>
            </ul>
          </div>
        </div>

        {error && (
          <div className="error-message">
            <strong>Error:</strong> {error}
          </div>
        )}

        {result && (
          <div className="results-section">
            <h3>Primer Designs: {Object.keys(result.primer_designs).length} Genes</h3>
            
            {Object.entries(result.primer_designs).map(([geneSymbol, design]) => (
              <div key={geneSymbol} style={{ marginBottom: "24px", paddingBottom: "24px", borderBottom: "1px solid var(--border-subtle)" }}>
                <h4 style={{ marginBottom: "12px", color: "var(--accent)" }}>
                  {geneSymbol} {design.ensembl_id && `(${design.ensembl_id})`}
                </h4>
                
                {design.primer_pairs.map((pair, idx) => (
                  <div key={idx} style={{ marginBottom: "16px", padding: "12px", background: "rgba(15, 23, 42, 0.5)", borderRadius: "0.75rem" }}>
                    <div style={{ display: "grid", gridTemplateColumns: "1fr 1fr", gap: "12px", marginBottom: "8px" }}>
                      <div>
                        <strong style={{ fontSize: "0.85rem", color: "var(--muted)" }}>Forward Primer:</strong>
                        <div style={{ fontFamily: "monospace", fontSize: "0.9rem", marginTop: "4px" }}>
                          {pair.forward_primer}
                        </div>
                        <div style={{ fontSize: "0.8rem", color: "var(--muted)", marginTop: "2px" }}>
                          Tm: {pair.forward_tm.toFixed(1)}°C
                        </div>
                      </div>
                      <div>
                        <strong style={{ fontSize: "0.85rem", color: "var(--muted)" }}>Reverse Primer:</strong>
                        <div style={{ fontFamily: "monospace", fontSize: "0.9rem", marginTop: "4px" }}>
                          {pair.reverse_primer}
                        </div>
                        <div style={{ fontSize: "0.8rem", color: "var(--muted)", marginTop: "2px" }}>
                          Tm: {pair.reverse_tm.toFixed(1)}°C
                        </div>
                      </div>
                    </div>
                    
                    <div style={{ display: "flex", gap: "16px", fontSize: "0.85rem", color: "var(--muted)", flexWrap: "wrap" }}>
                      <span>Amplicon: {pair.amplicon_size} bp</span>
                      <span>GC: {pair.gc_content.toFixed(1)}%</span>
                      {pair.exon_junction && <span>Junction: {pair.exon_junction}</span>}
                      {pair.isoform_specific && <span style={{ color: "var(--accent)" }}>Isoform-specific</span>}
                      <span>Off-targets: {pair.off_target_hits}</span>
                      {pair.primer3_score !== null && pair.primer3_score !== undefined && (
                        <span>Score: {pair.primer3_score.toFixed(2)}</span>
                      )}
                    </div>
                  </div>
                ))}
              </div>
            ))}
          </div>
        )}
      </div>
    </div>
  );
};

