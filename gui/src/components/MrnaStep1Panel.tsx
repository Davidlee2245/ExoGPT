import React, { useState, useEffect, useRef } from "react";

interface MrnaBiomarker {
  gene_symbol: string;
  ensembl_id?: string;
  ev_logfc?: number;
  tcga_logfc?: number;
  geo_logfc?: number;
  ev_presence_score: number;
  ev_evidence_count: number;
  literature_hits: number;
  score: number;
  data_sources?: string[];
}

interface MrnaStep1Result {
  disease: string;
  biofluid: string;
  ev_mrna_biomarkers: MrnaBiomarker[];
  data_sources?: string[];
}

// Valid HPA API cancer search terms (categorized)
// Using HPA's exact terminology as required by their API
const HPA_CANCER_TYPES = {
  "Common Cancers": [
    "Colorectal cancer",
    "Breast cancer",
    "Lung cancer",
    "Prostate cancer",
    "Skin Cutaneous Melanoma",  // HPA exact term for melanoma
    "Pancreatic cancer",
    "Gastric cancer",
    "Liver cancer",
    "Ovarian cancer",
  ],
  "Other Cancers": [
    "Cervical cancer",
    "Endometrial cancer",
    "Renal cancer",
    "Bladder cancer",
    "Thyroid cancer",
    "Brain cancer",
  ],
};

// Flatten for autocomplete
const ALL_HPA_CANCER_TYPES = [
  ...HPA_CANCER_TYPES["Common Cancers"],
  ...HPA_CANCER_TYPES["Other Cancers"],
];

export const MrnaStep1Panel: React.FC = () => {
  const [disease, setDisease] = useState("Metastatic melanoma");
  const [biofluid, setBiofluid] = useState("Plasma");
  const [loading, setLoading] = useState(false);
  const [result, setResult] = useState<MrnaStep1Result | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [markdownTable, setMarkdownTable] = useState<string | null>(null);
  const [showHpaSuggestions, setShowHpaSuggestions] = useState(false);
  const [useExorbase2, setUseExorbase2] = useState(true);
  const [useTCGA, setUseTCGA] = useState(true);
  const diseaseInputRef = useRef<HTMLInputElement>(null);
  const suggestionsRef = useRef<HTMLDivElement>(null);

  // Restore last results
  useEffect(() => {
    try {
      const savedResult = localStorage.getItem("mrna_step1_result");
      const savedDisease = localStorage.getItem("mrna_step1_disease");
      const savedBiofluid = localStorage.getItem("mrna_step1_biofluid");
      const savedMarkdown = localStorage.getItem("mrna_step1_markdown_table");

      if (savedDisease) {
        setDisease(savedDisease);
      }
      if (savedBiofluid) {
        setBiofluid(savedBiofluid);
      }
      if (savedResult) {
        const parsed: MrnaStep1Result = JSON.parse(savedResult);
        setResult(parsed);
      }
      if (savedMarkdown) {
        setMarkdownTable(savedMarkdown || null);
      }
    } catch (e) {
      console.warn("Failed to restore mRNA Step 1 results from localStorage:", e);
    }
  }, []);

  const handleRun = async () => {
    setLoading(true);
    setError(null);
    setResult(null);
    setMarkdownTable(null);

    try {
      const response = await fetch("http://localhost:5000/api/mrna/step1/run", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
        body: JSON.stringify({
          disease,
          biofluid,
          top_n: 50,
          use_exorbase2: useExorbase2,
          use_tcga: useTCGA,
        }),
      });

      const data = await response.json();

      if (!data.success) {
        throw new Error(data.error || "Failed to run mRNA biomarker discovery");
      }

      setResult(data.result);
      setMarkdownTable(data.markdown_table || null);

      // Save to localStorage
      localStorage.setItem("mrna_step1_result", JSON.stringify(data.result));
      localStorage.setItem("mrna_step1_disease", disease);
      localStorage.setItem("mrna_step1_biofluid", biofluid);
      if (data.markdown_table) {
        localStorage.setItem("mrna_step1_markdown_table", data.markdown_table);
      }
    } catch (err) {
      setError(err instanceof Error ? err.message : "An error occurred");
    } finally {
      setLoading(false);
    }
  };

  return (
    <div className="panel">
      <div className="panel-header">
        <h2>mRNA Step 1: EV Biomarker Finder</h2>
        <p>Discover EV-derived mRNA biomarkers from ExoRBase, ExoCarta, Vesiclepedia, TCGA, GEO, and PubMed</p>
      </div>

      <div className="panel-body">
        <div className="panel-grid">
          <div className="panel-card">
            <h3>Input Parameters</h3>
            <div className="field" style={{ position: "relative" }}>
              <label className="field-label">
                Disease
                <span style={{ fontSize: "0.75rem", color: "var(--muted)", marginLeft: "6px" }}>
                  (HPA-compatible terms available)
                </span>
              </label>
              <input
                ref={diseaseInputRef}
                type="text"
                value={disease}
                onChange={(e) => {
                  setDisease(e.target.value);
                  setShowHpaSuggestions(e.target.value.length > 0);
                }}
                onFocus={() => setShowHpaSuggestions(true)}
                onBlur={() => {
                  // Delay hiding suggestions to allow clicking on them
                  setTimeout(() => setShowHpaSuggestions(false), 200);
                }}
                placeholder="e.g., Metastatic melanoma or select from HPA types"
                disabled={loading}
                style={{ position: "relative" }}
              />
              {showHpaSuggestions && (
                <div
                  ref={suggestionsRef}
                  className="hpa-suggestions"
                  style={{
                    position: "absolute",
                    top: "100%",
                    left: 0,
                    right: 0,
                    zIndex: 1000,
                    backgroundColor: "var(--bg-panel)",
                    border: "1px solid var(--border-subtle)",
                    borderRadius: "0.5rem",
                    marginTop: "4px",
                    maxHeight: "300px",
                    overflowY: "auto",
                    boxShadow: "0 4px 6px rgba(0, 0, 0, 0.1)",
                  }}
                >
                  <div
                    style={{
                      padding: "8px 12px",
                      fontSize: "0.75rem",
                      color: "var(--muted)",
                      fontWeight: "bold",
                      borderBottom: "1px solid var(--border-subtle)",
                    }}
                  >
                    HPA API Compatible Cancer Types:
                  </div>
                  {Object.entries(HPA_CANCER_TYPES).map(([category, types]) => (
                    <div key={category}>
                      <div
                        style={{
                          padding: "6px 12px",
                          fontSize: "0.75rem",
                          color: "var(--muted)",
                          fontWeight: 600,
                          backgroundColor: "rgba(100, 116, 139, 0.1)",
                        }}
                      >
                        {category}
                      </div>
                      {types
                        .filter(
                          (type) =>
                            disease.length === 0 ||
                            type.toLowerCase().includes(disease.toLowerCase())
                        )
                        .map((type) => (
                          <div
                            key={type}
                            onClick={() => {
                              setDisease(type);
                              setShowHpaSuggestions(false);
                              diseaseInputRef.current?.focus();
                            }}
                            style={{
                              padding: "8px 16px",
                              cursor: "pointer",
                              fontSize: "0.85rem",
                              borderBottom: "1px solid rgba(100, 116, 139, 0.1)",
                              transition: "background 0.15s",
                            }}
                            onMouseEnter={(e) => {
                              e.currentTarget.style.backgroundColor = "var(--accent-soft)";
                            }}
                            onMouseLeave={(e) => {
                              e.currentTarget.style.backgroundColor = "transparent";
                            }}
                          >
                            {type}
                          </div>
                        ))}
                    </div>
                  ))}
                  {disease.length > 0 &&
                    !ALL_HPA_CANCER_TYPES.some((type) =>
                      type.toLowerCase().includes(disease.toLowerCase())
                    ) && (
                      <div
                        style={{
                          padding: "8px 12px",
                          fontSize: "0.8rem",
                          color: "#f59e0b",
                          fontStyle: "italic",
                        }}
                      >
                        ⚠️ "{disease}" may not be HPA-compatible. Try selecting from suggestions above.
                      </div>
                    )}
                </div>
              )}
            </div>
            <div className="field">
              <label className="field-label">Biofluid</label>
              <input
                type="text"
                value={biofluid}
                onChange={(e) => setBiofluid(e.target.value)}
                placeholder="e.g., Plasma"
                disabled={loading}
              />
            </div>
            <div className="field">
              <label className="field-label">Data Sources</label>
              <div style={{ display: "flex", flexDirection: "column", gap: "10px", marginTop: "8px" }}>
                <label style={{ display: "flex", alignItems: "center", gap: "10px", cursor: "pointer", userSelect: "none" }}>
                  <input
                    type="checkbox"
                    checked={useExorbase2}
                    onChange={(e) => setUseExorbase2(e.target.checked)}
                    disabled={loading}
                    style={{ width: "18px", height: "18px", cursor: "pointer", accentColor: "var(--accent)" }}
                  />
                  <span style={{ fontSize: "0.9rem", color: "var(--text)" }}>ExoRBase2 (EV RNA-seq data)</span>
                </label>
                <label style={{ display: "flex", alignItems: "center", gap: "10px", cursor: "pointer", userSelect: "none" }}>
                  <input
                    type="checkbox"
                    checked={useTCGA}
                    onChange={(e) => setUseTCGA(e.target.checked)}
                    disabled={loading}
                    style={{ width: "18px", height: "18px", cursor: "pointer", accentColor: "var(--accent)" }}
                  />
                  <span style={{ fontSize: "0.9rem", color: "var(--text)" }}>TCGA (Tumor vs Normal expression)</span>
                </label>
              </div>
              <p className="field-help" style={{ marginTop: "6px", fontSize: "0.8rem" }}>
                Select which data sources to use for biomarker discovery. At least one must be selected.
              </p>
            </div>
            <button
              className="run-button"
              onClick={handleRun}
              disabled={loading || !disease || !biofluid || (!useExorbase2 && !useTCGA)}
            >
              {loading ? "Running..." : "Run mRNA Step 1"}
            </button>
          </div>

          <div className="panel-card">
            <h3>Data Sources</h3>
            <p style={{ fontSize: "0.85rem", color: "var(--muted)" }}>
              The pipeline queries:
            </p>
            <ul style={{ fontSize: "0.85rem", color: "var(--muted)", marginTop: "8px" }}>
              <li>ExoRBase/ExoRBase2: EV RNA-seq metrics</li>
              <li>ExoCarta/Vesiclepedia: EV mRNA evidence</li>
              <li>TCGA: Tumor vs normal mRNA expression</li>
              <li>GEO: Disease vs control RNA-seq</li>
              <li>PubMed: Literature evidence</li>
            </ul>
            {result && result.data_sources && result.data_sources.length > 0 && (
              <div style={{ marginTop: "12px" }}>
                <p style={{ fontSize: "0.85rem", fontWeight: 600 }}>Sources used:</p>
                <p style={{ fontSize: "0.8rem", color: "var(--muted)" }}>
                  {result.data_sources.join(", ")}
                </p>
              </div>
            )}
          </div>
        </div>

        {error && (
          <div className="error-message">
            <strong>Error:</strong> {error}
          </div>
        )}

        {result && (
          <div className="results-section">
            <h3>Results: {result.ev_mrna_biomarkers.length} mRNA Biomarkers Found</h3>
            
            {result.ev_mrna_biomarkers.length > 0 && (
              <div className="biomarker-table-container">
                <table className="biomarker-table">
                  <thead>
                    <tr>
                      <th>Rank</th>
                      <th>Gene</th>
                      <th>Ensembl ID</th>
                      <th>EV logFC</th>
                      <th>Disease logFC</th>
                      <th>EV Presence</th>
                      <th>Literature</th>
                      <th>Score</th>
                    </tr>
                  </thead>
                  <tbody>
                    {result.ev_mrna_biomarkers.slice(0, 50).map((bm, idx) => (
                      <tr key={idx}>
                        <td>{idx + 1}</td>
                        <td>{bm.gene_symbol}</td>
                        <td>{bm.ensembl_id || "N/A"}</td>
                        <td>{bm.ev_logfc !== null && bm.ev_logfc !== undefined ? bm.ev_logfc.toFixed(2) : "N/A"}</td>
                        <td>
                          {bm.tcga_logfc !== null && bm.tcga_logfc !== undefined
                            ? bm.tcga_logfc.toFixed(2)
                            : bm.geo_logfc !== null && bm.geo_logfc !== undefined
                            ? bm.geo_logfc.toFixed(2)
                            : "N/A"}
                        </td>
                        <td>{bm.ev_presence_score}</td>
                        <td>{bm.literature_hits}</td>
                        <td>{bm.score.toFixed(2)}</td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            )}

            {markdownTable && (
              <div className="markdown-output">
                <h4>Markdown Output</h4>
                <pre className="markdown-block">{markdownTable}</pre>
              </div>
            )}
          </div>
        )}
      </div>
    </div>
  );
};

