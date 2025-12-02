import React, { useState, useEffect } from "react";

export const Step3AntibodyPanel: React.FC = () => {
  const [step2JsonPath, setStep2JsonPath] = useState("./epitopes/mm_plasma_epitopes.json");
  const [searchMode, setSearchMode] = useState<"commercial" | "structural" | "hybrid">("hybrid");
  const [commercialSources, setCommercialSources] = useState<string[]>(["all"]);
  const [epitopeOverlapThreshold, setEpitopeOverlapThreshold] = useState(0.3);
  const [maxResultsPerEpitope, setMaxResultsPerEpitope] = useState(50);
  const [loading, setLoading] = useState(false);
  const [results, setResults] = useState<any>(null);
  const [error, setError] = useState<string | null>(null);
  const [step2JsonData, setStep2JsonData] = useState<any>(null);

  // Load Step 2 results from localStorage if available
  useEffect(() => {
    try {
      const savedResult = localStorage.getItem('step2_result');
      const savedPath = localStorage.getItem('step2_json_path');
      
      if (savedResult) {
        const result = JSON.parse(savedResult);
        setStep2JsonData(result);
        if (savedPath) {
          setStep2JsonPath(savedPath);
        }
      }
    } catch (e) {
      console.warn('Failed to load Step 2 results from localStorage:', e);
    }
  }, []);

  const handleRun = async () => {
    setLoading(true);
    setError(null);
    setResults(null);

    try {
      const requestBody: any = {
        search_mode: searchMode,
        commercial_sources: commercialSources,
        epitope_overlap_threshold: epitopeOverlapThreshold,
        max_results_per_epitope: maxResultsPerEpitope,
      };

      // If step2JsonData is provided, use it; otherwise use path
      if (step2JsonData) {
        requestBody.step2_json_data = step2JsonData;
      } else {
        requestBody.step2_json_path = step2JsonPath;
      }

      const response = await fetch("http://localhost:5000/api/step3/antibody-finder", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
        body: JSON.stringify(requestBody),
      });

      if (!response.ok) {
        const errorData = await response.json();
        throw new Error(errorData.error || "Failed to run antibody finder");
      }

      const data = await response.json();
      if (data.success) {
        setResults(data.result);
      } else {
        throw new Error(data.error || "Unknown error occurred");
      }
    } catch (err) {
      setError(err instanceof Error ? err.message : "Failed to run antibody finder");
    } finally {
      setLoading(false);
    }
  };

  return (
    <section className="panel">
      <header className="panel-header">
        <h2>Step 3-1 · Antibody Recommendation</h2>
        <p>
          Find existing commercial antibodies and structural antibody-antigen complexes that match your curated targets and epitopes from Step 2.
        </p>
      </header>
      <div className="panel-body">
        <div className="panel-grid">
          <div className="panel-card">
            <h3>Inputs</h3>
            <label className="field">
              <span className="field-label">Step 2 JSON Result Path</span>
              <input 
                value={step2JsonPath} 
                onChange={(e) => {
                  setStep2JsonPath(e.target.value);
                  setStep2JsonData(null);
                }}
                placeholder="./epitopes/..._epitopes.json"
                disabled={loading || !!step2JsonData}
              />
            </label>
            {step2JsonData && (
              <p className="field-help" style={{ color: "green" }}>
                ✓ Using Step 2 results directly (from Step 2 panel)
              </p>
            )}
            <label className="field">
              <span className="field-label">Search Mode</span>
              <select
                value={searchMode}
                onChange={(e) => setSearchMode(e.target.value as "commercial" | "structural" | "hybrid")}
                disabled={loading}
                style={{
                  borderRadius: "0.55rem",
                  border: "1px solid var(--border-subtle)",
                  background: "rgba(255, 255, 255, 0.8)",
                  color: "var(--text)",
                  padding: "6px 8px",
                  fontSize: "0.84rem",
                  width: "100%",
                }}
              >
                <option value="commercial">Commercial Antibodies Only</option>
                <option value="structural">Structural Matches Only (PDB)</option>
                <option value="hybrid">Hybrid (Commercial + Structural)</option>
              </select>
            </label>
            <label className="field">
              <span className="field-label">Epitope Overlap Threshold</span>
              <input 
                type="number"
                step="0.1"
                min="0"
                max="1"
                value={epitopeOverlapThreshold} 
                onChange={(e) => setEpitopeOverlapThreshold(parseFloat(e.target.value) || 0.3)}
                disabled={loading || searchMode === "commercial"}
              />
              <div className="field-help">
                Minimum Jaccard index for structural epitope matches (0.0 - 1.0). Only applies to structural/hybrid modes.
              </div>
            </label>
            <label className="field">
              <span className="field-label">Max Results per Epitope</span>
              <input 
                type="number"
                value={maxResultsPerEpitope} 
                onChange={(e) => setMaxResultsPerEpitope(parseInt(e.target.value) || 50)}
                min={1}
                max={200}
                disabled={loading}
              />
            </label>
            <button 
              className="run-button" 
              onClick={handleRun}
              disabled={loading || (!step2JsonPath && !step2JsonData)}
            >
              {loading ? "Searching..." : "Find Existing Antibodies"}
            </button>
          </div>
          <div className="panel-card">
            <h3>Results</h3>
            {error && (
              <div className="error-message">
                <strong>Error:</strong> {error}
              </div>
            )}
            {loading && (
              <div className="field-help">
                Searching for existing antibodies... This may take a few moments.
              </div>
            )}
            {results && (
              <div>
                <div className="result-summary">
                  <p><strong>Target:</strong> {results.target_name} ({results.gene_symbol})</p>
                  <p><strong>Epitopes with matches:</strong> {results.epitope_matches?.length || 0}</p>
                </div>
                {results.epitope_matches?.map((match: any, idx: number) => (
                  <div key={idx} style={{ marginTop: "20px" }}>
                    <h4 style={{ marginTop: 0, marginBottom: "12px" }}>
                      Epitope: {match.epitope_id}
                      <span style={{ marginLeft: "12px", fontSize: "0.85rem", fontWeight: "normal", color: "var(--muted)" }}>
                        ({match.match_confidence || "N/A"} confidence)
                      </span>
                    </h4>
                    
                    {match.commercial_antibodies && match.commercial_antibodies.length > 0 && (
                      <div style={{ marginTop: "12px" }}>
                        <h5 style={{ marginBottom: "8px", fontSize: "0.9rem", color: "var(--muted)" }}>
                          Commercial Antibodies ({match.commercial_antibodies.length})
                        </h5>
                        <div className="biomarker-table-container">
                          <table className="biomarker-table">
                            <thead>
                              <tr>
                                <th>Provider</th>
                                <th>Catalog Number</th>
                                <th>Type</th>
                                <th>Clone</th>
                                <th>Applications</th>
                                <th>Actions</th>
                              </tr>
                            </thead>
                            <tbody>
                              {match.commercial_antibodies.map((ab: any, abIdx: number) => (
                                <tr key={abIdx}>
                                  <td><strong>{ab.vendor || "Unknown"}</strong></td>
                                  <td>
                                    {ab.product_id && ab.product_id !== "N/A" && ab.product_id !== "multiple" ? (
                                      <code style={{ fontSize: "0.85rem", background: "rgba(20, 184, 166, 0.1)", padding: "2px 6px", borderRadius: "4px" }}>
                                        {ab.product_id}
                                      </code>
                                    ) : (
                                      <span style={{ color: "var(--muted)", fontStyle: "italic" }}>N/A</span>
                                    )}
                                  </td>
                                  <td>
                                    {ab.isotype ? (
                                      <span style={{ 
                                        fontSize: "0.85rem",
                                        padding: "2px 8px",
                                        borderRadius: "4px",
                                        background: ab.isotype.toLowerCase().includes("mono") 
                                          ? "rgba(20, 184, 166, 0.15)" 
                                          : "rgba(148, 163, 184, 0.15)",
                                        color: ab.isotype.toLowerCase().includes("mono")
                                          ? "var(--accent-strong)"
                                          : "var(--muted)"
                                      }}>
                                        {ab.isotype}
                                      </span>
                                    ) : (
                                      <span style={{ color: "var(--muted)", fontStyle: "italic" }}>—</span>
                                    )}
                                  </td>
                                  <td>
                                    {ab.clone_name && ab.clone_name !== "N/A" ? (
                                      <span style={{ fontSize: "0.85rem" }}>{ab.clone_name}</span>
                                    ) : (
                                      <span style={{ color: "var(--muted)", fontStyle: "italic" }}>—</span>
                                    )}
                                  </td>
                                  <td>
                                    {ab.applications && ab.applications.length > 0 ? (
                                      <div style={{ display: "flex", flexWrap: "wrap", gap: "4px" }}>
                                        {ab.applications.map((app: string, appIdx: number) => (
                                          <span
                                            key={appIdx}
                                            style={{
                                              fontSize: "0.75rem",
                                              padding: "2px 6px",
                                              borderRadius: "4px",
                                              background: "rgba(20, 184, 166, 0.15)",
                                              color: "var(--accent-strong)",
                                              fontWeight: "500"
                                            }}
                                          >
                                            {app}
                                          </span>
                                        ))}
                                      </div>
                                    ) : (
                                      <span style={{ color: "var(--muted)", fontStyle: "italic" }}>—</span>
                                    )}
                                  </td>
                                  <td>
                                    {ab.url && ab.url !== "N/A" && !ab.url.includes("antibodypedia.com/gene") ? (
                                      <a
                                        href={ab.url}
                                        target="_blank"
                                        rel="noopener noreferrer"
                                        style={{
                                          fontSize: "0.85rem",
                                          color: "var(--accent)",
                                          textDecoration: "none",
                                          fontWeight: "500"
                                        }}
                                        onMouseEnter={(e) => e.currentTarget.style.textDecoration = "underline"}
                                        onMouseLeave={(e) => e.currentTarget.style.textDecoration = "none"}
                                      >
                                        View Product →
                                      </a>
                                    ) : (
                                      <span style={{ color: "var(--muted)", fontSize: "0.85rem" }}>—</span>
                                    )}
                                  </td>
                                </tr>
                              ))}
                            </tbody>
                          </table>
                        </div>
                      </div>
                    )}
                    
                    {match.structural_antibodies && match.structural_antibodies.length > 0 && (
                      <div style={{ marginTop: "16px" }}>
                        <h5 style={{ marginBottom: "8px", fontSize: "0.9rem", color: "var(--muted)" }}>
                          Structural Matches ({match.structural_antibodies.length})
                        </h5>
                        <div className="biomarker-table-container">
                          <table className="biomarker-table">
                            <thead>
                              <tr>
                                <th>PDB ID</th>
                                <th>Antibody Name</th>
                                <th>Epitope Overlap</th>
                                <th>Resolution</th>
                                <th>Actions</th>
                              </tr>
                            </thead>
                            <tbody>
                              {match.structural_antibodies.map((ab: any, abIdx: number) => (
                                <tr key={abIdx}>
                                  <td><strong>{ab.pdb_id}</strong></td>
                                  <td>{ab.antibody_name || ab.clone_name || "—"}</td>
                                  <td>
                                    <span style={{ 
                                      fontSize: "0.85rem",
                                      padding: "2px 8px",
                                      borderRadius: "4px",
                                      background: ab.epitope_overlap_score > 0.5 
                                        ? "rgba(34, 197, 94, 0.15)" 
                                        : ab.epitope_overlap_score > 0.3
                                        ? "rgba(251, 191, 36, 0.15)"
                                        : "rgba(148, 163, 184, 0.15)",
                                      color: ab.epitope_overlap_score > 0.5
                                        ? "#22c55e"
                                        : ab.epitope_overlap_score > 0.3
                                        ? "#fbbf24"
                                        : "var(--muted)",
                                      fontWeight: "600"
                                    }}>
                                      {(ab.epitope_overlap_score * 100).toFixed(1)}%
                                    </span>
                                  </td>
                                  <td>
                                    {ab.resolution ? `${ab.resolution} Å` : "—"}
                                  </td>
                                  <td>
                                    <a
                                      href={`https://www.rcsb.org/structure/${ab.pdb_id}`}
                                      target="_blank"
                                      rel="noopener noreferrer"
                                      style={{
                                        fontSize: "0.85rem",
                                        color: "var(--accent)",
                                        textDecoration: "none",
                                        fontWeight: "500"
                                      }}
                                      onMouseEnter={(e) => e.currentTarget.style.textDecoration = "underline"}
                                      onMouseLeave={(e) => e.currentTarget.style.textDecoration = "none"}
                                    >
                                      View Structure →
                                    </a>
                                  </td>
                                </tr>
                              ))}
                            </tbody>
                          </table>
                        </div>
                      </div>
                    )}
                    
                    {(!match.commercial_antibodies || match.commercial_antibodies.length === 0) &&
                     (!match.structural_antibodies || match.structural_antibodies.length === 0) && (
                      <p style={{ color: "var(--muted)", fontStyle: "italic", marginTop: "12px" }}>
                        No antibodies found for this epitope.
                      </p>
                    )}
                  </div>
                ))}
                {(!results.epitope_matches || results.epitope_matches.length === 0) && (
                  <p style={{ color: "var(--muted)", fontStyle: "italic" }}>
                    No matches found. Try adjusting search parameters or check your Step 2 input.
                  </p>
                )}
              </div>
            )}
            {!loading && !results && !error && (
              <div className="field-help">
                Enter Step 2 JSON path and click "Find Existing Antibodies" to search.
              </div>
            )}
          </div>
        </div>
      </div>
    </section>
  );
};

