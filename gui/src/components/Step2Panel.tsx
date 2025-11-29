import React, { useMemo, useState } from "react";

interface EpitopePatch {
  epitope_id: string;
  epitope_residues: number[];
  epitope_description: string;
  selection_method: string;
  notes?: string;
}

interface CuratedTarget {
  target_name: string;
  gene_symbol: string;
  uniprot?: string;
  structure_type: string;
  model_path: string;
  pdb_id?: string;
  domain: string;
  domain_residues?: [number, number];
  epitopes: EpitopePatch[];
  selection_rationale: string;
  surface_likelihood?: number;
  biomarker_score?: number;
}

interface Step2Result {
  disease: string;
  biofluid: string;
  user_specified_target?: string;
  curated_targets: CuratedTarget[];
  structure_mapping_file: string;
  num_targets: number;
  methodology_notes?: string;
}

export const Step2Panel: React.FC = () => {
  const [step1JsonPath, setStep1JsonPath] = useState("./scores/metastatic_melanoma_plasma_biomarkers.json");
  const [target, setTarget] = useState("");
  const [maxTargets, setMaxTargets] = useState(5);
  const [loading, setLoading] = useState(false);
  const [result, setResult] = useState<Step2Result | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [markdownTable, setMarkdownTable] = useState<string | null>(null);
  const [step1JsonData, setStep1JsonData] = useState<any>(null); // For passing Step 1 results directly

  const handleRun = async () => {
    setLoading(true);
    setError(null);
    setResult(null);
    setMarkdownTable(null);

    try {
      const requestBody: any = {
        max_targets: maxTargets,
      };

      // If step1JsonData is provided, use it; otherwise use path
      if (step1JsonData) {
        requestBody.step1_json_data = step1JsonData;
      } else {
        requestBody.step1_json_path = step1JsonPath;
      }

      if (target) {
        requestBody.target = target;
      }

      const response = await fetch("http://localhost:5000/api/step2/run", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
        body: JSON.stringify(requestBody),
      });

      const data = await response.json();

      if (!response.ok || !data.success) {
        throw new Error(data.error || "Failed to run Step 2");
      }

      setResult(data.result);
      setMarkdownTable(data.markdown_table || null);
      
      // Save results to localStorage for Step 3 to use
      try {
        localStorage.setItem('step2_result', JSON.stringify(data.result));
        // Also save a default path for reference
        const defaultPath = `./epitopes/${data.result.disease?.toLowerCase().replace(/\s+/g, "_")}_${data.result.biofluid?.toLowerCase()}_epitopes.json`;
        localStorage.setItem('step2_json_path', defaultPath);
      } catch (e) {
        console.warn('Failed to save Step 2 results to localStorage:', e);
      }
    } catch (err) {
      setError(err instanceof Error ? err.message : "Unknown error occurred");
    } finally {
      setLoading(false);
    }
  };

  // Load Step 1 results from localStorage if available
  React.useEffect(() => {
    try {
      const savedResult = localStorage.getItem('step1_result');
      const savedPath = localStorage.getItem('step1_json_path');
      const savedDisease = localStorage.getItem('step1_disease');
      const savedBiofluid = localStorage.getItem('step1_biofluid');
      
      if (savedResult) {
        const result = JSON.parse(savedResult);
        setStep1JsonData(result);
        // Use saved path if available, otherwise generate from disease/biofluid
        if (savedPath) {
          setStep1JsonPath(savedPath);
        } else if (savedDisease && savedBiofluid) {
          const path = `./scores/${savedDisease.toLowerCase().replace(/\s+/g, "_")}_${savedBiofluid.toLowerCase()}_biomarkers.json`;
          setStep1JsonPath(path);
        }
      }
    } catch (e) {
      console.warn('Failed to load Step 1 results from localStorage:', e);
    }
  }, []);

  const cli = [
    "python -m exo_gpt.step2_target_epitope_curator",
    `  --step1_json "${step1JsonPath}"`,
    target ? `  --target "${target}"` : "",
    `  --max_targets ${maxTargets}`,
    `  --out_json "./epitopes/mm_plasma_epitopes.json"`,
    `  --out_md "./epitopes/mm_plasma_epitopes.md"`
  ].filter(Boolean).join(" \\\n");

  const structureBreakdown = useMemo(() => {
    if (!result) return null;
    return result.curated_targets.reduce<Record<string, number>>((acc, target) => {
      const rawType = target.structure_type || "unknown";
      const label = rawType.toLowerCase() === "pdb" ? "PDB" : rawType.toLowerCase() === "alphafold" ? "AlphaFold" : rawType;
      acc[label] = (acc[label] || 0) + 1;
      return acc;
    }, {});
  }, [result]);

  return (
    <section className="panel">
      <header className="panel-header">
        <h2>Step 2 · Target & Epitope Curator</h2>
        <p>
          Select nanobinder-suitable exosomal surface proteins from Step 1 results. Map to structures, extract extracellular domains, 
          and propose 1–3 solvent-exposed loop epitopes per target.
        </p>
      </header>
      <div className="panel-body">
        <div className="panel-grid">
          <div className="panel-card">
            <h3>Inputs</h3>
            <label className="field">
              <span className="field-label">Step 1 JSON Result Path</span>
              <input 
                value={step1JsonPath} 
                onChange={(e) => {
                  setStep1JsonPath(e.target.value);
                  setStep1JsonData(null); // Clear direct data when path is used
                }}
                placeholder="./scores/..._biomarkers.json"
                disabled={loading || !!step1JsonData}
              />
            </label>
            {step1JsonData && (
              <p className="field-help" style={{ color: "green" }}>
                ✓ Using Step 1 results directly (from Step 1 panel)
              </p>
            )}
            <label className="field">
              <span className="field-label">Target (Optional)</span>
              <input 
                value={target} 
                onChange={(e) => setTarget(e.target.value)}
                placeholder="e.g., CD274 or PD-L1"
                disabled={loading}
              />
            </label>
            <label className="field">
              <span className="field-label">Max Targets</span>
              <input 
                type="number"
                value={maxTargets} 
                onChange={(e) => setMaxTargets(parseInt(e.target.value) || 5)}
                min={1}
                max={20}
                disabled={loading}
              />
            </label>
            <p className="field-help">
              Uses structure mapping from <code>./data/structure_mapping.csv</code>
            </p>
            <button 
              className="run-button" 
              onClick={handleRun}
              disabled={loading || (!step1JsonPath && !step1JsonData)}
            >
              {loading ? "Running..." : "Run Step 2"}
            </button>
          </div>
          <div className="panel-card">
            <h3>Command</h3>
            <p className="field-help">Copy-paste into your Exosome-GPT Python environment.</p>
            <pre className="code-block">
              <code>{cli}</code>
            </pre>
          </div>
        </div>

        {error && (
          <div className="error-message">
            <strong>Error:</strong> {error}
          </div>
        )}

        {result && (
          <div className="results-section">
            <h3>Results</h3>
            <div className="result-summary">
              <p>
                <strong>Disease:</strong> {result.disease} | <strong>Biofluid:</strong> {result.biofluid}
              </p>
              {result.user_specified_target && (
                <p>
                  <strong>Specified Target:</strong> {result.user_specified_target}
                </p>
              )}
              <p>
                <strong>Curated {result.num_targets} target(s)</strong>
              </p>
              {structureBreakdown && (
                <p className="structure-summary">
                  <strong>3D Models:</strong>{" "}
                  {Object.entries(structureBreakdown).map(([label, count], idx, arr) => (
                    <React.Fragment key={label}>
                      <span className={`structure-badge badge-${label.toLowerCase()}`}>
                        {label} ({count})
                      </span>
                      {idx < arr.length - 1 && <span className="structure-divider">·</span>}
                    </React.Fragment>
                        ))}
                </p>
              )}
            </div>

            {result.curated_targets.length > 0 && (
              <div className="targets-container">
                <h4>Curated Targets & Epitopes</h4>
                <div className="biomarker-table-container" style={{ marginTop: "15px" }}>
                  <table className="biomarker-table">
                    <thead>
                      <tr>
                        <th>Target</th>
                        <th>Gene</th>
                        <th>UniProt</th>
                        <th>Structure</th>
                        <th>Domain</th>
                        <th>Domain Residues</th>
                        <th>Surface Likelihood</th>
                        <th>Biomarker Score</th>
                        <th>Epitopes</th>
                        <th>Epitope Residues</th>
                      </tr>
                    </thead>
                    <tbody>
                      {result.curated_targets.map((target, idx) => {
                        const structureDisplay = target.pdb_id 
                          ? `${target.structure_type.toUpperCase()} (${target.pdb_id})`
                          : target.structure_type.toUpperCase();
                        const domainResiduesDisplay = target.domain_residues
                          ? `${target.domain_residues[0]}-${target.domain_residues[1]}`
                          : "—";
                        
                        // Create rows for each epitope, or one row if no epitopes
                        if (target.epitopes.length === 0) {
                          return (
                            <tr key={idx}>
                              <td><strong>{target.target_name}</strong></td>
                              <td>{target.gene_symbol}</td>
                              <td>{target.uniprot || "—"}</td>
                              <td>{structureDisplay}</td>
                              <td>{target.domain}</td>
                              <td>{domainResiduesDisplay}</td>
                              <td>
                                {target.surface_likelihood !== null && target.surface_likelihood !== undefined
                                  ? target.surface_likelihood.toFixed(2)
                                  : "—"}
                              </td>
                              <td>
                                {target.biomarker_score !== null && target.biomarker_score !== undefined
                                  ? target.biomarker_score.toFixed(2)
                                  : "—"}
                              </td>
                              <td>None</td>
                              <td>—</td>
                            </tr>
                          );
                        }
                        
                        // Multiple rows if multiple epitopes - show all target info in each row
                        return target.epitopes.map((ep, epIdx) => (
                          <tr key={`${idx}-${epIdx}`}>
                            <td><strong>{target.target_name}</strong></td>
                            <td>{target.gene_symbol}</td>
                            <td>{target.uniprot || "—"}</td>
                            <td>{structureDisplay}</td>
                            <td>{target.domain}</td>
                            <td>{domainResiduesDisplay}</td>
                            <td>
                              {target.surface_likelihood !== null && target.surface_likelihood !== undefined
                                ? target.surface_likelihood.toFixed(2)
                                : "—"}
                            </td>
                            <td>
                              {target.biomarker_score !== null && target.biomarker_score !== undefined
                                ? target.biomarker_score.toFixed(2)
                                : "—"}
                            </td>
                            <td>
                              <strong>{ep.epitope_id}</strong>
                              <br />
                              <span style={{ fontSize: "0.85em", color: "#666" }}>
                                {ep.epitope_description}
                              </span>
                            </td>
                            <td>{ep.epitope_residues.join(", ")}</td>
                          </tr>
                        ));
                      })}
                    </tbody>
                  </table>
                </div>
              </div>
            )}

            {markdownTable && (
              <div className="markdown-output">
                <h4>Markdown Summary</h4>
                <pre className="markdown-block">
                  <code>{markdownTable}</code>
                </pre>
              </div>
            )}

            {result.debug_info && (
              <div
                className="debug-info"
                style={{
                  marginTop: "18px",
                  padding: "14px 16px",
                  backgroundColor: "#f8fafc",
                  borderRadius: "10px",
                  fontSize: "0.9em",
                  color: "#0f172a",
                  border: "1px solid #cbd5f5",
                }}
              >
                <strong style={{ display: "block", marginBottom: "6px" }}>Debug Info:</strong>
                <ul style={{ margin: "5px 0", paddingLeft: "22px", color: "#0f172a" }}>
                  <li>Total biomarkers: {result.debug_info.total_biomarkers}</li>
                  <li>Suitable biomarkers: {result.debug_info.suitable_biomarkers}</li>
                  <li>Selected biomarkers: {result.debug_info.selected_biomarkers}</li>
                  <li>Structure mapping file: {result.debug_info.structure_mapping_file}</li>
                  <li>Structure mapping exists: {result.debug_info.structure_mapping_exists ? "Yes" : "No"}</li>
                  <li>Structure mapping loaded: {result.debug_info.structure_mapping_loaded ? "Yes" : "No"}</li>
                  <li>Structure mapping rows: {result.debug_info.structure_mapping_rows}</li>
                  <li>Skipped (no structure): {result.debug_info.skipped_no_structure}</li>
                </ul>
                {result.debug_info.skipped_biomarkers && result.debug_info.skipped_biomarkers.length > 0 && (
                  <div style={{ marginTop: "12px" }}>
                    <strong>Skipped Biomarkers (Detailed):</strong>
                    {result.debug_info.skipped_biomarkers.map((skipped: any, idx: number) => (
                      <div
                        key={idx}
                        style={{
                          marginTop: "8px",
                          padding: "6px 8px",
                          backgroundColor: "#ffffff",
                          border: "1px solid #d1d5db",
                          borderRadius: "6px",
                          color: "#111827",
                        }}
                      >
                        <strong>{skipped.gene_symbol}</strong> (UniProt: {skipped.uniprot || "N/A"})<br />
                        Normalized: gene={skipped.gene_normalized}, uniprot={skipped.uniprot_normalized}<br />
                        Gene in mapping: {skipped.gene_in_mapping ? "YES" : "NO"}<br />
                        UniProt in mapping: {skipped.uniprot_in_mapping ? "YES" : "NO"}<br />
                        Mapping genes: {skipped.mapping_genes_sample?.join(", ") || "N/A"}<br />
                        Mapping UniProts: {skipped.mapping_uniprots_sample?.join(", ") || "N/A"}<br />
                        Actual genes in DF: {skipped.actual_genes_in_df?.join(", ") || "N/A"}<br />
                        Actual UniProts in DF: {skipped.actual_uniprots_in_df?.join(", ") || "N/A"}<br />
                        DF columns: {skipped.mapping_df_columns?.join(", ") || "N/A"}
                      </div>
                    ))}
                  </div>
                )}
                {result.debug_info.selected_biomarker_details && result.debug_info.selected_biomarker_details.length > 0 && (
                  <div style={{ marginTop: "12px" }}>
                    <strong>Selected Biomarkers:</strong>
                    <ul style={{ margin: "6px 0", paddingLeft: "22px" }}>
                      {result.debug_info.selected_biomarker_details.map((bm: any, idx: number) => (
                        <li key={idx}>
                          {bm.gene_symbol} (UniProt: {bm.uniprot || "N/A"}, Type: {bm.analyte_type}, Score: {bm.score?.toFixed(2)})
                        </li>
                      ))}
                    </ul>
                  </div>
                )}
              </div>
            )}
          </div>
        )}
      </div>
    </section>
  );
};

