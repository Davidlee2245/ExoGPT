import React, { useState } from "react";

interface Biomarker {
  gene_symbol: string;
  uniprot?: string;
  analyte_type: string;
  surface_likelihood?: number;
  num_studies: number;
  mean_logfc?: number;
  min_p_value?: number;
  min_fdr?: number;
  evidence_level: string;
  score: number;
}

interface Step1Result {
  disease: string;
  biofluid: string;
  normalized_disease_ids?: {
    doid?: string;
    mesh?: string;
  };
  ev_biomarkers: Biomarker[];
  data_sources?: string[];
}

interface DownloadProgress {
  message: string;
  percent: number;
}

export const Step1Panel: React.FC = () => {
  const [disease, setDisease] = useState("Metastatic melanoma");
  const [biofluid, setBiofluid] = useState("Plasma");
  const [loading, setLoading] = useState(false);
  const [downloading, setDownloading] = useState(false);
  const [downloadProgress, setDownloadProgress] = useState<DownloadProgress[]>([]);
  const [downloadPercent, setDownloadPercent] = useState(0);
  const [result, setResult] = useState<Step1Result | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [markdownTable, setMarkdownTable] = useState<string | null>(null);

  const handleRun = async () => {
    setLoading(true);
    setDownloading(false);
    setError(null);
    setResult(null);
    setMarkdownTable(null);
    setDownloadProgress([]);
    setDownloadPercent(0);

    try {
      // First check if data exists
      const checkResponse = await fetch("http://localhost:5000/api/step1/check", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
        body: JSON.stringify({
          disease,
          biofluid,
        }),
      });

      const checkData = await checkResponse.json();
      
      if (!checkData.success) {
        throw new Error(checkData.error || "Failed to check data");
      }

      // If no local data, download it
      if (!checkData.has_local_data) {
        setDownloading(true);
        setDownloadProgress([{ message: "No local data found. Starting download...", percent: 0 }]);
        
        // Download data with progress tracking
        const downloadResponse = await fetch("http://localhost:5000/api/step1/download", {
          method: "POST",
          headers: {
            "Content-Type": "application/json",
          },
          body: JSON.stringify({
            disease,
            biofluid,
          }),
        });

        const downloadData = await downloadResponse.json();
        
        if (!downloadData.success) {
          throw new Error(downloadData.error || "Failed to download data");
        }

        // Update progress
        if (downloadData.progress) {
          setDownloadProgress(downloadData.progress);
          const lastProgress = downloadData.progress[downloadData.progress.length - 1];
          if (lastProgress) {
            setDownloadPercent(lastProgress.percent);
          }
        }

        // Small delay to show completion
        await new Promise(resolve => setTimeout(resolve, 500));
        setDownloading(false);
      }

      // Now run the biomarker finder
      const response = await fetch("http://localhost:5000/api/step1/run", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
        body: JSON.stringify({
          disease,
          biofluid,
          auto_download: false, // Already downloaded if needed
        }),
      });

      const data = await response.json();

      if (!response.ok || !data.success) {
        throw new Error(data.error || "Failed to run Step 1");
      }

      setResult(data.result);
      setMarkdownTable(data.markdown_table || null);
      
      // Save results to localStorage for Step 2 to use
      try {
        localStorage.setItem('step1_result', JSON.stringify(data.result));
        localStorage.setItem('step1_disease', disease);
        localStorage.setItem('step1_biofluid', biofluid);
        // Also save the JSON path if available
        if (data.saved_json_path) {
          // Store relative path for display
          const relativePath = data.saved_json_path.replace(/^.*\/scores\//, './scores/');
          localStorage.setItem('step1_json_path', relativePath);
        }
      } catch (e) {
        console.warn('Failed to save Step 1 results to localStorage:', e);
      }
      
      // Show download info if files were downloaded
      if (data.downloaded_files && data.downloaded_files.length > 0) {
        setDownloadProgress(prev => [
          ...prev,
          { message: `Downloaded ${data.downloaded_files.length} file(s) successfully!`, percent: 100 }
        ]);
      }
    } catch (err) {
      setError(err instanceof Error ? err.message : "Unknown error occurred");
      setDownloading(false);
    } finally {
      setLoading(false);
    }
  };

  const jsonOut = `./scores/${disease.toLowerCase().replace(/\s+/g, "_")}_${biofluid.toLowerCase()}_biomarkers.json`;
  const mdOut = jsonOut.replace(".json", ".md");

  const cli = [
    "python -m exo_gpt.step1_ev_biomarker_finder",
    `  --disease "${disease}"`,
    `  --biofluid "${biofluid}"`,
    `  --out_json "${jsonOut}"`,
    `  --out_md "${mdOut}"`
  ].join(" \\\n");

  return (
    <section className="panel">
      <header className="panel-header">
        <h2>Step 1 · EV Biomarker Finder</h2>
        <p>
          Automatically searches across publications, databases (Vesiclepedia, ExoCarta, EVpedia, EVmiRNA, exoRBase), and user experiments. 
          Normalizes disease and biofluid, aggregates EV datasets, and ranks exosomal biomarkers with a bias toward surface proteins.
        </p>
      </header>
      <div className="panel-body">
        <div className="panel-grid">
          <div className="panel-card">
            <h3>Inputs</h3>
            <label className="field">
              <span className="field-label">Disease</span>
              <input 
                value={disease} 
                onChange={(e) => setDisease(e.target.value)}
                disabled={loading}
                placeholder="e.g., Metastatic melanoma"
              />
            </label>
            <label className="field">
              <span className="field-label">Biofluid</span>
              <input 
                value={biofluid} 
                onChange={(e) => setBiofluid(e.target.value)}
                disabled={loading}
                placeholder="e.g., Plasma"
              />
            </label>
            <p className="field-help">
              Automatically searches: <code>./data/publications/</code>, <code>./data/databases/</code>, <code>./data/experiments/</code>
            </p>
            <button 
              className="run-button" 
              onClick={handleRun}
              disabled={loading || downloading || !disease || !biofluid}
            >
              {downloading ? "Downloading..." : loading ? "Running..." : "Run Step 1"}
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

        {(downloading || downloadProgress.length > 0) && (
          <div className="download-progress-section">
            <h4>Download Progress</h4>
            <div className="progress-bar-container">
              <div className="progress-bar" style={{ width: `${downloadPercent}%` }}></div>
            </div>
            <div className="progress-log">
              {downloadProgress.map((item, idx) => (
                <div key={idx} className="progress-log-item">
                  <span className="progress-percent">{Math.round(item.percent)}%</span>
                  <span className="progress-message">{item.message}</span>
                </div>
              ))}
            </div>
          </div>
        )}

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
              {result.normalized_disease_ids && (
                <p className="field-help">
                  {result.normalized_disease_ids.doid && `DOID: ${result.normalized_disease_ids.doid}`}
                  {result.normalized_disease_ids.mesh && ` | MeSH: ${result.normalized_disease_ids.mesh}`}
                </p>
              )}
              <p>
                <strong>Found {result.ev_biomarkers.length} biomarkers</strong>
                {result.data_sources && result.data_sources.length > 0 && (
                  <> from <strong>{result.data_sources.length} data sources</strong></>
                )}
              </p>
              {result.data_sources && result.data_sources.length > 0 && (
                <p className="field-help">
                  Sources: {result.data_sources.slice(0, 3).join(", ")}
                  {result.data_sources.length > 3 && ` + ${result.data_sources.length - 3} more`}
                </p>
              )}
            </div>

            {result.ev_biomarkers.length > 0 && (
              <div className="biomarker-table-container">
                <table className="biomarker-table">
                  <thead>
                    <tr>
                      <th>Rank</th>
                      <th>Gene Symbol</th>
                      <th>UniProt</th>
                      <th>Type</th>
                      <th>Surface Likelihood</th>
                      <th>Studies</th>
                      <th>Mean logFC</th>
                      <th>Min p-value</th>
                      <th>Evidence</th>
                      <th>Score</th>
                    </tr>
                  </thead>
                  <tbody>
                    {result.ev_biomarkers.slice(0, 50).map((bm, idx) => (
                      <tr key={idx}>
                        <td>{idx + 1}</td>
                        <td><strong>{bm.gene_symbol}</strong></td>
                        <td>{bm.uniprot || "—"}</td>
                        <td>{bm.analyte_type}</td>
                        <td>
                          {bm.surface_likelihood !== null && bm.surface_likelihood !== undefined
                            ? bm.surface_likelihood.toFixed(2)
                            : "—"}
                        </td>
                        <td>{bm.num_studies}</td>
                        <td>
                          {bm.mean_logfc !== null && bm.mean_logfc !== undefined
                            ? bm.mean_logfc.toFixed(2)
                            : "—"}
                        </td>
                        <td>
                          {bm.min_p_value !== null && bm.min_p_value !== undefined
                            ? bm.min_p_value.toExponential(2)
                            : "—"}
                        </td>
                        <td>
                          <span className={`evidence-badge evidence-${bm.evidence_level}`}>
                            {bm.evidence_level}
                          </span>
                        </td>
                        <td><strong>{bm.score.toFixed(2)}</strong></td>
                      </tr>
                    ))}
                  </tbody>
                </table>
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
          </div>
        )}
      </div>
    </section>
  );
};

