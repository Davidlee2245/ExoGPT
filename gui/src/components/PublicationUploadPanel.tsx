import React, { useState, useEffect } from "react";

interface Publication {
  filename: string;
  size: number;
  size_mb: number;
  modified: number;
  extension: string;
  is_csv: boolean;
}

interface PublicationDetails {
  publication: {
    title?: string;
    journal?: string;
    authors?: string;
    year?: string;
    doi?: string;
    pmid?: string;
  };
  biomarkers: Array<{
    gene_symbol?: string;
    uniprot?: string;
    analyte_type?: string;
    disease?: string;
    biofluid?: string;
    logfc?: number;
    p_value?: number;
    fdr?: number;
    method?: string;
  }>;
}

interface PublicationUploadPanelProps {
  dataType?: "protein" | "mrna";  // Determines which data directory to use
}

export const PublicationUploadPanel: React.FC<PublicationUploadPanelProps> = ({ dataType = "protein" }) => {
  const [error, setError] = useState<string | null>(null);
  const [query, setQuery] = useState("");
  const [searching, setSearching] = useState(false);
  const [searchResults, setSearchResults] = useState<any[]>([]);
  const [publications, setPublications] = useState<Publication[]>([]);
  const [loadingPublications, setLoadingPublications] = useState(false);
  const [selectedPublication, setSelectedPublication] = useState<string | null>(null);
  const [publicationDetails, setPublicationDetails] = useState<PublicationDetails | null>(null);
  const [loadingDetails, setLoadingDetails] = useState(false);
  const [extracting, setExtracting] = useState<string | null>(null);
  const [extractionStatus, setExtractionStatus] = useState<{ [filename: string]: { success: boolean; message: string } }>({});

  // Load publications on component mount and after uploads
  const loadPublications = async () => {
    setLoadingPublications(true);
    try {
      const endpoint = dataType === "mrna" 
        ? "http://localhost:5000/api/mrna/step1/publication/list"
        : "http://localhost:5000/api/step1/publication/list";
      const response = await fetch(endpoint);
      const data = await response.json();
      
      if (data.success) {
        setPublications(data.publications || []);
      }
    } catch (err) {
      console.error("Failed to load publications:", err);
    } finally {
      setLoadingPublications(false);
    }
  };

  useEffect(() => {
    loadPublications();
  }, []);

  const loadPublicationDetails = async (filename: string) => {
    if (selectedPublication === filename && publicationDetails) {
      return; // Already loaded
    }
    
    setLoadingDetails(true);
    setSelectedPublication(filename);
    
    try {
      const encodedFilename = encodeURIComponent(filename);
      const endpoint = dataType === "mrna"
        ? `http://localhost:5000/api/mrna/step1/publication/details/${encodedFilename}`
        : `http://localhost:5000/api/step1/publication/details/${encodedFilename}`;
      const response = await fetch(endpoint);
      const data = await response.json();
      
      if (data.success) {
        setPublicationDetails(data);
      } else {
        setPublicationDetails(null);
      }
    } catch (err) {
      console.error("Failed to load publication details:", err);
      setPublicationDetails(null);
    } finally {
      setLoadingDetails(false);
    }
  };

  const handleExtract = async (filename: string) => {
    setExtracting(filename);
    setError(null);
    setExtractionStatus(prev => ({ ...prev, [filename]: { success: false, message: "Extracting..." } }));
    
    try {
      const endpoint = dataType === "mrna"
        ? "http://localhost:5000/api/mrna/step1/publication/extract"
        : "http://localhost:5000/api/step1/publication/extract";
      const response = await fetch(endpoint, {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
        body: JSON.stringify({ filename }),
      });
      
      const data = await response.json();
      
      if (!response.ok || !data.success) {
        throw new Error(data.error || "Extraction failed");
      }
      
      setExtractionStatus(prev => ({
        ...prev,
        [filename]: {
          success: true,
          message: `✓ Extracted ${data.extracted_data?.biomarkers_found || 0} biomarkers successfully!`
        }
      }));
      
      // Reload publication details to show new data
      if (selectedPublication === filename) {
        await loadPublicationDetails(filename);
      }
      
      // Reload publications list
      await loadPublications();
      
    } catch (err) {
      const errorMsg = err instanceof Error ? err.message : "Extraction failed";
      setExtractionStatus(prev => ({
        ...prev,
        [filename]: {
          success: false,
          message: `❌ ${errorMsg}`
        }
      }));
      setError(errorMsg);
    } finally {
      setExtracting(null);
    }
  };


  const handleRAGSearch = async () => {
    if (!query.trim()) {
      setError("Please enter a search query");
      return;
    }

    setSearching(true);
    setError(null);

    try {
      const endpoint = dataType === "mrna"
        ? "http://localhost:5000/api/mrna/step1/publication/search"
        : "http://localhost:5000/api/step1/publication/search";
      const response = await fetch(endpoint, {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
        body: JSON.stringify({
          query: query,
        }),
      });

      const data = await response.json();

      if (!response.ok || !data.success) {
        throw new Error(data.error || "Search failed");
      }

      setSearchResults(data.biomarkers || []);
    } catch (err) {
      setError(err instanceof Error ? err.message : "Search failed");
    } finally {
      setSearching(false);
    }
  };

  return (
    <section className="panel">
      <header className="panel-header">
        <h2>Step 1 · Publication Upload (RAG-based Biomarker Search)</h2>
        <p>
          View uploaded publications and their extracted biomarkers. Use RAG (Retrieval-Augmented Generation) to search across publications.
          Publications are automatically processed when uploaded to <code>data/publications/</code> directory.
        </p>
      </header>
      <div className="panel-body">
        <div className="panel-grid" style={{ gridTemplateColumns: "1fr 1fr" }}>
          <div className="panel-card" style={{ gridColumn: "span 1" }}>
            <div style={{ display: "flex", justifyContent: "space-between", alignItems: "center", marginBottom: "12px" }}>
              <h3 style={{ margin: 0 }}>Uploaded Publications</h3>
              <button
                onClick={loadPublications}
                disabled={loadingPublications}
                style={{
                  padding: "6px 12px",
                  fontSize: "0.85rem",
                  backgroundColor: "transparent",
                  border: "1px solid var(--border-subtle)",
                  borderRadius: "6px",
                  color: "var(--text)",
                  cursor: loadingPublications ? "not-allowed" : "pointer",
                  opacity: loadingPublications ? 0.5 : 1
                }}
              >
                {loadingPublications ? "Loading..." : "🔄 Refresh"}
              </button>
            </div>
            {loadingPublications ? (
              <p>Loading publications...</p>
            ) : publications.length === 0 ? (
              <p className="field-help">No publications found in data/publications/</p>
            ) : (
              <div className="publications-list">
                <p className="field-help" style={{ marginBottom: "8px" }}>
                  Found {publications.length} publication(s) in <code>data/publications/</code>
                </p>
                <div style={{ maxHeight: "400px", overflowY: "auto" }}>
                  {publications.map((pub, idx) => (
                    <div
                      key={idx}
                      onClick={() => loadPublicationDetails(pub.filename)}
                      style={{
                        padding: "12px",
                        marginBottom: "8px",
                        border: selectedPublication === pub.filename ? "2px solid var(--accent)" : "1px solid var(--border-subtle)",
                        borderRadius: "8px",
                        cursor: "pointer",
                        backgroundColor: selectedPublication === pub.filename ? "rgba(56, 189, 248, 0.1)" : "transparent",
                        transition: "all 0.2s"
                      }}
                    >
                      <div style={{ display: "flex", justifyContent: "space-between", alignItems: "start", marginBottom: "8px" }}>
                        <strong style={{ fontSize: "0.95rem" }}>{pub.filename}</strong>
                        <span style={{
                          padding: "2px 8px",
                          borderRadius: "4px",
                          fontSize: "0.75rem",
                          backgroundColor: pub.is_csv ? "rgba(56, 189, 248, 0.2)" : "rgba(148, 163, 184, 0.2)",
                          color: pub.is_csv ? "var(--accent)" : "var(--muted)"
                        }}>
                          {pub.is_csv ? "CSV/TSV" : pub.extension.toUpperCase().replace(".", "")}
                        </span>
                      </div>
                      <div style={{ display: "flex", justifyContent: "space-between", alignItems: "center", marginTop: "8px" }}>
                        <div style={{ display: "flex", gap: "16px", fontSize: "0.85rem", color: "var(--muted)" }}>
                          <span>{pub.size_mb} MB</span>
                          <span>{new Date(pub.modified * 1000).toLocaleDateString()}</span>
                        </div>
                        <button
                          onClick={(e) => {
                            e.stopPropagation(); // Prevent card selection
                            handleExtract(pub.filename);
                          }}
                          disabled={extracting === pub.filename || extracting !== null}
                          style={{
                            padding: "6px 12px",
                            fontSize: "0.8rem",
                            backgroundColor: extracting === pub.filename ? "var(--accent-soft)" : "rgba(56, 189, 248, 0.2)",
                            border: "1px solid var(--accent)",
                            borderRadius: "6px",
                            color: "var(--accent)",
                            cursor: extracting === pub.filename || extracting !== null ? "not-allowed" : "pointer",
                            opacity: extracting === pub.filename || extracting !== null ? 0.6 : 1,
                            fontWeight: "600"
                          }}
                        >
                          {extracting === pub.filename ? "⏳ Extracting..." : "🤖 AI Extractor"}
                        </button>
                      </div>
                      {extractionStatus[pub.filename] && (
                        <div style={{
                          marginTop: "8px",
                          padding: "8px",
                          borderRadius: "4px",
                          fontSize: "0.85rem",
                          backgroundColor: extractionStatus[pub.filename].success 
                            ? "rgba(34, 197, 94, 0.1)" 
                            : "rgba(239, 68, 68, 0.1)",
                          color: extractionStatus[pub.filename].success 
                            ? "rgb(34, 197, 94)" 
                            : "rgb(239, 68, 68)",
                          border: `1px solid ${extractionStatus[pub.filename].success ? "rgb(34, 197, 94)" : "rgb(239, 68, 68)"}`
                        }}>
                          {extractionStatus[pub.filename].message}
                        </div>
                      )}
                      {selectedPublication === pub.filename && (
                        <div style={{ marginTop: "12px", paddingTop: "12px", borderTop: "1px solid var(--border-subtle)" }}>
                          {loadingDetails ? (
                            <p style={{ fontSize: "0.85rem", color: "var(--muted)" }}>Loading details...</p>
                          ) : publicationDetails ? (
                            <>
                              {publicationDetails.publication.title && (
                                <p style={{ margin: "4px 0", fontSize: "0.85rem" }}>
                                  <strong>Title:</strong> {publicationDetails.publication.title}
                                </p>
                              )}
                              {publicationDetails.publication.journal && (
                                <p style={{ margin: "4px 0", fontSize: "0.85rem" }}>
                                  <strong>Journal:</strong> {publicationDetails.publication.journal}
                                </p>
                              )}
                              {publicationDetails.publication.authors && (
                                <p style={{ margin: "4px 0", fontSize: "0.85rem" }}>
                                  <strong>Authors:</strong> {publicationDetails.publication.authors}
                                </p>
                              )}
                              {publicationDetails.publication.year && (
                                <p style={{ margin: "4px 0", fontSize: "0.85rem" }}>
                                  <strong>Year:</strong> {publicationDetails.publication.year}
                                </p>
                              )}
                              {publicationDetails.biomarkers && publicationDetails.biomarkers.length > 0 && (
                                <div style={{ marginTop: "12px" }}>
                                  <strong style={{ fontSize: "0.85rem" }}>
                                    Biomarkers ({publicationDetails.biomarkers.length}):
                                  </strong>
                                  <div style={{ marginTop: "8px", maxHeight: "200px", overflowY: "auto" }}>
                                    <table className="biomarker-table" style={{ fontSize: "0.8rem" }}>
                                      <thead>
                                        <tr>
                                          <th>Gene</th>
                                          <th>UniProt</th>
                                          <th>Type</th>
                                          <th>Disease</th>
                                          <th>Biofluid</th>
                                        </tr>
                                      </thead>
                                      <tbody>
                                        {publicationDetails.biomarkers.map((bm: any, bmIdx: number) => (
                                          <tr key={bmIdx}>
                                            <td><strong>{bm.gene_symbol || "—"}</strong></td>
                                            <td>{bm.uniprot || "—"}</td>
                                            <td>{bm.analyte_type || "protein"}</td>
                                            <td>{bm.disease || "—"}</td>
                                            <td>{bm.biofluid || "—"}</td>
                                          </tr>
                                        ))}
                                      </tbody>
                                    </table>
                                  </div>
                                </div>
                              )}
                            </>
                          ) : (
                            <p style={{ fontSize: "0.85rem", color: "var(--muted)", fontStyle: "italic" }}>
                              No extracted data found in Excel. Upload and process this file to extract biomarkers.
                            </p>
                          )}
                        </div>
                      )}
                    </div>
                  ))}
                </div>
              </div>
            )}
          </div>

          <div className="panel-card" style={{ gridColumn: "span 1" }}>
            <h3>RAG-based Search</h3>
            <label className="field">
              <span className="field-label">Search Query</span>
              <textarea
                value={query}
                onChange={(e) => setQuery(e.target.value)}
                placeholder="e.g., Find biomarkers for melanoma in plasma samples"
                rows={8}
                disabled={searching}
                style={{ width: "100%", minHeight: "200px", fontFamily: "inherit" }}
              />
              <p className="field-help">
                Enter a natural language query to search across uploaded publications
              </p>
            </label>
            <button
              className="run-button"
              onClick={handleRAGSearch}
              disabled={!query.trim() || searching}
            >
              {searching ? "Searching..." : "Search Publications"}
            </button>
          </div>
        </div>

        {error && (
          <div className="error-message">
            <strong>Error:</strong> {error}
          </div>
        )}

        {searchResults.length > 0 && (
          <div className="results-section">
            <h3>Search Results</h3>
            <div className="result-summary">
              <p>
                <strong>Found {searchResults.length} biomarker(s)</strong> matching your query
              </p>
            </div>
            <div className="biomarker-table-container">
              <table className="biomarker-table">
                <thead>
                  <tr>
                    <th>Gene Symbol</th>
                    <th>UniProt</th>
                    <th>Type</th>
                    <th>Log2FC</th>
                    <th>P-value</th>
                    <th>Relevance Score</th>
                    <th>Source</th>
                  </tr>
                </thead>
                <tbody>
                  {searchResults.map((bm: any, idx: number) => (
                    <tr key={idx}>
                      <td><strong>{bm.gene_symbol || "—"}</strong></td>
                      <td>{bm.uniprot || "—"}</td>
                      <td>{bm.analyte_type || "protein"}</td>
                      <td>{bm.logfc ? bm.logfc.toFixed(2) : "—"}</td>
                      <td>{bm.p_value ? bm.p_value.toExponential(2) : "—"}</td>
                      <td>{bm.relevance_score ? bm.relevance_score.toFixed(2) : "—"}</td>
                      <td>{bm.source || "Publication"}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </div>
        )}
      </div>
    </section>
  );
};

