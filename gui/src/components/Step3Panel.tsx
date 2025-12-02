import React, { useState } from "react";

export const Step3Panel: React.FC = () => {
  const [step2JsonPath, setStep2JsonPath] = useState("./epitopes/mm_plasma_epitopes.json");
  const [outputBaseDir, setOutputBaseDir] = useState("./workflows");
  const [binderLength, setBinderLength] = useState(90);
  const [numBackbones, setNumBackbones] = useState(10);
  const [numSequencesPerBackbone, setNumSequencesPerBackbone] = useState(8);
  const [numAfModels, setNumAfModels] = useState(5);
  const [numAfRecycles, setNumAfRecycles] = useState(3);
  const [useColabFold, setUseColabFold] = useState(true);
  const [execute, setExecute] = useState(false);
  const [loading, setLoading] = useState(false);
  const [plan, setPlan] = useState<any>(null);
  const [error, setError] = useState<string | null>(null);
  const [step2JsonData, setStep2JsonData] = useState<any>(null);
  const [progressMessages, setProgressMessages] = useState<Array<{message: string, percent: number}>>([]);
  const [selectedFile, setSelectedFile] = useState<string | null>(null);
  const [fileContents, setFileContents] = useState<string | null>(null);
  const [fileName, setFileName] = useState<string>("");
  const [loadingFile, setLoadingFile] = useState(false);
  const progressRef = React.useRef<HTMLDivElement>(null);

  // Auto-scroll progress to bottom when new messages arrive
  React.useEffect(() => {
    if (progressRef.current) {
      progressRef.current.scrollTop = progressRef.current.scrollHeight;
    }
  }, [progressMessages]);

  const handleRun = async () => {
    setLoading(true);
    setError(null);
    setPlan(null);
    setProgressMessages([]);

    try {
      const requestBody: any = {
        output_base_dir: outputBaseDir,
        binder_length: binderLength,
        num_backbones: numBackbones,
        num_sequences_per_backbone: numSequencesPerBackbone,
        num_af_models: numAfModels,
        num_af_recycles: numAfRecycles,
        use_colabfold: useColabFold,
        execute: execute,
      };

      // If step2JsonData is provided, use it; otherwise use path
      if (step2JsonData) {
        requestBody.step2_json_data = step2JsonData;
      } else {
        requestBody.step2_json_path = step2JsonPath;
      }

      // Use SSE streaming for real-time progress when execute is true
      if (execute) {
        // Use fetch with streaming for Server-Sent Events
        let response: Response;
        try {
          response = await fetch("http://localhost:5000/api/step3/stream", {
            method: "POST",
            headers: {
              "Content-Type": "application/json",
            },
            body: JSON.stringify(requestBody),
          });
        } catch (fetchError) {
          console.error("Fetch error:", fetchError);
          // Fallback to regular endpoint if streaming fails
          console.warn("Streaming endpoint failed, falling back to regular endpoint");
          // Continue to the else block below by setting execute to false temporarily
          // Actually, let's just throw a more helpful error and suggest using non-execute mode
          throw new Error(
            `Failed to connect to streaming endpoint: ${fetchError instanceof Error ? fetchError.message : "Unknown error"}. ` +
            `Make sure the backend is running on http://localhost:5000. ` +
            `You can try unchecking "Execute Tools" to use the regular endpoint.`
          );
        }

        if (!response.ok) {
          let errorData: any;
          try {
            // Try to parse as JSON first (for immediate errors)
            const contentType = response.headers.get("content-type");
            if (contentType && contentType.includes("application/json")) {
              errorData = await response.json();
            } else {
              const text = await response.text();
              errorData = { error: text || `Server returned status ${response.status}: ${response.statusText}` };
            }
          } catch (parseError) {
            errorData = { error: `Server returned status ${response.status}: ${response.statusText}` };
          }
          throw new Error(errorData.error || "Failed to start Step 3");
        }

        // Check if response is actually a stream
        const contentType = response.headers.get("content-type");
        if (!contentType || !contentType.includes("text/event-stream")) {
          // Not a stream, try to parse as JSON
          try {
            const data = await response.json();
            if (data.error) {
              throw new Error(data.error);
            }
            setPlan(data.plan);
            return;
          } catch (jsonError) {
            throw new Error("Unexpected response format from server");
          }
        }

        // Read the stream
        const reader = response.body?.getReader();
        const decoder = new TextDecoder();

        if (!reader) {
          throw new Error("Failed to get response stream");
        }

        let buffer = '';
        let resultReceived = false;
        
        try {
          while (true) {
            const { done, value } = await reader.read();
            if (done) break;

            buffer += decoder.decode(value, { stream: true });
            const lines = buffer.split('\n');
            buffer = lines.pop() || ''; // Keep incomplete line in buffer

            for (const line of lines) {
              if (line.trim() === '' || line.startsWith(':')) {
                // Skip empty lines and comments (keepalive)
                continue;
              }
              
              if (line.startsWith('data: ')) {
                try {
                  const data = JSON.parse(line.slice(6)); // Remove 'data: ' prefix
                  
                  if (data.type === 'progress') {
                    // Add progress message in real-time
                    setProgressMessages(prev => [...prev, {
                      message: data.message,
                      percent: data.percent
                    }]);
                  } else if (data.type === 'result') {
                    // Final result received
                    setPlan(data.plan);
                    resultReceived = true;
                  } else if (data.type === 'error') {
                    // Error received
                    const errorMsg = data.error?.message || data.error || "Unknown error occurred";
                    setError(errorMsg);
                    throw new Error(errorMsg);
                  }
                } catch (err) {
                  if (err instanceof Error && err.message !== "Unknown error occurred") {
                    throw err;
                  }
                  console.error('Error parsing SSE data:', err);
                }
              }
            }
          }
          
          if (!resultReceived) {
            throw new Error("Stream ended without result");
          }
        } finally {
          reader.releaseLock();
        }
      } else {
        // Non-execute mode: use regular API endpoint
        const response = await fetch("http://localhost:5000/api/step3/run", {
          method: "POST",
          headers: {
            "Content-Type": "application/json",
          },
          body: JSON.stringify(requestBody),
        });

        const data = await response.json();

        if (!response.ok || !data.success) {
          throw new Error(data.error || "Failed to run Step 3");
        }

        setPlan(data.plan);
        // Update progress messages - they come all at once but we display them
        if (data.progress && data.progress.length > 0) {
          setProgressMessages(data.progress);
        }
      }
    } catch (err) {
      console.error("Step 3 error:", err);
      const errorMessage = err instanceof Error ? err.message : "Unknown error occurred";
      setError(errorMessage);
      // Also log to console for debugging
      if (err instanceof TypeError && err.message.includes("fetch")) {
        console.error("Network error - check if backend is running on http://localhost:5000");
      }
    } finally {
      setLoading(false);
    }
  };

  // Load Step 2 results from localStorage if available
  React.useEffect(() => {
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

  const handleFileClick = async (filePath: string) => {
    setLoadingFile(true);
    setSelectedFile(filePath);
    setFileContents(null);
    
    try {
      const response = await fetch("http://localhost:5000/api/file/read", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
        body: JSON.stringify({ file_path: filePath }),
      });

      const data = await response.json();

      if (!response.ok || !data.success) {
        throw new Error(data.error || "Failed to read file");
      }

      setFileContents(data.contents);
      setFileName(data.file_name || filePath);
    } catch (err) {
      setError(err instanceof Error ? err.message : "Failed to load file");
      setFileContents(null);
    } finally {
      setLoadingFile(false);
    }
  };

  const cli = [
    "python -m exo_gpt.step3_nanobinder_orchestrator",
    `  --step2_json "${step2JsonPath}"`,
    `  --output_base_dir "${outputBaseDir}"`,
    `  --binder_length ${binderLength}`,
    `  --num_backbones ${numBackbones}`,
    `  --num_sequences_per_backbone ${numSequencesPerBackbone}`,
    `  --num_af_models ${numAfModels}`,
    `  --num_af_recycles ${numAfRecycles}`,
    useColabFold ? "  --use_colabfold" : "  --no-use_colabfold",
    `  --out_plan "./workflows/design_plan.json"`
  ].filter(Boolean).join(" \\\n");

  return (
    <section className="panel">
      <header className="panel-header">
        <h2>Step 3-2 · Nanobinder Design Orchestrator</h2>
        <p>
          Orchestrates RFdiffusion backbone generation, ProteinMPNN sequence design, and AlphaFold-Multimer (ColabFold) 
          structure prediction. Generates reusable scripts and configs for execution on external machines/HPC.
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
              <span className="field-label">Output Base Directory</span>
              <input 
                value={outputBaseDir} 
                onChange={(e) => setOutputBaseDir(e.target.value)}
                placeholder="./workflows"
                disabled={loading}
              />
            </label>
            <label className="field">
              <span className="field-label">Binder Length (aa)</span>
              <input 
                type="number"
                value={binderLength} 
                onChange={(e) => setBinderLength(parseInt(e.target.value) || 90)}
                min={50}
                max={200}
                disabled={loading}
              />
            </label>
            <label className="field">
              <span className="field-label">Number of Backbones</span>
              <input 
                type="number"
                value={numBackbones} 
                onChange={(e) => setNumBackbones(parseInt(e.target.value) || 10)}
                min={1}
                max={100}
                disabled={loading}
              />
            </label>
            <label className="field">
              <span className="field-label">Sequences per Backbone</span>
              <input 
                type="number"
                value={numSequencesPerBackbone} 
                onChange={(e) => setNumSequencesPerBackbone(parseInt(e.target.value) || 8)}
                min={1}
                max={50}
                disabled={loading}
              />
            </label>
            <label className="field">
              <span className="field-label">AlphaFold Models per Sequence</span>
              <input 
                type="number"
                value={numAfModels} 
                onChange={(e) => setNumAfModels(parseInt(e.target.value) || 5)}
                min={1}
                max={20}
                disabled={loading}
              />
            </label>
            <label className="field">
              <span className="field-label">AlphaFold Recycles</span>
              <input 
                type="number"
                value={numAfRecycles} 
                onChange={(e) => setNumAfRecycles(parseInt(e.target.value) || 3)}
                min={1}
                max={10}
                disabled={loading}
              />
            </label>
            <label className="field">
              <span className="field-label">Use ColabFold</span>
              <input 
                type="checkbox"
                checked={useColabFold} 
                onChange={(e) => setUseColabFold(e.target.checked)}
                disabled={loading}
              />
            </label>
            <label className="field">
              <span className="field-label">Execute Tools (if available)</span>
              <input 
                type="checkbox"
                checked={execute} 
                onChange={(e) => setExecute(e.target.checked)}
                disabled={loading}
              />
              <p className="field-help" style={{ fontSize: "0.85em", marginTop: "5px" }}>
                If checked, will attempt to execute RFdiffusion, ProteinMPNN, and ColabFold if they are installed.
                <strong> Scripts are always generated regardless</strong> - you can run them manually on HPC/cluster.
              </p>
            </label>
            <button 
              className="run-button" 
              onClick={handleRun}
              disabled={loading || (!step2JsonPath && !step2JsonData)}
            >
              {loading ? "Generating..." : "Generate Design Plan"}
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

        {progressMessages.length > 0 && (
          <div className="progress-section" style={{ 
            marginTop: "20px", 
            padding: "15px", 
            backgroundColor: "#1e1e1e", 
            borderRadius: "4px",
            border: "1px solid #444",
            color: "#e0e0e0"
          }}>
            <h4 style={{ color: "#fff", marginBottom: "10px" }}>Execution Progress</h4>
            <div 
              ref={progressRef}
              style={{ 
                maxHeight: "400px", 
                overflowY: "auto",
                fontFamily: "monospace",
                fontSize: "0.85em",
                lineHeight: "1.6"
              }}
            >
              {progressMessages.map((msg, idx) => (
                <div 
                  key={idx} 
                  style={{ 
                    marginBottom: "3px", 
                    padding: "2px 0",
                    color: msg.message.toLowerCase().includes("failed") || msg.message.toLowerCase().includes("error") 
                      ? "#ff6b6b" 
                      : msg.message.toLowerCase().includes("success") || msg.message.toLowerCase().includes("completed")
                      ? "#51cf66"
                      : "#e0e0e0"
                  }}
                >
                  <span style={{ color: "#4dabf7", fontWeight: "bold" }}>
                    [{msg.percent.toFixed(1)}%]
                  </span>{" "}
                  <span>{msg.message}</span>
                </div>
              ))}
            </div>
          </div>
        )}

        {plan && plan.success && (
          <div className="results-section">
            <h3>Design Plan Generated</h3>
            <div className="result-summary">
              <p>
                <strong>Disease:</strong> {plan.disease} | <strong>Biofluid:</strong> {plan.biofluid}
              </p>
              <p>
                <strong>{plan.num_designs} design(s) generated</strong>
              </p>
              <p>
                <strong>Output Base Directory:</strong> {plan.output_base_dir}
              </p>
              <p>
                <strong>Generated {plan.generated_files?.length || 0} files</strong>
              </p>
            </div>

            {plan.designs && plan.designs.length > 0 && (
              <div className="designs-container">
                <h4>Design Summary</h4>
                <div className="biomarker-table-container" style={{ marginTop: "15px" }}>
                  <table className="biomarker-table">
                    <thead>
                      <tr>
                        <th>Target</th>
                        <th>Gene</th>
                        <th>Epitope ID</th>
                        <th>Epitope Residues</th>
                        <th>Binder Length (aa)</th>
                        <th>Backbones</th>
                        <th>Seqs/Backbone</th>
                        <th>Total Sequences</th>
                        <th>AF Models/Seq</th>
                        <th>Total AF Predictions</th>
                        <th>Output Directory</th>
                      </tr>
                    </thead>
                    <tbody>
                      {plan.designs.map((design: any, idx: number) => (
                        <tr key={idx}>
                          <td><strong>{design.target_name}</strong></td>
                          <td>{design.gene_symbol}</td>
                          <td>{design.epitope_id}</td>
                          <td>{design.epitope_residues.join(", ")}</td>
                          <td>{design.binder_length}</td>
                          <td>{design.num_backbones}</td>
                          <td>{design.num_sequences_per_backbone}</td>
                          <td><strong>{design.total_sequences}</strong></td>
                          <td>{design.num_af_models}</td>
                          <td><strong>{design.total_af_predictions}</strong></td>
                          <td>
                            <code style={{ fontSize: "0.85em" }}>{design.output_directory}</code>
                          </td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
              </div>
            )}

            {plan.execution_mode && plan.tool_status && (
              <div className="tool-status" style={{ 
                marginTop: "20px", 
                padding: "15px", 
                backgroundColor: "#f9f9f9", 
                borderRadius: "4px",
                border: "1px solid #ddd"
              }}>
                <h4>Tool Availability</h4>
                <p className="field-help" style={{ marginBottom: "15px", fontStyle: "italic" }}>
                  Note: Scripts are still generated even if tools are not installed. You can run them manually on a machine with the tools installed.
                </p>
                {Object.entries(plan.tool_status).map(([tool, status]: [string, any]) => {
                  const toolInfo: Record<string, {name: string, install: string, docs: string}> = {
                    rfdiffusion: {
                      name: "RFdiffusion",
                      install: "git clone https://github.com/RosettaCommons/RFdiffusion.git",
                      docs: "https://github.com/RosettaCommons/RFdiffusion"
                    },
                    proteinmpnn: {
                      name: "ProteinMPNN",
                      install: "git clone https://github.com/dauparas/ProteinMPNN.git",
                      docs: "https://github.com/dauparas/ProteinMPNN"
                    },
                    colabfold: {
                      name: "ColabFold",
                      install: "pip install colabfold",
                      docs: "https://github.com/sokrypton/ColabFold"
                    },
                    alphafold: {
                      name: "AlphaFold-Multimer",
                      install: "Follow installation guide at DeepMind",
                      docs: "https://github.com/deepmind/alphafold"
                    }
                  };
                  
                  const info = toolInfo[tool] || { name: tool, install: "See documentation", docs: "" };
                  
                  return (
                    <div key={tool} style={{ 
                      marginBottom: "15px", 
                      padding: "10px", 
                      backgroundColor: status.available ? "#e8f5e9" : "#fff3e0",
                      borderRadius: "4px",
                      border: `1px solid ${status.available ? "#4caf50" : "#ff9800"}`
                    }}>
                      {status.available ? (
                        <div>
                          <p style={{ color: "green", margin: "5px 0", fontWeight: "bold" }}>
                            ✓ <strong>{info.name}</strong>: Available
                          </p>
                          <p style={{ fontSize: "0.9em", color: "#666", margin: "5px 0" }}>
                            Command: <code>{status.command}</code>
                          </p>
                        </div>
                      ) : (
                        <div>
                          <p style={{ color: "#ff6f00", margin: "5px 0", fontWeight: "bold" }}>
                            ⚠ <strong>{info.name}</strong>: Not installed
                          </p>
                          <p style={{ fontSize: "0.9em", color: "#666", margin: "5px 0" }}>
                            {status.error || "Tool not found"}
                          </p>
                          <details style={{ marginTop: "8px" }}>
                            <summary style={{ cursor: "pointer", color: "#1976d2", fontSize: "0.9em" }}>
                              Installation Instructions
                            </summary>
                            <div style={{ marginTop: "8px", padding: "8px", backgroundColor: "#fff", borderRadius: "4px" }}>
                              <p style={{ margin: "5px 0", fontSize: "0.85em" }}>
                                <strong>Install:</strong> <code>{info.install}</code>
                              </p>
                              <p style={{ margin: "5px 0", fontSize: "0.85em" }}>
                                <strong>Documentation:</strong>{" "}
                                <a href={info.docs} target="_blank" rel="noopener noreferrer" style={{ color: "#1976d2" }}>
                                  {info.docs}
                                </a>
                              </p>
                            </div>
                          </details>
                        </div>
                      )}
                    </div>
                  );
                })}
              </div>
            )}

            {plan.execution_mode && plan.execution_results && (
              <div className="execution-results" style={{ 
                marginTop: "20px", 
                padding: "15px", 
                backgroundColor: "#fff9e6", 
                borderRadius: "4px",
                border: "1px solid #ffd700"
              }}>
                <h4>Execution Results</h4>
                {Object.entries(plan.execution_results).map(([designKey, results]: [string, any]) => (
                  <div key={designKey} style={{ marginBottom: "15px", padding: "10px", backgroundColor: "#fff", borderRadius: "4px" }}>
                    <h5>{designKey}</h5>
                    {Object.entries(results).map(([step, result]: [string, any]) => (
                      <div key={step} style={{ marginLeft: "15px", marginBottom: "5px" }}>
                        {result.success ? (
                          <p style={{ color: "green", margin: "2px 0" }}>
                            ✓ <strong>{step}</strong>: {result.generated_files?.length || 0} files generated
                          </p>
                        ) : (
                          <p style={{ color: "red", margin: "2px 0" }}>
                            ✗ <strong>{step}</strong>: {result.error || "Failed"}
                          </p>
                        )}
                      </div>
                    ))}
                  </div>
                ))}
              </div>
            )}

            {plan.pipeline_steps && (
              <div className="pipeline-steps">
                <h4>Pipeline Steps</h4>
                <ol>
                  {plan.pipeline_steps.map((step: string, idx: number) => (
                    <li key={idx}>{step}</li>
                  ))}
                </ol>
              </div>
            )}

            {plan.generated_files && plan.generated_files.length > 0 && (
              <div className="generated-files">
                <h4>Generated Files ({plan.generated_files.length})</h4>
                <p className="field-help" style={{ marginBottom: "10px" }}>
                  Click on any file to view its contents
                </p>
                <div style={{ 
                  maxHeight: "300px", 
                  overflowY: "auto", 
                  border: "1px solid #ddd", 
                  borderRadius: "4px",
                  padding: "10px",
                  backgroundColor: "#fff"
                }}>
                  <ul style={{ listStyle: "none", padding: 0, margin: 0 }}>
                    {plan.generated_files.map((file: string, idx: number) => {
                      const fileName = file.split('/').pop() || file;
                      const fileDir = file.substring(0, file.lastIndexOf('/'));
                      
                      return (
                        <li 
                          key={idx} 
                          style={{ 
                            marginBottom: "5px",
                            cursor: "pointer",
                            padding: "8px",
                            borderRadius: "4px",
                            backgroundColor: selectedFile === file ? "#e3f2fd" : "#f5f5f5",
                            border: selectedFile === file ? "2px solid #2196f3" : "1px solid #ddd",
                            transition: "background-color 0.2s"
                          }}
                          onClick={() => handleFileClick(file)}
                          onMouseEnter={(e) => {
                            if (selectedFile !== file) {
                              e.currentTarget.style.backgroundColor = "#f0f0f0";
                            }
                          }}
                          onMouseLeave={(e) => {
                            if (selectedFile !== file) {
                              e.currentTarget.style.backgroundColor = "#f5f5f5";
                            }
                          }}
                          title={file}
                        >
                          <div style={{ display: "flex", alignItems: "center", gap: "8px" }}>
                            <span style={{ 
                              fontSize: "0.85em", 
                              color: "#666",
                              flexShrink: 0,
                              minWidth: "20px"
                            }}>
                              {idx + 1}.
                            </span>
                            <code style={{ 
                              fontSize: "0.9em",
                              fontWeight: selectedFile === file ? "bold" : "normal",
                              color: selectedFile === file ? "#1976d2" : "#333"
                            }}>
                              {fileName}
                            </code>
                            <span style={{ 
                              fontSize: "0.75em", 
                              color: "#999",
                              marginLeft: "auto",
                              fontStyle: "italic"
                            }}>
                              {fileDir}
                            </span>
                          </div>
                        </li>
                      );
                    })}
                  </ul>
                </div>
              </div>
            )}

            {selectedFile && (
              <div className="file-viewer" style={{ 
                marginTop: "20px", 
                border: "1px solid #ddd", 
                borderRadius: "4px",
                padding: "15px",
                backgroundColor: "#fafafa"
              }}>
                <div style={{ 
                  display: "flex", 
                  justifyContent: "space-between", 
                  alignItems: "center",
                  marginBottom: "10px",
                  paddingBottom: "10px",
                  borderBottom: "1px solid #ddd"
                }}>
                  <h4 style={{ margin: 0 }}>File: {fileName}</h4>
                  <div style={{ display: "flex", gap: "10px" }}>
                    {fileContents && (
                      <button 
                        onClick={() => {
                          navigator.clipboard.writeText(fileContents).then(() => {
                            alert("File contents copied to clipboard!");
                          }).catch(() => {
                            alert("Failed to copy to clipboard");
                          });
                        }}
                        style={{
                          padding: "5px 10px",
                          backgroundColor: "#4caf50",
                          color: "white",
                          border: "none",
                          borderRadius: "4px",
                          cursor: "pointer"
                        }}
                      >
                        Copy
                      </button>
                    )}
                    <button 
                      onClick={() => {
                        setSelectedFile(null);
                        setFileContents(null);
                      }}
                      style={{
                        padding: "5px 10px",
                        backgroundColor: "#f44336",
                        color: "white",
                        border: "none",
                        borderRadius: "4px",
                        cursor: "pointer"
                      }}
                    >
                      Close
                    </button>
                  </div>
                </div>
                {loadingFile ? (
                  <p>Loading file...</p>
                ) : fileContents !== null ? (
                  <pre style={{
                    backgroundColor: "#fff",
                    padding: "15px",
                    borderRadius: "4px",
                    overflow: "auto",
                    maxHeight: "600px",
                    fontSize: "0.85em",
                    lineHeight: "1.4",
                    border: "1px solid #ddd",
                    whiteSpace: "pre-wrap",
                    wordWrap: "break-word",
                    fontFamily: "monospace"
                  }}>
                    <code>{fileContents}</code>
                  </pre>
                ) : (
                  <p>Failed to load file contents</p>
                )}
              </div>
            )}
          </div>
        )}
      </div>
    </section>
  );
};
