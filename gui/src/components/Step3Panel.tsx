import React, { useState } from "react";

export const Step3Panel: React.FC = () => {
  const [configPath, setConfigPath] = useState("./workflows/example_pd_l1_config.json");

  const cli = [
    "python -m exo_gpt.step3_nanobinder_orchestrator",
    `  --config "${configPath}"`,
    `  --out_plan "./workflows/design_plan.json"`
  ].join(" \\\n");

  return (
    <section className="panel">
      <header className="panel-header">
        <h2>Step 3 · Nanobinder Design Orchestrator</h2>
        <p>
          Orchestrates RFdiffusion backbone generation, ProteinMPNN sequence design, and AlphaFold-Multimer (ColabFold) 
          structure prediction. Generates design plan or executes full pipeline.
        </p>
      </header>
      <div className="panel-body">
        <div className="panel-grid">
          <div className="panel-card">
            <h3>Inputs</h3>
            <label className="field">
              <span className="field-label">Config JSON Path</span>
              <input 
                value={configPath} 
                onChange={(e) => setConfigPath(e.target.value)}
                placeholder="./workflows/example_pd_l1_config.json"
              />
            </label>
            <p className="field-help">
              Config specifies: target PDB, epitope residues, binder length, number of designs, output directories.
            </p>
          </div>
          <div className="panel-card">
            <h3>Command</h3>
            <p className="field-help">Copy-paste into your Exosome-GPT Python environment.</p>
            <pre className="code-block">
              <code>{cli}</code>
            </pre>
          </div>
        </div>
      </div>
    </section>
  );
};

