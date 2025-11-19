import React, { useState } from "react";

export const Step2Panel: React.FC = () => {
  const [step1Json, setStep1Json] = useState("./scores/metastatic_melanoma_plasma_biomarkers.json");

  const cli = [
    "python -m exo_gpt.step2_target_epitope_curator",
    `  --step1_json "${step1Json}"`,
    `  --out_json "./epitopes/mm_plasma_epitopes.json"`,
    `  --out_md "./epitopes/mm_plasma_epitopes.md"`
  ].join(" \\\n");

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
              <span className="field-label">Step 1 JSON Result</span>
              <input 
                value={step1Json} 
                onChange={(e) => setStep1Json(e.target.value)}
                placeholder="./scores/..._biomarkers.json"
              />
            </label>
            <p className="field-help">
              Uses structure mapping from <code>./data/structure_mapping.csv</code>
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

