import React, { useState } from "react";

export const Step4Panel: React.FC = () => {
  const [afOutDir, setAfOutDir] = useState("./af_out/pd_l1");
  const [offtargetDir, setOfftargetDir] = useState("./offtarget");

  const cli = [
    "python -m exo_gpt.step4_scoring_offtarget",
    `  --af_out "${afOutDir}"`,
    `  --offtarget "${offtargetDir}"`,
    `  --out_json "./scores/pd_l1_binder_scores.json"`,
    `  --out_md "./scores/pd_l1_binder_scores_topN.md"`
  ].join(" \\\n");

  return (
    <section className="panel">
      <header className="panel-header">
        <h2>Step 4 · In-silico Scoring & Off-target Filtering</h2>
        <p>
          Parse AlphaFold output metrics (ipTM, ipAE, pLDDT), evaluate off-target binding, apply pass/fail criteria, 
          and rank binders by composite score.
        </p>
      </header>
      <div className="panel-body">
        <div className="panel-grid">
          <div className="panel-card">
            <h3>Inputs</h3>
            <label className="field">
              <span className="field-label">AlphaFold Output Directory</span>
              <input 
                value={afOutDir} 
                onChange={(e) => setAfOutDir(e.target.value)}
                placeholder="./af_out/pd_l1"
              />
            </label>
            <label className="field">
              <span className="field-label">Off-target Directory</span>
              <input 
                value={offtargetDir} 
                onChange={(e) => setOfftargetDir(e.target.value)}
                placeholder="./offtarget"
              />
            </label>
            <p className="field-help">
              Pass criteria: ipTM_on ≥ 0.80, ipAE_on ≤ 9, pLDDT ≥ 80, off-target max ipTM ≤ 0.70
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

