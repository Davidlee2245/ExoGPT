import React, { useState } from "react";

export const Step5Panel: React.FC = () => {
  const [scoresJson, setScoresJson] = useState("./scores/pd_l1_binder_scores.json");
  const [disease, setDisease] = useState("Metastatic melanoma");
  const [biofluid, setBiofluid] = useState("Plasma");
  const [target, setTarget] = useState("PD-L1");
  const [topN, setTopN] = useState("5");

  const cli = [
    "python -m exo_gpt.step5_wetlab_protocol_designer",
    `  --scores "${scoresJson}"`,
    `  --disease "${disease}"`,
    `  --biofluid "${biofluid}"`,
    `  --target "${target}"`,
    `  --top_n ${topN}`,
    `  --out_md "./protocols/pd_l1_mm_exosome_protocol.md"`,
    `  --out_json "./protocols/pd_l1_mm_exosome_protocol.json"`
  ].join(" \\\n");

  return (
    <section className="panel">
      <header className="panel-header">
        <h2>Step 5 · Wet-lab Protocol Designer</h2>
        <p>
          Generate comprehensive wet-lab validation protocols for top-ranked binders: expression (E. coli BL21), 
          purification (Ni-NTA → SEC), binding assays (BLI/SPR, ELISA), and exosome capture experiments.
        </p>
      </header>
      <div className="panel-body">
        <div className="panel-grid">
          <div className="panel-card">
            <h3>Inputs</h3>
            <label className="field">
              <span className="field-label">Scores JSON</span>
              <input 
                value={scoresJson} 
                onChange={(e) => setScoresJson(e.target.value)}
                placeholder="./scores/..._binder_scores.json"
              />
            </label>
            <label className="field">
              <span className="field-label">Disease</span>
              <input 
                value={disease} 
                onChange={(e) => setDisease(e.target.value)}
                placeholder="Metastatic melanoma"
              />
            </label>
            <label className="field">
              <span className="field-label">Biofluid</span>
              <input 
                value={biofluid} 
                onChange={(e) => setBiofluid(e.target.value)}
                placeholder="Plasma"
              />
            </label>
            <label className="field">
              <span className="field-label">Target</span>
              <input 
                value={target} 
                onChange={(e) => setTarget(e.target.value)}
                placeholder="PD-L1"
              />
            </label>
            <label className="field">
              <span className="field-label">Top N Binders</span>
              <input 
                type="number"
                value={topN} 
                onChange={(e) => setTopN(e.target.value)}
                placeholder="5"
              />
            </label>
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

