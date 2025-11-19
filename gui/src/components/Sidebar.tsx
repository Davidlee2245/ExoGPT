import React from "react";

export type StepId = "step1" | "step2" | "step3" | "step4" | "step5";

interface SidebarProps {
  activeStep: StepId;
  onChange: (step: StepId) => void;
}

export const Sidebar: React.FC<SidebarProps> = ({ activeStep, onChange }) => {
  const steps = [
    { id: "step1" as StepId, label: "Step 1", desc: "EV Biomarker Finder" },
    { id: "step2" as StepId, label: "Step 2", desc: "Target & Epitope Curator" },
    { id: "step3" as StepId, label: "Step 3", desc: "Nanobinder Orchestrator" },
    { id: "step4" as StepId, label: "Step 4", desc: "Scoring & Off-target" },
    { id: "step5" as StepId, label: "Step 5", desc: "Wet-lab Protocol" },
  ];

  return (
    <aside className="sidebar">
      <div className="sidebar-header">
        <h1 className="brand-title">Exosome-GPT</h1>
        <p className="brand-subtitle">Gemini-style nanobinder designer</p>
      </div>
      <nav className="sidebar-nav">
        {steps.map((step) => (
          <button
            key={step.id}
            className={`sidebar-item ${activeStep === step.id ? "active" : ""}`}
            onClick={() => onChange(step.id)}
          >
            <span className="sidebar-item-label">{step.label}</span>
            <span className="sidebar-item-desc">{step.desc}</span>
          </button>
        ))}
      </nav>
      <div className="sidebar-footer">
        <span className="sidebar-footer-title">Pipeline</span>
        <span className="sidebar-footer-text">
          Disease → Biomarker → Target → Design → Score → Protocol
        </span>
      </div>
    </aside>
  );
};

