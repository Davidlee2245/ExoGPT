import React from "react";

export type StepId = "step1" | "step2" | "step3" | "step4" | "step5";
export type Step1SubPanel = "main" | "publication-upload";

interface SidebarProps {
  activeStep: StepId;
  activeSubPanel?: Step1SubPanel;
  onChange: (step: StepId) => void;
  onSubPanelChange?: (subPanel: Step1SubPanel) => void;
}

export const Sidebar: React.FC<SidebarProps> = ({ 
  activeStep, 
  activeSubPanel = "main",
  onChange, 
  onSubPanelChange 
}) => {
  const steps = [
    { id: "step1" as StepId, label: "Step 1", desc: "EV Biomarker Finder" },
    { id: "step2" as StepId, label: "Step 2", desc: "Target & Epitope Curator" },
    { id: "step3" as StepId, label: "Step 3", desc: "Nanobinder Orchestrator" },
    { id: "step4" as StepId, label: "Step 4", desc: "Scoring & Off-target" },
    { id: "step5" as StepId, label: "Step 5", desc: "Wet-lab Protocol" },
  ];

  const step1SubPanels: { id: Step1SubPanel; label: string; desc: string }[] = [
    { id: "main", label: "Database Search", desc: "Search existing databases" },
    { id: "publication-upload", label: "Publication Upload", desc: "RAG-based biomarker search" },
  ];

  return (
    <aside className="sidebar">
      <div className="sidebar-header">
        <h1 className="brand-title">Exosome-GPT</h1>
        <p className="brand-subtitle">Gemini-style nanobinder designer</p>
      </div>
      <nav className="sidebar-nav">
        {steps.map((step) => (
          <div key={step.id}>
            <button
              className={`sidebar-item ${activeStep === step.id ? "active" : ""}`}
              onClick={() => onChange(step.id)}
            >
              <span className="sidebar-item-label">{step.label}</span>
              <span className="sidebar-item-desc">{step.desc}</span>
            </button>
            {step.id === "step1" && activeStep === "step1" && onSubPanelChange && (
              <div className="sidebar-subnav">
                {step1SubPanels.map((subPanel) => (
                  <button
                    key={subPanel.id}
                    className={`sidebar-subitem ${activeSubPanel === subPanel.id ? "active" : ""}`}
                    onClick={() => onSubPanelChange(subPanel.id)}
                  >
                    <span className="sidebar-subitem-label">{subPanel.label}</span>
                    <span className="sidebar-subitem-desc">{subPanel.desc}</span>
                  </button>
                ))}
              </div>
            )}
          </div>
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

