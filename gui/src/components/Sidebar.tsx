import React, { useState, useEffect } from "react";
import exosomeIcon from "../../icon/exosome.png";
import proteinIcon from "../../icon/protein.png";
import mrnaIcon from "../../icon/mRNA.png";

export type StepId = "protein-agent" | "step1" | "step2" | "step3" | "step3-1" | "step3-2" | "step4" | "step5" | "mrna-agent" | "mrna-step1" | "mrna-step2";
export type Step1SubPanel = "main" | "publication-upload";
export type AgentType = "protein" | "mrna";

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
  // Track which accordion section is expanded (only one at a time)
  const [expandedSection, setExpandedSection] = useState<"protein" | "mrna" | null>(null);

  // Auto-expand the section containing the active step
  useEffect(() => {
    const proteinSteps: StepId[] = ["protein-agent", "step1", "step2", "step3", "step3-1", "step3-2", "step4", "step5"];
    const mrnaSteps: StepId[] = ["mrna-agent", "mrna-step1", "mrna-step2"];
    
    if (proteinSteps.includes(activeStep)) {
      setExpandedSection("protein");
    } else if (mrnaSteps.includes(activeStep)) {
      setExpandedSection("mrna");
    }
  }, [activeStep]);

  const toggleSection = (section: "protein" | "mrna") => {
    setExpandedSection(expandedSection === section ? null : section);
  };
  const proteinAgent = { id: "protein-agent" as StepId, label: "Protein Agent", desc: "Chat with Protein Agent" };
  
  const proteinSteps = [
    { id: "step1" as StepId, label: "Step 1", desc: "EV Biomarker Finder" },
    { id: "step2" as StepId, label: "Step 2", desc: "Target & Epitope Curator" },
    { id: "step4" as StepId, label: "Step 4", desc: "Binding Score" },
    { id: "step5" as StepId, label: "Step 5", desc: "Protocol Generation" },
  ];

  const step3SubSteps = [
    { id: "step3-1" as StepId, label: "Step 3-1", desc: "Antibody Recommendation" },
    { id: "step3-2" as StepId, label: "Step 3-2", desc: "Nanobinder Design" },
  ];

  const mrnaAgent = { id: "mrna-agent" as StepId, label: "mRNA Agent", desc: "Chat with mRNA Agent" };
  
  const mrnaSteps = [
    { id: "mrna-step1" as StepId, label: "Step 1", desc: "EV Biomarker Finder" },
    { id: "mrna-step2" as StepId, label: "Step 2", desc: "Primer Design" },
  ];

  const proteinStep1SubPanels: { id: Step1SubPanel; label: string; desc: string }[] = [
    { id: "main", label: "Database Search (Cancer/Normal)", desc: "Search existing databases" },
    { id: "publication-upload", label: "Publication Upload", desc: "RAG-based biomarker search" },
  ];

  const mrnaStep1SubPanels: { id: Step1SubPanel; label: string; desc: string }[] = [
    { id: "main", label: "Database Search", desc: "Search existing databases" },
    { id: "publication-upload", label: "Publication Upload", desc: "RAG-based biomarker search" },
  ];

  const renderStep = (step: { id: StepId; label: string; desc: string }, agentType: AgentType) => {
    const isActive = activeStep === step.id;
    const showSubPanels = step.id === "step1" && isActive && onSubPanelChange;
    const showMrnaSubPanels = step.id === "mrna-step1" && isActive && onSubPanelChange;
    const step1SubPanels = step.id === "step1" ? proteinStep1SubPanels : mrnaStep1SubPanels;

  return (
          <div key={step.id}>
            <button
          className={`sidebar-item ${isActive ? "active" : ""}`}
              onClick={() => onChange(step.id)}
            >
              <span className="sidebar-item-label">{step.label}</span>
          {step.desc && <span className="sidebar-item-desc">{step.desc}</span>}
            </button>
        {(showSubPanels || showMrnaSubPanels) && (
              <div className="sidebar-subnav">
                {step1SubPanels.map((subPanel) => (
                  <button
                    key={subPanel.id}
                    className={`sidebar-subitem ${activeSubPanel === subPanel.id ? "active" : ""}`}
                onClick={() => onSubPanelChange?.(subPanel.id)}
                  >
                    <span className="sidebar-subitem-label">{subPanel.label}</span>
                    <span className="sidebar-subitem-desc">{subPanel.desc}</span>
                  </button>
                ))}
              </div>
            )}
          </div>
    );
  };

  const renderStep3 = () => {
    const isStep3Active = activeStep === "step3" || activeStep === "step3-1" || activeStep === "step3-2";
    const showStep3SubSteps = isStep3Active;

    return (
      <div key="step3">
        <button
          className={`sidebar-item ${isStep3Active ? "active" : ""}`}
          onClick={() => onChange("step3-1")}
        >
          <span className="sidebar-item-label">Step 3</span>
          <span className="sidebar-item-desc">Binder optimization</span>
        </button>
        {showStep3SubSteps && (
          <div className="sidebar-subnav">
            {step3SubSteps.map((subStep) => {
              const isSubStepActive = activeStep === subStep.id;
              return (
                <button
                  key={subStep.id}
                  className={`sidebar-subitem ${isSubStepActive ? "active" : ""}`}
                  onClick={() => onChange(subStep.id)}
                >
                  <span className="sidebar-subitem-label">{subStep.label}</span>
                  <span className="sidebar-subitem-desc">{subStep.desc}</span>
                </button>
              );
            })}
          </div>
        )}
      </div>
    );
  };

  return (
    <aside className="sidebar">
      <div className="sidebar-header">
        <div className="brand-title-wrapper">
          <img src={exosomeIcon} alt="Exosome" className="brand-icon" />
          <h1 className="brand-title">Exosome-GPT</h1>
        </div>
        <p className="brand-subtitle">Gemini-style nanobinder designer</p>
      </div>
      <nav className="sidebar-nav">
        <div className="sidebar-section">
          <button
            className="sidebar-section-header-button"
            onClick={() => toggleSection("protein")}
          >
            <span>Protein Agents</span>
            <svg
              className={`sidebar-chevron ${expandedSection === "protein" ? "expanded" : ""}`}
              width="12"
              height="12"
              viewBox="0 0 12 12"
              fill="none"
              xmlns="http://www.w3.org/2000/svg"
            >
              <path
                d="M3 4.5L6 7.5L9 4.5"
                stroke="currentColor"
                strokeWidth="1.5"
                strokeLinecap="round"
                strokeLinejoin="round"
              />
            </svg>
          </button>
          {expandedSection === "protein" && (
            <div className="sidebar-section-content">
              <button
                className={`sidebar-item ${activeStep === proteinAgent.id ? "active" : ""}`}
                onClick={() => onChange(proteinAgent.id)}
              >
                <img src={proteinIcon} alt="Protein" className="sidebar-item-icon" />
                <span className="sidebar-item-label">{proteinAgent.label}</span>
                {proteinAgent.desc && <span className="sidebar-item-desc">{proteinAgent.desc}</span>}
              </button>
              {proteinSteps.slice(0, 2).map((step) => renderStep(step, "protein"))}
              {renderStep3()}
              {proteinSteps.slice(2).map((step) => renderStep(step, "protein"))}
            </div>
          )}
        </div>
        <div className="sidebar-section">
          <button
            className="sidebar-section-header-button"
            onClick={() => toggleSection("mrna")}
          >
            <span>mRNA Agents</span>
            <svg
              className={`sidebar-chevron ${expandedSection === "mrna" ? "expanded" : ""}`}
              width="12"
              height="12"
              viewBox="0 0 12 12"
              fill="none"
              xmlns="http://www.w3.org/2000/svg"
            >
              <path
                d="M3 4.5L6 7.5L9 4.5"
                stroke="currentColor"
                strokeWidth="1.5"
                strokeLinecap="round"
                strokeLinejoin="round"
              />
            </svg>
          </button>
          {expandedSection === "mrna" && (
            <div className="sidebar-section-content">
              <button
                className={`sidebar-item ${activeStep === mrnaAgent.id ? "active" : ""}`}
                onClick={() => onChange(mrnaAgent.id)}
              >
                <img src={mrnaIcon} alt="mRNA" className="sidebar-item-icon" />
                <span className="sidebar-item-label">{mrnaAgent.label}</span>
                {mrnaAgent.desc && <span className="sidebar-item-desc">{mrnaAgent.desc}</span>}
              </button>
              {mrnaSteps.map((step) => renderStep(step, "mrna"))}
            </div>
          )}
        </div>
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

