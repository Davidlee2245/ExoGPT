import React, { useState } from "react";
import { Sidebar, Step1SubPanel, StepId } from "./components/Sidebar";
import { Step1Panel } from "./components/Step1Panel";
import { PublicationUploadPanel } from "./components/PublicationUploadPanel";
import { Step2Panel } from "./components/Step2Panel";
import { Step3Panel } from "./components/Step3Panel";
import { Step3AntibodyPanel } from "./components/Step3AntibodyPanel";
import { Step4Panel } from "./components/Step4Panel";
import { Step5Panel } from "./components/Step5Panel";
import { MrnaStep1Panel } from "./components/MrnaStep1Panel";
import { MrnaStep2Panel } from "./components/MrnaStep2Panel";
import { ProteinAgentPanel } from "./components/ProteinAgentPanel";
import { MrnaAgentPanel } from "./components/MrnaAgentPanel";

const App: React.FC = () => {
  const [activeStep, setActiveStep] = useState<StepId>("step1");
  const [activeSubPanel, setActiveSubPanel] = useState<Step1SubPanel>("main");

  const handleStepChange = (step: StepId) => {
    setActiveStep(step);
    // Reset sub-panel when switching away from step1 or mrna-step1
    if (step !== "step1" && step !== "mrna-step1") {
      setActiveSubPanel("main");
    }
  };

  const renderPanel = () => {
    switch (activeStep) {
      case "protein-agent":
        return <ProteinAgentPanel />;
      case "step1":
        if (activeSubPanel === "publication-upload") {
          return <PublicationUploadPanel dataType="protein" />;
        }
        return <Step1Panel />;
      case "step2":
        return <Step2Panel />;
      case "step3":
      case "step3-1":
        return <Step3AntibodyPanel />;
      case "step3-2":
        return <Step3Panel />;
      case "step4":
        return <Step4Panel />;
      case "step5":
        return <Step5Panel />;
      case "mrna-agent":
        return <MrnaAgentPanel />;
      case "mrna-step1":
        if (activeSubPanel === "publication-upload") {
          return <PublicationUploadPanel dataType="mrna" />;
        }
        return <MrnaStep1Panel />;
      case "mrna-step2":
        return <MrnaStep2Panel />;
      default:
        return null;
    }
  };

  return (
    <div className="app-root">
      <Sidebar 
        activeStep={activeStep} 
        activeSubPanel={activeSubPanel}
        onChange={handleStepChange}
        onSubPanelChange={setActiveSubPanel}
      />
      <main className="app-main">{renderPanel()}</main>
    </div>
  );
};

export default App;

