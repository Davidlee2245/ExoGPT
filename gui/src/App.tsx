import React, { useState } from "react";
import { Sidebar } from "./components/Sidebar";
import { Step1Panel } from "./components/Step1Panel";
import { Step2Panel } from "./components/Step2Panel";
import { Step3Panel } from "./components/Step3Panel";
import { Step4Panel } from "./components/Step4Panel";
import { Step5Panel } from "./components/Step5Panel";

export type StepId = "step1" | "step2" | "step3" | "step4" | "step5";

const App: React.FC = () => {
  const [activeStep, setActiveStep] = useState<StepId>("step1");

  const renderPanel = () => {
    switch (activeStep) {
      case "step1":
        return <Step1Panel />;
      case "step2":
        return <Step2Panel />;
      case "step3":
        return <Step3Panel />;
      case "step4":
        return <Step4Panel />;
      case "step5":
        return <Step5Panel />;
      default:
        return null;
    }
  };

  return (
    <div className="app-root">
      <Sidebar activeStep={activeStep} onChange={setActiveStep} />
      <main className="app-main">{renderPanel()}</main>
    </div>
  );
};

export default App;

