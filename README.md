# Genetic-vcf-Pharmacokinetic-Pharmacodynamic-Polygenic-Risk-Score-Integration-tool

> **CLI Bioinformatics Engine:** An end-to-end precision medicine pipeline integrating VCF variant parsing, CPIC pharmacogenomics (PGx), PK/PD mechanics, Polygenic Risk Scoring (PRS), ACMG pathogenicity, and automated risk-stratified drug selection across psychiatric, cardiovascular, and metabolic disorders.

---

## ⚠️ Important Legal & Regulatory Disclaimer (Research Use Only)

**THIS SOFTWARE IS FOR RESEARCH, EDUCATION, AND INFORMATIONAL PURPOSES ONLY.**

* **Not FDA Approved:** This tool, its underlying algorithms, and its outputs have **not** been cleared, reviewed, or approved by the United States Food and Drug Administration (FDA) or any other regulatory health authority. 
* **Not Medical Advice:** This software does **not** provide medical diagnoses, clinical treatment recommendations, or actionable healthcare guidance. It does not replace professional medical judgment, genetic counseling, or clinical diagnostics.
* **No FDA Violation Intent:** By utilizing, hosting, or distributing this software, you acknowledge that it is strictly intended as a bioinformatics research pipeline for academic, exploratory, and computational experimentation. It must **not** be used for direct-to-consumer clinical genetic testing, clinical diagnosis, or medical decision-making without formal regulatory validation and oversight by a licensed medical professional.

---

## 🧬 Pipeline Architecture & Features

The pipeline processes raw VCF files through an automated 7-engine workflow:

1. **Pharmacogenomic (PGx) Phenotyping Engine:** Evaluates star-allele metabolizer profiles across major Phase I/II enzymes, transporters, and immune markers (`CYP2D6`, `CYP2C9`, `CYP2C19`, `CYP3A5`, `SLCO1B1`, `DPYD`, `HLA-B*15:02`).
2. **Polygenic Risk Scoring (PRS) Engine:** Calculates percentile risk scores across complex polygenic traits:
   * **Psychiatric / Neurological:** Schizoaffective Disorder, Bipolar Disorder, Major Depressive Disorder, Alzheimer's Disease.
   * **Cardiovascular & Metabolic:** Coronary Artery Disease (CAD), Primary Hypertension, Hypercholesterolemia, Type 2 Diabetes.
   * **Oncology:** Hereditary Breast Cancer risk modifiers.
3. **PK / PD Mechanism Mapper:** Annotates bioactivation pathways, hepatic clearance kinetics, receptor sensitivity shifts (`VKORC1`, `ADRB1`), and immune-mediated cytotoxicity paths (`HLA-B*15:02` SJS/TEN risk).
4. **CPIC & Western Medicine Matrix:** Implements evidence-based clinical dosing guidelines, contraindications, and dosing adjustments.
5. **Drug-Drug Interaction (DDI) Screener:** Flags enzyme competition and dangerous co-prescriptions (e.g., Aripiprazole + Fluoxetine, Clopidogrel + Omeprazole, Warfarin + Amiodarone).
6. **ACMG Pathogenicity & Secondary Findings:** Detects actionable monogenic variants (`F5` Factor V Leiden, `BRCA1`, `SERPINA1`, `MTHFR`, `COMT`).
7. **Polygenic Risk & Disease-Targeted Therapy Bridge:** Automatically cross-references elevated polygenic burden (PRS) or pathogenic variants (ACMG) directly against individual PGx metabolism to recommend safe first-line or reassigned alternative therapies.

---

## ⚙️ Quick Start Guide

### 1. Repository Setup & Dependencies

Clone the repository and install required dependencies:

```bash
git clone [https://github.com/operationgenetics/Genetic-vcf-Pharmacokinetic-Pharmacodynamic-Polygenic-Risk-Score-Integration-tool.-.git](https://github.com/operationgenetics/Genetic-vcf-Pharmacokinetic-Pharmacodynamic-Polygenic-Risk-Score-Integration-tool.-.git)
cd Genetic-vcf-Pharmacokinetic-Pharmacodynamic-Polygenic-Risk-Score-Integration-tool.-
pip install -r requirements.txt
