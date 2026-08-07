import sqlite3
import os

DB_PATH = 'genomic_knowledgebase.db'

def init_db():
    if os.path.exists(DB_PATH):
        os.remove(DB_PATH)
        
    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()

    # 1. CPIC & Clinical Pharmacogenomics Rules
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS cpic_rules (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            drug_name TEXT UNIQUE,
            gene_symbol TEXT,
            phenotype TEXT,
            cpic_level TEXT,
            clinical_status TEXT,
            recommendation TEXT,
            target_disorder TEXT,
            therapeutic_class TEXT
        )
    ''')

    # 2. Pharmacokinetic (PK) & Pharmacodynamic (PD) Annotations
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS pk_pd_annotations (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            drug_name TEXT,
            gene_symbol TEXT,
            mechanism_type TEXT,
            biological_pathway TEXT,
            clinical_effect TEXT
        )
    ''')

    # 3. Polygenic Risk Score (PRS) Variant Weights
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS prs_weights (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            trait TEXT,
            rsid TEXT,
            weight REAL,
            risk_allele TEXT
        )
    ''')

    # 4. Drug-Drug Interaction (DDI) Rules
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS ddi_rules (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            drug_a TEXT,
            drug_b TEXT,
            interaction_severity TEXT,
            mechanism TEXT,
            clinical_guidance TEXT
        )
    ''')

    # 5. Pathogenicity & ACMG Secondary Findings
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS pathogenicity_db (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            rsid TEXT,
            gene_symbol TEXT,
            clinical_significance TEXT,
            associated_condition TEXT,
            acmg_actionable INTEGER
        )
    ''')

    # 6. Condition-to-Therapy Precision Matrix (PRS & Monogenic Bridge)
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS prs_condition_therapies (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            condition_or_trait TEXT,
            condition_type TEXT,
            first_line_drug TEXT,
            alternative_drug TEXT,
            target_gene_check TEXT,
            clinical_rationale TEXT
        )
    ''')

    # --- SEED DATA ---

    cpic_data = [
        # Cardiovascular & Anticoagulation
        ('Aspirin', 'CYP2C19', 'Poor Metabolizer', 'A', 'SUITABLE', 'Standard antiplatelet therapy.', 'Cardiovascular / CAD', 'Antiplatelet'),
        ('Clopidogrel', 'CYP2C19', 'Poor Metabolizer', 'A', 'CONTRAINDICATED', 'Avoid clopidogrel due to significantly reduced active metabolite formation. Switch to prasugrel or ticagrelor.', 'Thrombosis / CAD', 'Antiplatelet'),
        ('Simvastatin', 'SLCO1B1', 'Decreased Function', 'A', 'SUITABLE', 'Limit simvastatin dose to 20mg daily or switch to rosuvastatin/pravastatin.', 'Hypercholesterolemia / CAD', 'HMG-CoA Reductase Inhibitor'),
        ('Warfarin', 'CYP2C9', 'Poor Metabolizer', 'A', 'HIGH_RISK', 'Reduce initial dose by 50-80% due to severely reduced clearance.', 'Thromboembolism / Stroke', 'Anticoagulant'),
        
        # Psychiatric & Neurological (Mood & Psychotic Disorders)
        ('Aripiprazole', 'CYP2D6', 'Poor Metabolizer', 'A', 'HIGH_RISK', 'Reduce initial dose by 50% due to impaired clearance and elevated plasma levels.', 'Schizoaffective Disorder / Schizophrenia', 'Atypical Antipsychotic'),
        ('Risperidone', 'CYP2D6', 'Poor Metabolizer', 'A', 'HIGH_RISK', 'Titrate slowly or reduce dose by 50%; monitor for extrapyramidal symptoms.', 'Schizoaffective Disorder / Bipolar', 'Atypical Antipsychotic'),
        ('Clozapine', 'CYP1A2', 'Poor Metabolizer', 'A', 'SUITABLE', 'Monitor trough serum concentrations; lower maintenance doses required.', 'Treatment-Resistant Schizoaffective Disorder', 'Atypical Antipsychotic'),
        ('Carbamazepine', 'HLA-B*15:02', 'Positive', 'A', 'CONTRAINDICATED', 'Avoid due to high risk of Stevens-Johnson syndrome (SJS) and toxic epidermal necrolysis (TEN). Switch to Valproate or Lamotrigine.', 'Bipolar Disorder / Mood Stabilizer', 'Anticonvulsant / Mood Stabilizer'),
        ('Escitalopram', 'CYP2C19', 'Poor Metabolizer', 'A', 'CONTRAINDICATED', 'Reduce starting dose by 50% or select alternative drug not predominant on CYP2C19.', 'Major Depressive Disorder / MDD', 'SSRI Antidepressant'),
        ('Sertraline', 'CYP2C19', 'Poor Metabolizer', 'A', 'SUITABLE', 'Consider 50% dose reduction if co-administered with CYP2D6 inhibitors.', 'Major Depressive Disorder / MDD', 'SSRI Antidepressant'),
        
        # Immunology & Oncology
        ('Tacrolimus', 'CYP3A5', 'Poor Metabolizer', 'A', 'SUITABLE', 'Standard starting dose required for non-expressers.', 'Organ Transplant', 'Immunosuppressant'),
        ('Fluorouracil', 'DPYD', 'Poor Metabolizer', 'A', 'CONTRAINDICATED', 'Avoid use due to severe, potentially fatal toxicity.', 'Oncology', 'Antimetabolite')
    ]
    cursor.executemany("INSERT OR REPLACE INTO cpic_rules (drug_name, gene_symbol, phenotype, cpic_level, clinical_status, recommendation, target_disorder, therapeutic_class) VALUES (?,?,?,?,?,?,?,?)", cpic_data)

    pk_pd_data = [
        ('Warfarin', 'VKORC1', 'PD', 'Vitamin K Cycle / Clotting Cascade', 'Increased Sensitivity: Lower dosage required to achieve therapeutic INR target.'),
        ('Warfarin', 'CYP2C9', 'PK', 'Hepatic Phase I Clearance', 'Decreased Metabolism: Extended drug half-life and elevated systemic exposure.'),
        ('Clopidogrel', 'CYP2C19', 'PK', 'Hepatic Bioactivation', 'Impaired Conversion: Prodrug cannot be effectively converted to active thiol metabolite.'),
        ('Simvastatin', 'SLCO1B1', 'PK', 'OATP1B1 Hepatic Uptake', 'Transporter Deficiency: Decreased hepatic clearance, increasing risk of statin-induced myopathy.'),
        ('Aripiprazole', 'CYP2D6', 'PK', 'Hepatic Hydroxylation', 'Reduced Elimination: Drug accumulation increases sedation and extrapyramidal risk.'),
        ('Carbamazepine', 'HLA-B*15:02', 'PD', 'Immune System / MHC Class I Presentation', 'Immune Cytotoxicity: Direct activation of cytotoxic T-cells causing cutaneous necrosis.')
    ]
    cursor.executemany("INSERT OR REPLACE INTO pk_pd_annotations (drug_name, gene_symbol, mechanism_type, biological_pathway, clinical_effect) VALUES (?,?,?,?,?)", pk_pd_data)

    prs_data = [
        # Psychotic & Mood Traits
        ('Schizoaffective Disorder', 'rs1344706', 0.52, 'T'),
        ('Schizoaffective Disorder', 'rs1024582', 0.48, 'C'),
        ('Schizoaffective Disorder', 'rs1625579', 0.41, 'T'),
        ('Bipolar Disorder', 'rs10994336', 0.38, 'A'),
        ('Major Depressive Disorder', 'rs2522833', 0.35, 'C'),
        
        # Polygenic Chronic Conditions
        ('Coronary Artery Disease', 'rs10757278', 0.42, 'G'),
        ('Coronary Artery Disease', 'rs1333049', 0.35, 'C'),
        ('Hypercholesterolemia', 'rs629301', 0.28, 'T'),
        ('Type 2 Diabetes', 'rs7903146', 0.28, 'T'),
        ('Alzheimers Disease', 'rs429358', 0.85, 'C'),
        ('Hypertension', 'rs5186', 0.31, 'A'),
        ('Breast Cancer', 'rs11571833', 0.45, 'T')
    ]
    cursor.executemany("INSERT OR REPLACE INTO prs_weights (trait, rsid, weight, risk_allele) VALUES (?,?,?,?)", prs_data)

    ddi_data = [
        ('Clopidogrel', 'Omeprazole', 'Contraindicated', 'CYP2C19 Competitive Inhibition', 'Omeprazole inhibits CYP2C19, preventing Clopidogrel activation. Use Pantoprazole instead.'),
        ('Warfarin', 'Amiodarone', 'Major', 'CYP2C9 & CYP3A4 Inhibition', 'Amiodarone significantly increases Warfarin concentrations. Reduce Warfarin dose by 30-50%.'),
        ('Aripiprazole', 'Fluoxetine', 'Major', 'CYP2D6 Strong Inhibition', 'Fluoxetine doubles Aripiprazole exposure. Reduce Aripiprazole dose by 50%.')
    ]
    cursor.executemany("INSERT OR REPLACE INTO ddi_rules (drug_a, drug_b, interaction_severity, mechanism, clinical_guidance) VALUES (?,?,?,?,?)", ddi_data)

    path_data = [
        ('rs6025', 'F5', 'Pathogenic', 'Factor V Leiden Thrombophilia', 1),
        ('rs1801133', 'MTHFR', 'Likely_Pathogenic', 'Hyperhomocysteinemia', 0),
        ('rs28934571', 'SERPINA1', 'Pathogenic', 'Alpha-1 Antitrypsin Deficiency', 1),
        ('rs80357906', 'BRCA1', 'Pathogenic', 'Hereditary Breast and Ovarian Cancer', 1),
        ('rs4680', 'COMT', 'Pathogenic', 'Altered Prefrontal Dopamine Clearance', 0)
    ]
    cursor.executemany("INSERT OR REPLACE INTO pathogenicity_db (rsid, gene_symbol, clinical_significance, associated_condition, acmg_actionable) VALUES (?,?,?,?,?)", path_data)

    prs_therapy_data = [
        ('Schizoaffective Disorder', 'PRS_TRAIT', 'Aripiprazole', 'Clozapine / Olanzapine', 'CYP2D6', 'First-line atypical antipsychotic for elevated schizoaffective polygenic risk. Check CYP2D6 metabolizer status.'),
        ('Bipolar Disorder', 'PRS_TRAIT', 'Carbamazepine', 'Valproate / Lamotrigine', 'HLA-B*15:02', 'Mood stabilizer indicated. Contraindicated if HLA-B*15:02 positive due to SJS/TEN risk.'),
        ('Major Depressive Disorder', 'PRS_TRAIT', 'Escitalopram', 'Sertraline / Venlafaxine', 'CYP2C19', 'First-line SSRI for depressive burden. Re-assign if CYP2C19 poor metabolizer.'),
        ('Coronary Artery Disease', 'PRS_TRAIT', 'Simvastatin', 'Rosuvastatin / Atorvastatin', 'SLCO1B1', 'Primary prevention for elevated CAD burden. Adjust dose or switch to Rosuvastatin if SLCO1B1 impaired.'),
        ('Factor V Leiden Thrombophilia', 'ACMG_VARIANT', 'Warfarin', 'Direct Oral Anticoagulant (DOAC)', 'CYP2C9 / VKORC1', 'Anticoagulation indicated for thrombotic risk. Reduce Warfarin or use DOAC if CYP2C9 poor metabolizer.')
    ]
    cursor.executemany("INSERT OR REPLACE INTO prs_condition_therapies (condition_or_trait, condition_type, first_line_drug, alternative_drug, target_gene_check, clinical_rationale) VALUES (?,?,?,?,?,?)", prs_therapy_data)

    conn.commit()
    conn.close()
    print("✅ Master genomic knowledgebase updated with full psychiatric, cardiovascular, and polygenic rules.")

if __name__ == "__main__":
    init_db()
