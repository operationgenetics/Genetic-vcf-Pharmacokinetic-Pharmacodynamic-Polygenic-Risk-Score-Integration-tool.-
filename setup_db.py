import sqlite3
import os

DB_PATH = 'genomic_knowledgebase.db'

def init_db():
    if os.path.exists(DB_PATH):
        os.remove(DB_PATH)
        
    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()

    # 1. CPIC & Western Medicine Base Rules
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS cpic_rules (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            drug_name TEXT,
            gene_symbol TEXT,
            phenotype TEXT,
            cpic_level TEXT,
            clinical_status TEXT,
            recommendation TEXT,
            target_disorder TEXT,
            therapeutic_class TEXT
        )
    ''')

    # 2. PK / PD Engine Table
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS pk_pd_annotations (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            drug_name TEXT,
            gene_symbol TEXT,
            mechanism_type TEXT, -- 'PK' (Absorption/Metabolism) or 'PD' (Target/Receptor)
            biological_pathway TEXT,
            clinical_effect TEXT
        )
    ''')

    # 3. Polygenic Risk Scores (PRS) Expansion
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS prs_weights (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            trait TEXT,
            rsid TEXT,
            weight REAL,
            risk_allele TEXT
        )
    ''')

    # 4. Drug-Drug Interaction (DDI) Matrix
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS ddi_rules (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            drug_a TEXT,
            drug_b TEXT,
            interaction_severity TEXT, -- 'Contraindicated', 'Major', 'Moderate'
            mechanism TEXT,
            clinical_guidance TEXT
        )
    ''')

    # 5. Pathogenicity Engine (ClinVar / ACMG)
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS pathogenicity_db (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            rsid TEXT,
            gene_symbol TEXT,
            clinical_significance TEXT, -- 'Pathogenic', 'Likely_Pathogenic', 'VUS'
            associated_condition TEXT,
            acmg_actionable INTEGER -- 1 if ACMG Secondary Finding, 0 otherwise
        )
    ''')

    # --- SEED DATA ---

    # CPIC Baseline Data
    cpic_data = [
        ('Aspirin', 'CYP2C19', 'Poor Metabolizer', 'A', 'SUITABLE', 'Standard antiplatelet therapy.', 'Cardiovascular', 'Antiplatelet'),
        ('Clopidogrel', 'CYP2C19', 'Poor Metabolizer', 'A', 'CONTRAINDICATED', 'Avoid clopidogrel due to significantly reduced active metabolite formation. Switch to prasugrel or ticagrelor.', 'Thrombosis', 'Antiplatelet'),
        ('Escitalopram', 'CYP2C19', 'Poor Metabolizer', 'A', 'CONTRAINDICATED', 'Reduce starting dose by 50% or select alternative drug not predominant on CYP2C19.', 'Major Depressive Disorder', 'SSRI Antidepressant'),
        ('Simvastatin', 'SLCO1B1', 'Decreased Function', 'A', 'SUITABLE', 'Limit simvastatin dose to 20mg daily or switch to rosuvastatin/pravastatin.', 'Hypercholesterolemia', 'HMG-CoA Reductase Inhibitor'),
        ('Warfarin', 'CYP2C9', 'Poor Metabolizer', 'A', 'HIGH_RISK', 'Reduce initial dose by 50-80% due to severely reduced clearance.', 'Thromboembolism', 'Anticoagulant'),
        ('Tacrolimus', 'CYP3A5', 'Poor Metabolizer', 'A', 'SUITABLE', 'Standard starting dose required for non-expressers.', 'Organ Transplant', 'Immunosuppressant'),
        ('Fluorouracil', 'DPYD', 'Poor Metabolizer', 'A', 'CONTRAINDICATED', 'Avoid use due to severe, potentially fatal toxicity.', 'Oncology', 'Antimetabolite')
    ]
    cursor.executemany("INSERT INTO cpic_rules (drug_name, gene_symbol, phenotype, cpic_level, clinical_status, recommendation, target_disorder, therapeutic_class) VALUES (?,?,?,?,?,?,?,?)", cpic_data)

    # PK/PD Annotations
    pk_pd_data = [
        ('Warfarin', 'VKORC1', 'PD', 'Vitamin K Cycle / Clotting Cascade', 'Increased Sensitivity: Lower dosage required to achieve therapeutic INR target.'),
        ('Warfarin', 'CYP2C9', 'PK', 'Hepatic Phase I Clearance', 'Decreased Metabolism: Extended drug half-life and elevated systemic exposure.'),
        ('Clopidogrel', 'CYP2C19', 'PK', 'Hepatic Bioactivation', 'Impaired Conversion: Prodrug cannot be effectively converted to active thiol metabolite.'),
        ('Simvastatin', 'SLCO1B1', 'PK', 'OATP1B1 Hepatic Uptake', 'Transporter Deficiency: Decreased hepatic clearance, increasing risk of statin-induced myopathy.'),
        ('Codeine', 'CYP2D6', 'PK', 'O-demethylation to Morphine', 'Poor Metabolism: Absence of analgesic effect due to failure to produce morphine.'),
        ('Metoprolol', 'ADRB1', 'PD', 'Beta-1 Adrenergic Receptor Signaling', 'Altered Sensitivity: Enhanced blood pressure and heart rate response.')
    ]
    cursor.executemany("INSERT INTO pk_pd_annotations (drug_name, gene_symbol, mechanism_type, biological_pathway, clinical_effect) VALUES (?,?,?,?,?)", pk_pd_data)

    # Multi-Trait PRS Weights
    prs_data = [
        ('Coronary Artery Disease', 'rs10757278', 0.42, 'G'),
        ('Coronary Artery Disease', 'rs1333049', 0.35, 'C'),
        ('Hypercholesterolemia', 'rs629301', 0.28, 'T'),
        ('Type 2 Diabetes', 'rs7903146', 0.28, 'T'),
        ('Alzheimers Disease', 'rs429358', 0.85, 'C'),
        ('Hypertension', 'rs5186', 0.31, 'A'),
        ('Breast Cancer', 'rs11571833', 0.45, 'T')
    ]
    cursor.executemany("INSERT INTO prs_weights (trait, rsid, weight, risk_allele) VALUES (?,?,?,?)", prs_data)

    # DDI Rules Matrix
    ddi_data = [
        ('Clopidogrel', 'Omeprazole', 'Contraindicated', 'CYP2C19 Competitive Inhibition', 'Omeprazole inhibits CYP2C19, preventing Clopidogrel activation. Use Pantoprazole instead.'),
        ('Warfarin', 'Amiodarone', 'Major', 'CYP2C9 & CYP3A4 Inhibition', 'Amiodarone significantly increases Warfarin concentrations. Reduce Warfarin dose by 30-50%.'),
        ('Simvastatin', 'Amlodipine', 'Moderate', 'CYP3A4 Metabolic Competition', 'Increased risk of myopathy. Do not exceed 20mg/day of Simvastatin.')
    ]
    cursor.executemany("INSERT INTO ddi_rules (drug_a, drug_b, interaction_severity, mechanism, clinical_guidance) VALUES (?,?,?,?,?)", ddi_data)

    # Pathogenicity Engine Data
    path_data = [
        ('rs6025', 'F5', 'Pathogenic', 'Factor V Leiden Thrombophilia', 1),
        ('rs1801133', 'MTHFR', 'Likely_Pathogenic', 'Hyperhomocysteinemia', 0),
        ('rs28934571', 'SERPINA1', 'Pathogenic', 'Alpha-1 Antitrypsin Deficiency', 1),
        ('rs80357906', 'BRCA1', 'Pathogenic', 'Hereditary Breast and Ovarian Cancer', 1)
    ]
    cursor.executemany("INSERT INTO pathogenicity_db (rsid, gene_symbol, clinical_significance, associated_condition, acmg_actionable) VALUES (?,?,?,?,?)", path_data)

    conn.commit()
    conn.close()
    print("✅ Full genomic knowledgebase initialized with PK, PD, PRS, DDI, and Pathogenicity tables.")

if __name__ == "__main__":
    init_db()
