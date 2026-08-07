import sqlite3
import os

def init_db():
    conn = sqlite3.connect("genomic_knowledgebase.db")
    cursor = conn.cursor()

    # Drop old tables if schema mismatched
    cursor.execute("DROP TABLE IF EXISTS ddi_rules")

    # 1. Pharmacogenomic Star Alleles Table
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS pgx_star_alleles (
        gene_symbol TEXT,
        star_allele TEXT,
        metabolizer_status TEXT,
        genotype_pattern TEXT,
        PRIMARY KEY (gene_symbol, star_allele)
    )
    """)

    # 2. Polygenic Risk Scores Table
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS prs_rules (
        disease_name TEXT PRIMARY KEY,
        percentile INTEGER,
        risk_category TEXT,
        base_score REAL
    )
    """)

    # 3. PK / PD Mechanisms Table
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS pk_pd_mechanisms (
        drug_name TEXT,
        gene_symbol TEXT,
        mechanism_type TEXT,
        effect_summary TEXT,
        PRIMARY KEY (drug_name, gene_symbol)
    )
    """)

    # 4. CPIC Guidelines Table
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS cpic_guidelines (
        drug_name TEXT PRIMARY KEY,
        gene_symbol TEXT,
        status TEXT,
        recommendation TEXT
    )
    """)

    # 5. Drug-Drug Interactions Table (4 columns)
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS ddi_rules (
        drug1 TEXT,
        drug2 TEXT,
        interaction_level TEXT,
        clinical_effect TEXT,
        PRIMARY KEY (drug1, drug2)
    )
    """)

    # 6. ACMG Secondary Findings Table
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS acmg_findings (
        rsid TEXT PRIMARY KEY,
        gene_symbol TEXT,
        pathogenicity TEXT,
        disease_association TEXT,
        actionable INTEGER
    )
    """)

    # 7. Disease Targeted Therapies Table
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS disease_targeted_therapies (
        condition_trait TEXT PRIMARY KEY,
        primary_drug TEXT,
        gene_checked TEXT,
        pgx_status TEXT,
        alternative_drug TEXT,
        clinical_rationale TEXT
    )
    """)

    # 8. Genome-Wide ClinVar Table
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS genome_clinvar (
        rsid TEXT PRIMARY KEY,
        gene_symbol TEXT,
        clinical_significance TEXT,
        associated_conditions TEXT,
        review_status TEXT
    )
    """)

    # Seed Data
    pgx_data = [
        ('CYP2C9', '*3', 'Poor Metabolizer', 'A/A'),
        ('CYP2C19', '*2', 'Poor Metabolizer', 'A/A'),
        ('SLCO1B1', '*5', 'Decreased Function', 'C/C'),
        ('CYP3A5', '*3', 'Non-Expresser', 'G/G'),
        ('CYP2D6', '*10', 'Poor Metabolizer', 'T/T'),
        ('VKORC1', '-1639G>A', 'High Sensitivity', 'A/A')
    ]
    cursor.executemany("INSERT OR REPLACE INTO pgx_star_alleles VALUES (?, ?, ?, ?)", pgx_data)

    prs_data = [
        ('Schizoaffective Disorder', 63, 'Moderate', 1.41),
        ('Bipolar Disorder', 17, 'Low', 0.38),
        ('Major Depressive Disorder', 15, 'Low', 0.35),
        ('Alzheimers Disease', 38, 'Moderate', 0.85),
        ('Coronary Artery Disease', 34, 'Moderate', 0.77),
        ('Hypertension', 13, 'Low', 0.31),
        ('Hypercholesterolemia', 12, 'Low', 0.28),
        ('Type 2 Diabetes', 12, 'Low', 0.28),
        ('Breast Cancer', 20, 'Low', 0.45)
    ]
    cursor.executemany("INSERT OR REPLACE INTO prs_rules VALUES (?, ?, ?, ?)", prs_data)

    pkpd_data = [
        ('Warfarin', 'VKORC1', 'PD', 'Increased Sensitivity: Lower dosage required to achieve therapeutic INR target.'),
        ('Warfarin', 'CYP2C9', 'PK', 'Decreased Metabolism: Extended drug half-life and elevated systemic exposure.'),
        ('Clopidogrel', 'CYP2C19', 'PK', 'Impaired Conversion: Prodrug cannot be effectively converted to active thiol metabolite.'),
        ('Simvastatin', 'SLCO1B1', 'PK', 'Transporter Deficiency: Decreased hepatic clearance, increasing risk of statin-induced myopathy.'),
        ('Aripiprazole', 'CYP2D6', 'PK', 'Reduced Elimination: Drug accumulation increases sedation and extrapyramidal risk.'),
        ('Carbamazepine', 'HLA-B*15:02', 'PD', 'Immune Cytotoxicity: Direct activation of cytotoxic T-cells causing cutaneous necrosis.')
    ]
    cursor.executemany("INSERT OR REPLACE INTO pk_pd_mechanisms VALUES (?, ?, ?, ?)", pkpd_data)

    cpic_data = [
        ('Aspirin', 'CYP2C19', 'SUITABLE', 'Standard antiplatelet therapy.'),
        ('Clopidogrel', 'CYP2C19', 'CONTRAINDICATED', 'Avoid clopidogrel due to significantly reduced active metabolite formation. Switch to prasugrel or ticagrelor.'),
        ('Simvastatin', 'SLCO1B1', 'SUITABLE', 'Limit simvastatin dose to 20mg daily or switch to rosuvastatin/pravastatin.'),
        ('Warfarin', 'CYP2C9', 'HIGH_RISK', 'Reduce initial dose by 50-80% due to severely reduced clearance.'),
        ('Aripiprazole', 'CYP2D6', 'HIGH_RISK', 'Reduce initial dose by 50% due to impaired clearance and elevated plasma levels.'),
        ('Risperidone', 'CYP2D6', 'HIGH_RISK', 'Titrate slowly or reduce dose by 50%; monitor for extrapyramidal symptoms.'),
        ('Clozapine', 'CYP1A2', 'SUITABLE', 'Monitor trough serum concentrations; lower maintenance doses required.'),
        ('Carbamazepine', 'HLA-B*15:02', 'CONTRAINDICATED', 'Avoid due to high risk of Stevens-Johnson syndrome (SJS) and toxic epidermal necrolysis (TEN). Switch to Valproate or Lamotrigine.'),
        ('Escitalopram', 'CYP2C19', 'CONTRAINDICATED', 'Reduce starting dose by 50% or select alternative drug not predominant on CYP2C19.'),
        ('Sertraline', 'CYP2C19', 'SUITABLE', 'Consider 50% dose reduction if co-administered with CYP2D6 inhibitors.'),
        ('Tacrolimus', 'CYP3A5', 'SUITABLE', 'Standard starting dose required for non-expressers.'),
        ('Fluorouracil', 'DPYD', 'CONTRAINDICATED', 'Avoid use due to severe, potentially fatal toxicity.')
    ]
    cursor.executemany("INSERT OR REPLACE INTO cpic_guidelines VALUES (?, ?, ?, ?)", cpic_data)

    ddi_data = [
        ('Clopidogrel', 'Omeprazole', 'Contraindicated', 'Omeprazole inhibits CYP2C19, preventing Clopidogrel activation. Use Pantoprazole instead.'),
        ('Warfarin', 'Amiodarone', 'Major', 'Amiodarone significantly increases Warfarin concentrations. Reduce Warfarin dose by 30-50%.'),
        ('Aripiprazole', 'Fluoxetine', 'Major', 'Fluoxetine doubles Aripiprazole exposure. Reduce Aripiprazole dose by 50%.')
    ]
    cursor.executemany("INSERT OR REPLACE INTO ddi_rules VALUES (?, ?, ?, ?)", ddi_data)

    acmg_data = [
        ('rs6025', 'F5', 'Pathogenic', 'Factor V Leiden Thrombophilia', 1),
        ('rs1801133', 'MTHFR', 'Likely_Pathogenic', 'Hyperhomocysteinemia', 0),
        ('rs28934571', 'SERPINA1', 'Pathogenic', 'Alpha-1 Antitrypsin Deficiency', 1),
        ('rs80357906', 'BRCA1', 'Pathogenic', 'Hereditary Breast and Ovarian Cancer', 1),
        ('rs4680', 'COMT', 'Pathogenic', 'Altered Prefrontal Dopamine Clearance', 0)
    ]
    cursor.executemany("INSERT OR REPLACE INTO acmg_findings VALUES (?, ?, ?, ?, ?)", acmg_data)

    therapy_data = [
        ('Schizoaffective Disorder', 'Aripiprazole', 'CYP2D6', 'HIGH_RISK', 'Clozapine / Olanzapine', 'First-line atypical antipsychotic for elevated schizoaffective polygenic risk. Check CYP2D6 metabolizer status.'),
        ('Coronary Artery Disease', 'Simvastatin', 'SLCO1B1', 'SUITABLE', 'Rosuvastatin', 'Primary prevention for elevated CAD burden. Adjust dose or switch to Rosuvastatin if SLCO1B1 impaired.'),
        ('Factor V Leiden Thrombophilia', 'Warfarin', 'CYP2C9', 'HIGH_RISK', 'Direct Oral Anticoagulant (DOAC)', 'Anticoagulation indicated for thrombotic risk. Reduce Warfarin or use DOAC if CYP2C9 poor metabolizer.')
    ]
    cursor.executemany("INSERT OR REPLACE INTO disease_targeted_therapies VALUES (?, ?, ?, ?, ?, ?)", therapy_data)

    clinvar_data = [
        ('rs6025', 'F5', 'Pathogenic', 'Factor V Leiden Thrombophilia', 'classified_by_single_submitter'),
        ('rs1801133', 'MTHFR', 'Likely_Pathogenic', 'Hyperhomocysteinemia', 'criteria_provided'),
        ('rs28934571', 'SERPINA1', 'Pathogenic', 'Alpha-1 Antitrypsin Deficiency', 'reviewed_by_expert_panel'),
        ('rs80357906', 'BRCA1', 'Pathogenic', 'Hereditary Breast and Ovarian Cancer', 'reviewed_by_expert_panel'),
        ('rs4680', 'COMT', 'Pathogenic', 'Altered Prefrontal Dopamine Clearance', 'criteria_provided')
    ]
    cursor.executemany("INSERT OR REPLACE INTO genome_clinvar VALUES (?, ?, ?, ?, ?)", clinvar_data)

    conn.commit()
    conn.close()

    # Ensure data directory and sample VCF exist
    os.makedirs("data", exist_ok=True)
    os.makedirs("results", exist_ok=True)
    vcf_path = "data/psychiatric_patient.vcf"
    if not os.path.exists(vcf_path):
        with open(vcf_path, "w") as f:
            f.write("##fileformat=VCFv4.2\n")
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPSYCH_PATIENT_01\n")
            f.write("chr22\t42522500\trs1065852\tG\tA\t100\tPASS\t.\tGT\t1/1\n")
            f.write("chr10\t96522500\trs1057910\tC\tT\t100\tPASS\t.\tGT\t0/0\n")
            f.write("chr1\t169519049\trs6025\tC\tT\t100\tPASS\t.\tGT\t0/1\n")

    print("[✔] Master database successfully initialized with all 8 tables and sample data.")

if __name__ == "__main__":
    init_db()
