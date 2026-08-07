import sqlite3

def setup_database():
    conn = sqlite3.connect("genomic_knowledgebase.db")
    cursor = conn.cursor()

    # Drop existing tables to enforce full clean schema
    tables = [
        "pgx_star_alleles", "prs_traits", "prs_weights", 
        "pk_pd_mechanisms", "cpic_guidelines", "ddi_rules", 
        "acmg_findings", "disease_targeted_therapies"
    ]
    for table in tables:
        cursor.execute(f"DROP TABLE IF EXISTS {table}")

    # 1. PGx Star Alleles with Haplotype Support
    cursor.execute("""
    CREATE TABLE pgx_star_alleles (
        gene_symbol TEXT,
        diplotype TEXT,
        metabolizer_status TEXT,
        required_variants TEXT,
        PRIMARY KEY (gene_symbol, diplotype)
    )""")

    # 2. PRS Trait Parameters (Mean & StdDev for Z-score calc)
    cursor.execute("""
    CREATE TABLE prs_traits (
        trait TEXT PRIMARY KEY,
        mean_score REAL NOT NULL,
        std_score REAL NOT NULL
    )""")

    # 3. PRS Weights (PGS Catalog style)
    cursor.execute("""
    CREATE TABLE prs_weights (
        trait TEXT NOT NULL,
        rsid TEXT NOT NULL,
        risk_allele TEXT NOT NULL,
        effect_weight REAL NOT NULL,
        PRIMARY KEY (trait, rsid)
    )""")

    # 4. PK/PD Mechanisms
    cursor.execute("""
    CREATE TABLE pk_pd_mechanisms (
        drug_name TEXT,
        gene_symbol TEXT,
        mechanism_type TEXT,
        effect_summary TEXT,
        PRIMARY KEY (drug_name, gene_symbol)
    )""")

    # 5. CPIC Guidelines
    cursor.execute("""
    CREATE TABLE cpic_guidelines (
        drug_name TEXT PRIMARY KEY,
        gene_symbol TEXT,
        status TEXT,
        recommendation TEXT
    )""")

    # 6. Drug-Drug Interactions
    cursor.execute("""
    CREATE TABLE ddi_rules (
        drug1 TEXT,
        drug2 TEXT,
        interaction_level TEXT,
        clinical_effect TEXT,
        PRIMARY KEY (drug1, drug2)
    )""")

    # 7. Disease Targeted Therapies
    cursor.execute("""
    CREATE TABLE disease_targeted_therapies (
        condition_trait TEXT PRIMARY KEY,
        primary_drug TEXT,
        gene_checked TEXT,
        pgx_status TEXT,
        alternative_drug TEXT,
        clinical_rationale TEXT
    )""")

    # --- SEED DATA INSERTION ---

    # Star Allele Diplotype Patterns
    cursor.executemany("INSERT INTO pgx_star_alleles VALUES (?, ?, ?, ?)", [
        ("CYP2D6", "*10/*10", "Poor Metabolizer", "rs1065852:T:T"),
        ("CYP2D6", "*4/*4", "Poor Metabolizer", "rs3892097:A:A"),
        ("CYP2D6", "*1/*1", "Normal Metabolizer", "DEFAULT"),
        ("CYP2C19", "*2/*2", "Poor Metabolizer", "rs4244285:A:A"),
        ("CYP2C19", "*2/*17", "Intermediate Metabolizer", "rs4244285:G:A"),
        ("CYP2C19", "*1/*1", "Normal Metabolizer", "DEFAULT"),
        ("CYP2C9", "*3/*3", "Poor Metabolizer", "rs1057910:C:C"),
        ("CYP2C9", "*1/*3", "Intermediate Metabolizer", "rs1057910:A:C"),
        ("CYP2C9", "*1/*1", "Normal Metabolizer", "DEFAULT"),
        ("SLCO1B1", "*5/*5", "Poor Function", "rs4149056:C:C"),
        ("SLCO1B1", "*1/*1", "Normal Function", "DEFAULT"),
        ("CYP3A5", "*1/*1", "Expresser", "DEFAULT")
    ])

    # PRS Trait Normalization Parameters (Population Baseline)
    cursor.executemany("INSERT INTO prs_traits VALUES (?, ?, ?)", [
        ("Schizoaffective Disorder", 1.20, 0.45),
        ("Bipolar Disorder", 1.10, 0.40),
        ("Major Depressive Disorder", 0.50, 0.25),
        ("Coronary Artery Disease", 2.10, 0.80),
        ("Type 2 Diabetes", 1.80, 0.60)
    ])

    # PRS Variants & Effect Weights
    cursor.executemany("INSERT INTO prs_weights VALUES (?, ?, ?, ?)", [
        ("Schizoaffective Disorder", "rs9272219", "T", 0.42),
        ("Schizoaffective Disorder", "rs1024611", "A", 0.38),
        ("Schizoaffective Disorder", "rs1065852", "T", 0.25),
        ("Bipolar Disorder", "rs1800497", "T", 0.51),
        ("Bipolar Disorder", "rs1024611", "A", 0.35),
        ("Bipolar Disorder", "rs4244285", "A", 0.22),
        ("Major Depressive Disorder", "rs1057910", "C", 0.31),
        ("Major Depressive Disorder", "rs1800497", "T", 0.28)
    ])

    # PK / PD Mechanisms
    cursor.executemany("INSERT INTO pk_pd_mechanisms VALUES (?, ?, ?, ?)", [
        ("Warfarin", "VKORC1", "PD", "Increased Sensitivity: Lower dosage required to achieve therapeutic INR target."),
        ("Warfarin", "CYP2C9", "PK", "Decreased Metabolism: Extended drug half-life and elevated systemic exposure."),
        ("Clopidogrel", "CYP2C19", "PK", "Impaired Conversion: Prodrug cannot be effectively converted to active thiol metabolite."),
        ("Simvastatin", "SLCO1B1", "PK", "Transporter Deficiency: Decreased hepatic clearance, increasing risk of statin-induced myopathy."),
        ("Aripiprazole", "CYP2D6", "PK", "Reduced Elimination: Drug accumulation increases sedation and extrapyramidal risk."),
        ("Carbamazepine", "HLA-B*15:02", "PD", "Immune Cytotoxicity: Direct activation of cytotoxic T-cells causing cutaneous necrosis.")
    ])

    # CPIC Guidelines
    cursor.executemany("INSERT INTO cpic_guidelines VALUES (?, ?, ?, ?)", [
        ("Aspirin", "CYP2C19", "SUITABLE", "Standard antiplatelet therapy."),
        ("Clopidogrel", "CYP2C19", "CONTRAINDICATED", "Avoid clopidogrel due to significantly reduced active metabolite formation. Switch to prasugrel or ticagrelor."),
        ("Simvastatin", "SLCO1B1", "SUITABLE", "Limit simvastatin dose to 20mg daily or switch to rosuvastatin/pravastatin."),
        ("Warfarin", "CYP2C9", "HIGH_RISK", "Reduce initial dose by 50-80% due to severely reduced clearance."),
        ("Aripiprazole", "CYP2D6", "HIGH_RISK", "Reduce initial dose by 50% due to impaired clearance and elevated plasma levels."),
        ("Risperidone", "CYP2D6", "HIGH_RISK", "Titrate slowly or reduce dose by 50%; monitor for extrapyramidal symptoms."),
        ("Clozapine", "CYP1A2", "SUITABLE", "Monitor trough serum concentrations; lower maintenance doses required."),
        ("Carbamazepine", "HLA-B*15:02", "CONTRAINDICATED", "Avoid due to high risk of Stevens-Johnson syndrome (SJS) and toxic epidermal necrolysis (TEN). Switch to Valproate or Lamotrigine."),
        ("Escitalopram", "CYP2C19", "CONTRAINDICATED", "Reduce starting dose by 50% or select alternative drug not predominant on CYP2C19."),
        ("Sertraline", "CYP2C19", "SUITABLE", "Consider 50% dose reduction if co-administered with CYP2D6 inhibitors."),
        ("Tacrolimus", "CYP3A5", "SUITABLE", "Standard starting dose required for non-expressers."),
        ("Fluorouracil", "DPYD", "CONTRAINDICATED", "Avoid use due to severe, potentially fatal toxicity.")
    ])

    # Drug-Drug Interactions
    cursor.executemany("INSERT INTO ddi_rules VALUES (?, ?, ?, ?)", [
        ("Clopidogrel", "Omeprazole", "Contraindicated", "Omeprazole inhibits CYP2C19, preventing Clopidogrel activation. Use Pantoprazole instead."),
        ("Warfarin", "Amiodarone", "Major", "Amiodarone significantly increases Warfarin concentrations. Reduce Warfarin dose by 30-50%."),
        ("Aripiprazole", "Fluoxetine", "Major", "Fluoxetine doubles Aripiprazole exposure. Reduce Aripiprazole dose by 50%.")
    ])

    # Disease Targeted Therapies
    cursor.executemany("INSERT INTO disease_targeted_therapies VALUES (?, ?, ?, ?, ?, ?)", [
        ("Schizoaffective Disorder", "Aripiprazole / Risperidone", "CYP2D6", "HIGH_RISK", "Clozapine / Olanzapine", "First-line atypical antipsychotic for elevated schizoaffective polygenic risk. Check CYP2D6 metabolizer status."),
        ("Coronary Artery Disease", "Simvastatin", "SLCO1B1", "SUITABLE", "Rosuvastatin", "Primary prevention for elevated CAD burden. Adjust dose or switch to Rosuvastatin if SLCO1B1 impaired."),
        ("Factor V Leiden Thrombophilia", "Warfarin", "CYP2C9", "HIGH_RISK", "Direct Oral Anticoagulant (DOAC)", "Anticoagulation indicated for thrombotic risk. Reduce Warfarin or use DOAC if CYP2C9 poor metabolizer.")
    ])

    conn.commit()
    conn.close()
    print("[✔] Database successfully upgraded to Enterprise Schema.")

if __name__ == "__main__":
    setup_database()
