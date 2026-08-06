#!/usr/bin/env python3
"""
init_db.py - Unified Precision Medicine Knowledgebase Initializer
Combines Western Drugs, RxCUI/ATC Mappings, DDI Rules, CPIC DGIs, 
ClinVar Classifications, and PRS-to-Therapeutic Guidelines.
"""

import sqlite3
from pathlib import Path

DB_PATH = Path("genomic_knowledgebase.db")

def initialize_database():
    print(f"[+] Initializing genomic knowledgebase at: {DB_PATH.resolve()}")
    
    conn = sqlite3.connect(str(DB_PATH))
    cursor = conn.cursor()

    # Enable foreign keys
    cursor.execute("PRAGMA foreign_keys = ON;")

    # Drop old tables to guarantee clean alignment
    tables_to_drop = [
        "knowledgebase", "prs_therapeutic_guidelines", 
        "ddi_pair_rules", "ddi_class_rules", "dgi_rules", "clinvar_variants"
    ]
    for table in tables_to_drop:
        cursor.execute(f"DROP TABLE IF EXISTS {table};")

    # 1. CORE WESTERN MEDICINE KNOWLEDGEBASE
    cursor.execute("""
    CREATE TABLE knowledgebase (
        rxcui TEXT PRIMARY KEY,
        drug_name TEXT NOT NULL,
        therapeutic_class TEXT NOT NULL,
        atc_code TEXT,
        target_disorder TEXT NOT NULL,
        gene_symbol TEXT NOT NULL,
        evidence_tier TEXT NOT NULL,
        recommendation TEXT NOT NULL
    );
    """)

    # 2. PRS-TO-THERAPEUTICS CROSS-REFERENCE
    cursor.execute("""
    CREATE TABLE prs_therapeutic_guidelines (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        prs_id TEXT NOT NULL,
        condition_name TEXT NOT NULL,
        recommended_intervention_class TEXT NOT NULL,
        clinical_rationale TEXT NOT NULL
    );
    """)

    # 3. PAIRWISE DDI RULES (rxcui_a + rxcui_b)
    cursor.execute("""
    CREATE TABLE ddi_pair_rules (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        rxcui_a TEXT NOT NULL,
        rxcui_b TEXT NOT NULL,
        severity TEXT CHECK(severity IN ('Minor', 'Moderate', 'Major', 'Contraindicated')),
        mechanism TEXT,
        clinical_effect TEXT,
        FOREIGN KEY (rxcui_a) REFERENCES knowledgebase(rxcui),
        FOREIGN KEY (rxcui_b) REFERENCES knowledgebase(rxcui),
        UNIQUE(rxcui_a, rxcui_b)
    );
    """)

    # 4. CLASS-BASED DDI RULES (via ATC Classes)
    cursor.execute("""
    CREATE TABLE ddi_class_rules (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        class_a_code TEXT NOT NULL,
        class_b_code TEXT NOT NULL,
        severity TEXT CHECK(severity IN ('Minor', 'Moderate', 'Major', 'Contraindicated')),
        clinical_effect TEXT,
        UNIQUE(class_a_code, class_b_code)
    );
    """)

    # 5. PHARMGKB / CPIC DRUG-GENE INTERACTION (DGI) RULES
    cursor.execute("""
    CREATE TABLE dgi_rules (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        rxcui TEXT NOT NULL,
        gene_symbol TEXT NOT NULL,
        phenotype TEXT NOT NULL,
        recommendation TEXT NOT NULL,
        cpic_level TEXT,
        FOREIGN KEY (rxcui) REFERENCES knowledgebase(rxcui)
    );
    """)

    # 6. CLINVAR VARIANT CLASSIFICATIONS
    cursor.execute("""
    CREATE TABLE clinvar_variants (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        rsid TEXT,
        chrom TEXT NOT NULL,
        pos INTEGER NOT NULL,
        ref TEXT NOT NULL,
        alt TEXT NOT NULL,
        gene_symbol TEXT,
        clinical_significance TEXT,
        associated_trait TEXT,
        review_status TEXT,
        UNIQUE(chrom, pos, ref, alt)
    );
    """)

    # Performance Indexing
    cursor.execute("CREATE INDEX idx_ddi_pair ON ddi_pair_rules(rxcui_a, rxcui_b);")
    cursor.execute("CREATE INDEX idx_ddi_class ON ddi_class_rules(class_a_code, class_b_code);")
    cursor.execute("CREATE INDEX idx_dgi_lookup ON dgi_rules(rxcui, gene_symbol);")
    cursor.execute("CREATE INDEX idx_clinvar_pos ON clinvar_variants(chrom, pos, ref, alt);")

    # --- POPULATE SEED DATA ---

    # Core Western Meds with RxCUIs and ATC Codes
    western_meds = [
        ('36437', 'Escitalopram (Lexapro)', 'SSRI Antidepressant', 'N06AB10', 'Major Depressive Disorder / Anxiety', 'CYP2C19', '1A', 'Monitor plasma levels due to CYP2C19 pathway sensitivity; consider lower initial dose titration.'),
        ('99280', 'Sertraline (Zoloft)', 'SSRI Antidepressant', 'N06AB06', 'Major Depressive Disorder', 'CYP2C19', '1A', 'Processed normally via CYP2C19. Standard starting dose recommended.'),
        ('214159', 'Clopidogrel (Plavix)', 'Antiplatelet', 'B01AC04', 'Cardiovascular Disease / Thrombosis', 'CYP2C19', '1A', 'Poor metabolizers exhibit significantly reduced antiplatelet effect. Consider alternative therapy like prasugrel or ticagrelor.'),
        ('6918', 'Metoprolol Succinate', 'Beta Blocker', 'C07AB02', 'Hypertension / Heart Failure', 'CYP2D6', '1A', 'Normal clearance profile; standard titration guidelines apply for hypertension/heart failure.'),
        ('10689', 'Tramadol', 'Analgesic', 'N02AJ13', 'Pain Management', 'CYP2D6', '1A', 'Ultrarapid metabolizers risk opioid toxicity; poor metabolizers lack analgesic efficacy. Verify conversion efficiency.'),
        ('2670', 'Codeine', 'Analgesic', 'R05DA04', 'Pain Management', 'CYP2D6', '1A', 'Avoid use in ultra-rapid (respiratory depression risk) and poor metabolizers (lack of pain relief).'),
        ('36567', 'Simvastatin', 'Statin', 'C10AA01', 'Hypercholesterolemia', 'SLCO1B1', '1A', 'Increased risk of myopathy with variants. Consider lower dose or alternative statin like rosuvastatin.'),
        ('2034', 'Carbamazepine', 'Anticonvulsant', 'N03AF01', 'Epilepsy / Bipolar Disorder', 'HLA-B', '1A', 'Strong association with Stevens-Johnson syndrome/toxic epidermal necrolysis in HLA-B*15:02 carriers. Avoid if positive.'),
        ('1191', 'Aspirin', 'Antiplatelet', 'B01AC06', 'Coronary Artery Disease', 'PTGS1', '1A', 'Standard dosage protocol.'),
        ('11289', 'Warfarin', 'Anticoagulant', 'B01AA03', 'Thromboembolism', 'VKORC1', '1A', 'Adjust initial dose based on VKORC1 / CYP2C9 genotypes.')
    ]
    cursor.executemany("""
    INSERT INTO knowledgebase (rxcui, drug_name, therapeutic_class, atc_code, target_disorder, gene_symbol, evidence_tier, recommendation)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?)
    """, western_meds)

    # PRS Guidelines
    prs_guidelines = [
        ('PGS000018', 'Coronary Artery Disease', 'High-Potency Statins / Prevention', 'High polygenic cardiovascular load warrants aggressive lipid-lowering therapy, subject to SLCO1B1 clearance.'),
        ('PGS000034', 'Major Depressive Disorder', 'First-Line SSRIs (CYP2C19 guided)', 'Elevated depression risk combined with metabolic profile guides initial antidepressant selection.'),
        ('PGS000014', 'Type 2 Diabetes (T2D)', 'Metabolic & Glycemic Monitoring', 'Elevated polygenic diabetes risk warrants proactive glycemic tracking and lifestyle/pharmacological intervention.'),
        ('PGS000021', 'Atrial Fibrillation', 'Cardiovascular Monitoring', 'Low polygenic predisposition detected for rhythm anomalies.')
    ]
    cursor.executemany("""
    INSERT INTO prs_therapeutic_guidelines (prs_id, condition_name, recommended_intervention_class, clinical_rationale)
    VALUES (?, ?, ?, ?)
    """, prs_guidelines)

    # DDI Rules (Pairwise)
    ddi_pairs = [
        ('1191', '11289', 'Major', 'Pharmacodynamic Synergism', 'Increased risk of serious gastrointestinal and systemic bleeding.'),
        ('214159', '1191', 'Moderate', 'Additive Antiplatelet Effect', 'Increased bleeding risk; monitor dual antiplatelet therapy closely.')
    ]
    cursor.executemany("""
    INSERT INTO ddi_pair_rules (rxcui_a, rxcui_b, severity, mechanism, clinical_effect)
    VALUES (?, ?, ?, ?, ?)
    """, ddi_pairs)

    # DDI Rules (Class-Level)
    ddi_classes = [
        ('N06AB', 'B01AA', 'Moderate', 'SSRIs combined with oral anticoagulants increase mucosal bleeding risk.'),
        ('B01AC', 'B01AA', 'Major', 'Combined antiplatelet and anticoagulant therapy significantly elevates major hemorrhage risks.')
    ]
    cursor.executemany("""
    INSERT INTO ddi_class_rules (class_a_code, class_b_code, severity, clinical_effect)
    VALUES (?, ?, ?, ?)
    """, ddi_classes)

    # CPIC DGI Rules
    dgi_rules = [
        ('214159', 'CYP2C19', 'Poor Metabolizer', 'Switch to prasugrel or ticagrelor due to loss of efficacy.', 'A'),
        ('36437', 'CYP2C19', 'Poor Metabolizer', 'Reduce starting dose by 50% or select alternative antidepressant.', 'A'),
        ('36567', 'SLCO1B1', 'Decreased Function', 'Limit dose to 20mg daily or switch to rosuvastatin to lower myopathy risk.', 'A')
    ]
    cursor.executemany("""
    INSERT INTO dgi_rules (rxcui, gene_symbol, phenotype, recommendation, cpic_level)
    VALUES (?, ?, ?, ?, ?)
    """, dgi_rules)

    # ClinVar Seed Variants
    clinvar_seed = [
        ('rs429358', '19', 44908684, 'T', 'C', 'APOE', 'Pathogenic', 'Alzheimer Disease 2', 'practice guideline'),
        ('rs1057910', '10', 94942290, 'A', 'C', 'CYP2C9', 'Drug response', 'Warfarin response', 'reviewed by expert panel')
    ]
    cursor.executemany("""
    INSERT INTO clinvar_variants (rsid, chrom, pos, ref, alt, gene_symbol, clinical_significance, associated_trait, review_status)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
    """, clinvar_seed)

    conn.commit()
    conn.close()
    print("[✔] Knowledgebase schema initialized and populated successfully.")

if __name__ == "__main__":
    initialize_database()