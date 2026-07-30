#!/usr/bin/env python3
import sqlite3
from pathlib import Path

def initialize_database():
    db_path = Path("genomic_knowledgebase.db")
    print(f"[+] Initializing genomic knowledgebase at: {db_path.resolve()}")
    
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()

    # Drop old tables if they exist to ensure schema alignment
    cursor.execute("DROP TABLE IF EXISTS knowledgebase;")
    cursor.execute("DROP TABLE IF EXISTS prs_therapeutic_guidelines;")

    # 1. Create the primary knowledgebase table expected by bot.py
    cursor.execute("""
    CREATE TABLE knowledgebase (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        drug_name TEXT NOT NULL,
        therapeutic_class TEXT NOT NULL,
        target_disorder TEXT NOT NULL,
        gene_symbol TEXT NOT NULL,
        evidence_tier TEXT NOT NULL,
        recommendation TEXT NOT NULL
    );
    """)

    # 2. Create the new PRS-to-Therapeutics cross-reference table
    cursor.execute("""
    CREATE TABLE prs_therapeutic_guidelines (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        prs_id TEXT NOT NULL,
        condition_name TEXT NOT NULL,
        recommended_intervention_class TEXT NOT NULL,
        clinical_rationale TEXT NOT NULL
    );
    """)

    # Comprehensive Western medicine records matching bot.py query parameters
    sample_data = [
        ("Escitalopram (Lexapro)", "SSRI Antidepressant", "Major Depressive Disorder / Anxiety", "CYP2C19", "1A", "Monitor plasma levels due to CYP2C19 pathway sensitivity; consider lower initial dose titration."),
        ("Sertraline (Zoloft)", "SSRI Antidepressant", "Major Depressive Disorder", "CYP2C19", "1A", "Processed normally via CYP2C19. Standard starting dose recommended."),
        ("Clopidogrel (Plavix)", "Antiplatelet", "Cardiovascular Disease / Thrombosis", "CYP2C19", "1A", "Poor metabolizers exhibit significantly reduced antiplatelet effect. Consider alternative therapy like prasugrel or ticagrelor."),
        ("Metoprolol Succinate", "Beta Blocker", "Hypertension / Heart Failure", "CYP2D6", "1A", "Normal clearance profile; standard titration guidelines apply for hypertension/heart failure."),
        ("Tramadol", "Analgesic", "Pain Management", "CYP2D6", "1A", "Ultrarapid metabolizers risk opioid toxicity; poor metabolizers lack analgesic efficacy. Verify conversion efficiency."),
        ("Codeine", "Analgesic", "Pain Management", "CYP2D6", "1A", "Avoid use in ultra-rapid (respiratory depression risk) and poor metabolizers (lack of pain relief)."),
        ("Simvastatin", "Statin", "Hypercholesterolemia", "SLCO1B1", "1A", "Increased risk of myopathy with variants. Consider lower dose or alternative statin like rosuvastatin."),
        ("Carbamazepine", "Anticonvulsant", "Epilepsy / Bipolar Disorder", "HLA-B", "1A", "Strong association with Stevens-Johnson syndrome/toxic epidermal necrolysis in HLA-B*15:02 carriers. Avoid if positive.")
    ]

    cursor.executemany("""
    INSERT INTO knowledgebase (drug_name, therapeutic_class, target_disorder, gene_symbol, evidence_tier, recommendation)
    VALUES (?, ?, ?, ?, ?, ?)
    """, sample_data)

    # Sample PRS-to-Therapeutics mapping data
    sample_prs_guidelines = [
        ('PGS000018', 'Coronary Artery Disease', 'High-Potency Statins / Prevention', 'High polygenic cardiovascular load warrants aggressive lipid-lowering therapy, subject to SLCO1B1 clearance.'),
        ('PGS000034', 'Major Depressive Disorder', 'First-Line SSRIs (CYP2C19 guided)', 'Elevated depression risk combined with metabolic profile guides initial antidepressant selection.'),
        ('PGS000014', 'Type 2 Diabetes (T2D)', 'Metabolic & Glycemic Monitoring', 'Elevated polygenic diabetes risk warrants proactive glycemic tracking and lifestyle/pharmacological intervention.'),
        ('PGS000021', 'Atrial Fibrillation', 'Cardiovascular Monitoring', 'Low polygenic predisposition detected for rhythm anomalies.')
    ]

    cursor.executemany("""
    INSERT INTO prs_therapeutic_guidelines (prs_id, condition_name, recommended_intervention_class, clinical_rationale)
    VALUES (?, ?, ?, ?)
    """, sample_prs_guidelines)

    conn.commit()
    conn.close()
    print("[✔] 'genomic_knowledgebase.db' successfully generated and populated with matching schema and PRS guidelines.")

if __name__ == "__main__":
    initialize_database()
