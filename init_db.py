#!/usr/bin/env python3
import sqlite3
from pathlib import Path

def initialize_database():
    db_path = Path("genomic_knowledgebase.db")
    print(f"[+] Initializing genomic knowledgebase at: {db_path.resolve()}")
    
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()

    # Drop old table if it exists to ensure schema alignment with bot.py
    cursor.execute("DROP TABLE IF EXISTS knowledgebase;")

    # Create the exact schema expected by bot.py
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

    conn.commit()
    conn.close()
    print("[✔] 'genomic_knowledgebase.db' successfully generated and populated with matching schema.")

if __name__ == "__main__":
    initialize_database()
