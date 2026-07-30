#!/usr/bin/env python3
import sqlite3
from pathlib import Path

def initialize_database():
    db_path = Path("genomic_knowledgebase.db")
    print(f"[+] Initializing genomic knowledgebase at: {db_path.resolve()}")
    
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()

    # Create the schema expected by bot.py
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS clinical_guidelines (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        gene TEXT NOT NULL,
        drug_name TEXT NOT NULL,
        guideline_source TEXT NOT NULL,
        recommendation TEXT NOT NULL,
        evidence_level TEXT NOT NULL
    );
    """)

    # High-accuracy western medicine records (CPIC / PharmGKB guidelines)
    sample_data = [
        ("CYP2C19", "Escitalopram (Lexapro)", "CPIC Guideline", "Monitor plasma levels due to CYP2C19 pathway sensitivity; consider lower initial dose titration.", "1A"),
        ("CYP2C19", "Sertraline (Zoloft)", "CPIC Guideline", "Processed normally via CYP2C19. Standard starting dose recommended.", "1A"),
        ("CYP2C19", "Clopidogrel (Plavix)", "CPIC Guideline", "Poor metabolizers exhibit significantly reduced antiplatelet effect. Consider alternative therapy like prasugrel or ticagrelor.", "1A"),
        ("CYP2D6", "Metoprolol Succinate", "CPIC Guideline", "Normal clearance profile; standard titration guidelines apply for hypertension/heart failure.", "1A"),
        ("CYP2D6", "Tramadol", "CPIC Guideline", "Ultrarapid metabolizers risk opioid toxicity; poor metabolizers lack analgesic efficacy. Verify conversion efficiency.", "1A"),
        ("CYP2D6", "Codeine", "CPIC Guideline", "Avoid use in ultra-rapid (respiratory depression risk) and poor metabolizers (lack of pain relief).", "1A"),
        ("SLCO1B1", "Simvastatin", "CPIC Guideline", "Increased risk of myopathy with variants. Consider lower dose or alternative statin like rosuvastatin.", "1A"),
        ("HLA-B", "Carbamazepine", "CPIC Guideline", "Strong association with Stevens-Johnson syndrome/toxic epidermal necrolysis in HLA-B*15:02 carriers. Avoid if positive.", "1A")
    ]

    cursor.executemany("""
    INSERT INTO clinical_guidelines (gene, drug_name, guideline_source, recommendation, evidence_level)
    VALUES (?, ?, ?, ?, ?)
    """, sample_data)

    conn.commit()
    conn.close()
    print("[✔] 'genomic_knowledgebase.db' successfully generated and populated.")

if __name__ == "__main__":
    initialize_database()