#!/usr/bin/env python3
"""
init_db.py - Initializes and populates local genomic and drug knowledgebase.
"""

import sqlite3
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
DB_PATH = Path("genomic_knowledgebase.db")


def init_database(db_path: Path = DB_PATH) -> None:
    logging.info(f"Initializing genomic database at: {db_path}")
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # 1. Knowledgebase Table
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS knowledgebase (
        rxcui TEXT PRIMARY KEY,
        drug_name TEXT NOT NULL,
        therapeutic_class TEXT,
        target_disorder TEXT,
        atc_code TEXT,
        atc_5_prefix TEXT,
        gene_symbol TEXT,
        recommendation TEXT
    );
    """)

    # 2. ClinVar Variants Table
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS clinvar_variants (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        genome_build TEXT NOT NULL,
        chrom TEXT NOT NULL,
        pos INTEGER NOT NULL,
        ref TEXT NOT NULL,
        alt TEXT NOT NULL,
        rsid TEXT,
        gene_symbol TEXT,
        clinical_significance TEXT,
        associated_trait TEXT,
        review_status TEXT
    );
    """)

    # 3. DGI Rules Table
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS dgi_rules (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        rxcui TEXT NOT NULL,
        gene_symbol TEXT NOT NULL,
        phenotype TEXT NOT NULL,
        cpic_level TEXT,
        recommendation TEXT NOT NULL
    );
    """)

    # 4. DDI Pair Rules Table
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS v_ddi_pair_rules (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        rxcui_a TEXT NOT NULL,
        rxcui_b TEXT NOT NULL,
        severity TEXT NOT NULL,
        mechanism TEXT,
        clinical_effect TEXT NOT NULL
    );
    """)

    # 5. DDI Class Rules Table
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS ddi_class_rules (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        class_a_code TEXT NOT NULL,
        class_b_code TEXT NOT NULL,
        severity TEXT NOT NULL,
        clinical_effect TEXT NOT NULL
    );
    """)

    # Seed Knowledgebase Data
    kb_data = [
        ("1191", "Aspirin", "Antiplatelet", "Cardiovascular Disease", "B01AC06", "B01AC", "CYP2C19", "Standard dosing."),
        ("32968", "Clopidogrel", "Antiplatelet", "Thrombosis", "B01AC04", "B01AC", "CYP2C19", "Alternative antiplatelet therapy recommended if poor metabolizer."),
        ("36437", "Escitalopram", "SSRI", "Depression", "N06AB10", "N06AB", "CYP2C19", "Consider 50% reduction in dose for poor metabolizers."),
        ("36567", "Simvastatin", "HMG-CoA Reductase Inhibitor", "Hyperlipidemia", "C10AA01", "C10AA", "SLCO1B1", "Prescribe lower dose or alternative statin for decreased function.")
    ]
    cursor.executemany("""
    INSERT OR REPLACE INTO knowledgebase VALUES (?, ?, ?, ?, ?, ?, ?, ?);
    """, kb_data)

    # Seed ClinVar Data
    clinvar_data = [
        ("GRCh38", "19", 44908684, "T", "C", "rs429358", "APOE", "Pathogenic", "Alzheimer Disease", "practice guideline"),
        ("GRCh38", "10", 94942290, "A", "C", "rs1057910", "CYP2C9", "Pathogenic/Likely pathogenic", "Warfarin response", "criteria provided, multiple submitters")
    ]
    cursor.executemany("""
    INSERT OR REPLACE INTO clinvar_variants 
    (genome_build, chrom, pos, ref, alt, rsid, gene_symbol, clinical_significance, associated_trait, review_status)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?);
    """, clinvar_data)

    # Seed DGI Rules
    dgi_data = [
        ("32968", "CYP2C19", "Poor Metabolizer", "A", "Avoid clopidogrel due to significantly reduced active metabolite formation. Switch to prasugrel or ticagrelor."),
        ("36437", "CYP2C19", "Poor Metabolizer", "A", "Reduce starting dose by 50% or select alternative drug not predominant on CYP2C19."),
        ("36567", "SLCO1B1", "Decreased Function", "A", "Limit simvastatin dose to 20mg daily or switch to rosuvastatin/pravastatin.")
    ]
    cursor.executemany("""
    INSERT OR REPLACE INTO dgi_rules (rxcui, gene_symbol, phenotype, cpic_level, recommendation)
    VALUES (?, ?, ?, ?, ?);
    """, dgi_data)

    # Seed DDI Pairwise Rules
    ddi_pair_data = [
        ("1191", "32968", "Moderate", "Pharmacodynamic synergy", "Increased risk of bleeding when Aspirin is combined with Clopidogrel.")
    ]
    cursor.executemany("""
    INSERT OR REPLACE INTO v_ddi_pair_rules (rxcui_a, rxcui_b, severity, mechanism, clinical_effect)
    VALUES (?, ?, ?, ?, ?);
    """, ddi_pair_data)

    conn.commit()
    conn.close()
    logging.info("Database initialized successfully.")


if __name__ == "__main__":
    init_database()