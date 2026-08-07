#!/usr/bin/env python3
"""
init_db.py - Initializes and populates local genomic, ClinVar, PRS weights, and PGx knowledgebase.
"""

import sqlite3
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
DB_PATH = Path("genomic_knowledgebase.db")


def init_database(db_path: Path = DB_PATH) -> None:
    logging.info(f"Initializing comprehensive genomic database at: {db_path}")
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # 1. Drug Knowledgebase Table
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

    # 2. ClinVar Variants Table (Monogenic / Pathogenic & ACMG Secondary Findings)
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

    # 3. Polygenic Risk Score (PRS) Weights Table
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS prs_weights (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        trait_name TEXT NOT NULL,
        rsid TEXT,
        chrom TEXT,
        pos INTEGER,
        effect_allele TEXT,
        other_allele TEXT,
        weight REAL NOT NULL
    );
    """)

    # 4. CPIC DGI Rules Table (PGx Phenotyping Matrix)
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

    # 5. DDI Pair Rules Table
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

    # --- SEED KNOWLEDGEBASE ---
    kb_data = [
        ("1191", "Aspirin", "Antiplatelet", "Cardiovascular Disease", "B01AC06", "B01AC", "CYP2C19", "Standard dosing."),
        ("32968", "Clopidogrel", "Antiplatelet", "Thrombosis", "B01AC04", "B01AC", "CYP2C19", "Alternative antiplatelet therapy recommended if poor metabolizer."),
        ("36437", "Escitalopram", "SSRI", "Depression", "N06AB10", "N06AB", "CYP2C19", "Consider 50% reduction in dose for poor metabolizers."),
        ("36567", "Simvastatin", "HMG-CoA Reductase Inhibitor", "Hyperlipidemia", "C10AA01", "C10AA", "SLCO1B1", "Prescribe lower dose or alternative statin for decreased function."),
        ("11289", "Warfarin", "Anticoagulant", "Thromboembolism", "B01AA03", "B01AA", "CYP2C9", "Titrate dose carefully based on VKORC1 and CYP2C9 genotype.")
    ]
    cursor.executemany("INSERT OR REPLACE INTO knowledgebase VALUES (?, ?, ?, ?, ?, ?, ?, ?);", kb_data)

    # --- SEED CLINVAR DATA ---
    clinvar_data = [
        ("GRCh38", "19", 44908684, "T", "C", "rs429358", "APOE", "Pathogenic", "Alzheimer Disease", "practice guideline"),
        ("GRCh38", "10", 94942290, "A", "C", "rs1057910", "CYP2C9", "Pathogenic/Likely pathogenic", "Warfarin response", "criteria provided"),
        ("GRCh38", "13", 32315474, "C", "T", "rs80357906", "BRCA2", "Pathogenic", "Hereditary Breast and Ovarian Cancer", "reviewed by expert panel"),
        ("GRCh38", "1", 169519049, "A", "G", "rs1799853", "F5", "Pathogenic", "Factor V Leiden Thrombophilia", "criteria provided")
    ]
    cursor.executemany("""
    INSERT OR REPLACE INTO clinvar_variants 
    (genome_build, chrom, pos, ref, alt, rsid, gene_symbol, clinical_significance, associated_trait, review_status)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?);
    """, clinvar_data)

    # --- SEED PRS WEIGHTS ---
    prs_data = [
        ("Coronary Artery Disease", "rs1333049", "9", 22125503, "C", "T", 0.45),
        ("Coronary Artery Disease", "rs6922269", "6", 16091321, "A", "G", 0.31),
        ("Type 2 Diabetes", "rs7903146", "10", 114754029, "T", "C", 0.58),
        ("Type 2 Diabetes", "rs12255372", "10", 114757187, "T", "G", 0.42),
        ("Major Depressive Disorder", "rs2522122", "3", 51234500, "A", "G", 0.25),
        ("Alzheimer's Disease", "rs429358", "19", 44908684, "C", "T", 1.20)
    ]
    cursor.executemany("""
    INSERT OR REPLACE INTO prs_weights (trait_name, rsid, chrom, pos, effect_allele, other_allele, weight)
    VALUES (?, ?, ?, ?, ?, ?, ?);
    """, prs_data)

    # --- SEED DGI RULES ---
    dgi_data = [
        ("32968", "CYP2C19", "Poor Metabolizer", "A", "Avoid clopidogrel due to significantly reduced active metabolite formation. Switch to prasugrel or ticagrelor."),
        ("36437", "CYP2C19", "Poor Metabolizer", "A", "Reduce starting dose by 50% or select alternative drug not predominant on CYP2C19."),
        ("36567", "SLCO1B1", "Decreased Function", "A", "Limit simvastatin dose to 20mg daily or switch to rosuvastatin/pravastatin."),
        ("11289", "CYP2C9", "Poor Metabolizer", "A", "Initiate warfarin at significantly reduced maintenance dose (e.g., 1-3 mg/day) per pharmacogenetic dosing algorithm.")
    ]
    cursor.executemany("""
    INSERT OR REPLACE INTO dgi_rules (rxcui, gene_symbol, phenotype, cpic_level, recommendation)
    VALUES (?, ?, ?, ?, ?);
    """, dgi_data)

    # --- SEED DDI RULES ---
    ddi_pair_data = [
        ("1191", "32968", "Moderate", "Pharmacodynamic synergy", "Increased risk of bleeding when Aspirin is combined with Clopidogrel."),
        ("11289", "1191", "Major", "Additive antiplatelet/anticoagulant effects", "Severe bleeding risk elevation when Warfarin is combined with Aspirin.")
    ]
    cursor.executemany("""
    INSERT OR REPLACE INTO v_ddi_pair_rules (rxcui_a, rxcui_b, severity, mechanism, clinical_effect)
    VALUES (?, ?, ?, ?, ?);
    """, ddi_pair_data)

    conn.commit()
    conn.close()
    logging.info("Database initialized and populated successfully.")


if __name__ == "__main__":
    init_database()
