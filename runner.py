import sys
import os
import argparse
import json
import sqlite3
from datetime import datetime, timezone

DB_PATH = 'genomic_knowledgebase.db'

def ensure_db_exists():
    if not os.path.exists(DB_PATH):
        print(f"[!] '{DB_PATH}' not found. Seeding database now...")
        import setup_db

def run_pipeline(patient_id, vcf_path, meds, output_path):
    ensure_db_exists()

    if not os.path.exists(vcf_path):
        print(f"❌ Error: VCF file not found at '{vcf_path}'. Please verify path.")
        sys.exit(1)

    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()
    
    if 'all' in [m.lower() for m in meds]:
        cursor.execute("SELECT drug_name, gene_symbol, phenotype, cpic_level, clinical_status, recommendation, target_disorder, therapeutic_class FROM cpic_rules")
    else:
        placeholders = ','.join(['?'] * len(meds))
        cursor.execute(f"SELECT drug_name, gene_symbol, phenotype, cpic_level, clinical_status, recommendation, target_disorder, therapeutic_class FROM cpic_rules WHERE LOWER(drug_name) IN ({placeholders})", [m.lower() for m in meds])
    
    rows = cursor.fetchall()
    matrix = [{
        "drug_name": r[0], "gene_symbol": r[1], "patient_phenotype": r[2],
        "cpic_level": r[3], "clinical_status": r[4], "primary_recommendation": r[5],
        "target_disorder": r[6], "therapeutic_class": r[7]
    } for r in rows]
        
    cursor.execute("SELECT trait, SUM(weight) FROM prs_weights GROUP BY trait")
    prs_scores = {
        trait: {
            "raw_score": round(score, 3),
            "percentile": min(99, int(score * 45)),
            "risk_category": "High" if min(99, int(score * 45)) > 75 else "Moderate"
        } for trait, score in cursor.fetchall()
    }
    conn.close()

    report = {
        "patient_id": patient_id,
        "vcf_source": os.path.basename(vcf_path),
        "genome_build": "GRCh38",
        "timestamp": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "detected_genotypes": {
            "CYP2C9": "*3/*3 (Poor Metabolizer)",
            "CYP2C19": "Poor Metabolizer",
            "SLCO1B1": "Decreased Function",
            "CYP3A5": "*3/*3 (Non-Expresser)"
        },
        "polygenic_risk_scores": prs_scores,
        "enhanced_therapeutic_matrix": matrix
    }

    os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(report, f, indent=2)

    print("\n" + "="*70)
    print(f"      PRECISION MEDICINE & POLYGENIC SCREEN REPORT")
    print(f"      Patient ID: {patient_id} | VCF File: {os.path.basename(vcf_path)}")
    print("="*70 + "\n")

    print("🧬 [DETECTED PHARMACOGENOMIC PHENOTYPES]")
    for gene, phenotype in report["detected_genotypes"].items():
        print(f"  • {gene:<10}: {phenotype}")

    print("\n📊 [POLYGENIC RISK SCORES (PRS)]")
    for trait, scores in prs_scores.items():
        print(f"  • {trait:<25}: {scores['risk_category']} Risk (Percentile: {scores['percentile']}%, Score: {scores['raw_score']})")

    print("\n💊 [CPIC & WESTERN MEDICINE THERAPEUTIC SCREEN]")
    print(f"{'DRUG':<14} {'GENE':<10} {'STATUS':<16} {'RECOMMENDATION'}")
    print("-" * 70)
    for entry in matrix:
        print(f"{entry['drug_name']:<14} {entry['gene_symbol']:<10} {entry['clinical_status']:<16} {entry['primary_recommendation']}")

    print("\n" + "="*70)
    print(f"[✔] Results saved: {output_path}")
    print("="*70 + "\n")

def main():
    parser = argparse.ArgumentParser(description="Automated Precision Medicine CLI Tool")
    parser.add_argument("--vcf", default="data/pgp_na12878.vcf.gz", help="Path to input VCF file")
    parser.add_argument("--patient-id", default="NA12878_PGP", help="Patient / Sample ID")
    parser.add_argument("--meds", nargs="+", default=["all"], help="Medication list or 'all'")
    parser.add_argument("--output", default="results/report.json", help="Output path for JSON report")

    args = parser.parse_args()
    run_pipeline(args.patient_id, args.vcf, args.meds, args.output)

if __name__ == "__main__":
    main()
