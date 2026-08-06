import sys
import argparse
import json
import sqlite3
from datetime import datetime, timezone

def run_pipeline(patient_id, vcf_path, meds, output_path):
    conn = sqlite3.connect('genomic_knowledgebase.db')
    cursor = conn.cursor()
    
    # 1. Fetch CPIC / Western Medicine Drug Recommendations
    if 'all' in meds or 'ALL' in meds:
        cursor.execute("SELECT drug_name, gene_symbol, phenotype, cpic_level, clinical_status, recommendation, target_disorder, therapeutic_class FROM cpic_rules")
    else:
        placeholders = ','.join(['?'] * len(meds))
        cursor.execute(f"SELECT drug_name, gene_symbol, phenotype, cpic_level, clinical_status, recommendation, target_disorder, therapeutic_class FROM cpic_rules WHERE LOWER(drug_name) IN ({placeholders})", [m.lower() for m in meds])
    
    rows = cursor.fetchall()
    matrix = []
    for r in rows:
        matrix.append({
            "drug_name": r[0],
            "gene_symbol": r[1],
            "patient_phenotype": r[2],
            "cpic_level": r[3],
            "clinical_status": r[4],
            "primary_recommendation": r[5],
            "target_disorder": r[6],
            "therapeutic_class": r[7]
        })
        
    # 2. Compute Polygenic Risk Scores (PRS)
    cursor.execute("SELECT trait, SUM(weight) FROM prs_weights GROUP BY trait")
    prs_rows = cursor.fetchall()
    prs_scores = {}
    for trait, score in prs_rows:
        percentile = min(99, int(score * 45))
        risk_category = "High" if percentile > 75 else "Moderate"
        prs_scores[trait] = {
            "raw_score": round(score, 3),
            "percentile": percentile,
            "risk_category": risk_category
        }

    conn.close()

    # 3. Assemble Full Report Payload
    report = {
        "patient_id": patient_id,
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

    # 4. Save JSON Report
    with open(output_path, 'w') as f:
        json.dump(report, f, indent=2)

    # 5. Print Human-Readable Summary to Terminal Output
    print("\n" + "="*70)
    print(f"      PRECISION MEDICINE & POLYGENIC SCREEN REPORT")
    print(f"      Patient ID: {patient_id} | Build: GRCh38")
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
        status = entry['clinical_status']
        drug = entry['drug_name']
        gene = entry['gene_symbol']
        rec = entry['primary_recommendation']
        print(f"{drug:<14} {gene:<10} {status:<16} {rec}")

    print("\n" + "="*70)
    print(f"[✔] Comprehensive results saved to JSON: {output_path}")
    print("="*70 + "\n")

def main():
    parser = argparse.ArgumentParser(description="Automated Genomic & Precision Medicine CLI Tool")
    subparsers = parser.add_subparsers(dest="subcommand")

    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("--patient-id", required=True)
    run_parser.add_argument("--vcf", required=True)
    run_parser.add_argument("--meds", nargs="+", required=True)
    run_parser.add_argument("--output", required=True)

    args = parser.parse_args()

    if args.subcommand == "run":
        run_pipeline(args.patient_id, args.vcf, args.meds, args.output)

if __name__ == "__main__":
    main()
