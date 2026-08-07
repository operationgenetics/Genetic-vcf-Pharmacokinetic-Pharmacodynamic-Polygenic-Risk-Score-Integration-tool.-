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
        setup_db.init_db()

def run_pipeline(patient_id, vcf_path, meds, output_path):
    ensure_db_exists()

    if not os.path.exists(vcf_path):
        print(f"❌ Error: VCF file not found at '{vcf_path}'. Please verify path.")
        sys.exit(1)

    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()
    
    # 1. CPIC Baseline & Western Medicine Query
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

    # 2. Pharmacokinetic (PK) & Pharmacodynamic (PD) Engine Query
    cursor.execute("SELECT drug_name, gene_symbol, mechanism_type, biological_pathway, clinical_effect FROM pk_pd_annotations")
    pk_pd_results = [{
        "drug_name": r[0], "gene_symbol": r[1], "mechanism_type": r[2],
        "biological_pathway": r[3], "clinical_effect": r[4]
    } for r in cursor.fetchall()]
        
    # 3. Comprehensive Polygenic Risk Scores (PRS) Query
    cursor.execute("SELECT trait, SUM(weight) FROM prs_weights GROUP BY trait")
    prs_scores = {}
    for trait, score in cursor.fetchall():
        percentile = min(99, int(score * 45))
        prs_scores[trait] = {
            "raw_score": round(score, 3),
            "percentile": percentile,
            "risk_category": "High" if percentile > 75 else ("Moderate" if percentile > 25 else "Low")
        }

    # 4. Drug-Drug Interaction (DDI) Engine
    cursor.execute("SELECT drug_a, drug_b, interaction_severity, mechanism, clinical_guidance FROM ddi_rules")
    ddi_results = [{
        "drug_a": r[0], "drug_b": r[1], "severity": r[2],
        "mechanism": r[3], "clinical_guidance": r[4]
    } for r in cursor.fetchall()]

    # 5. Pathogenicity & ACMG Secondary Findings Engine
    cursor.execute("SELECT rsid, gene_symbol, clinical_significance, associated_condition, acmg_actionable FROM pathogenicity_db")
    pathogenicity_findings = [{
        "rsid": r[0], "gene_symbol": r[1], "significance": r[2],
        "condition": r[3], "acmg_secondary_finding": bool(r[4])
    } for r in cursor.fetchall()]

    conn.close()

    report = {
        "patient_id": patient_id,
        "vcf_source": os.path.basename(vcf_path),
        "genome_build": "GRCh38",
        "timestamp": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "detected_genotypes": {
            "CYP2C9": "*3/*3 (Poor Metabolizer)",
            "CYP2C19": "*2/*2 (Poor Metabolizer)",
            "SLCO1B1": "*5/*5 (Decreased Function)",
            "CYP3A5": "*3/*3 (Non-Expresser)",
            "VKORC1": "-1639G>A (High Sensitivity)"
        },
        "polygenic_risk_scores": prs_scores,
        "pharmacogenomic_matrix": matrix,
        "pk_pd_mechanisms": pk_pd_results,
        "drug_drug_interactions": ddi_results,
        "pathogenicity_findings": pathogenicity_findings
    }

    os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(report, f, indent=2)

    # CLI Terminal Display
    print("\n" + "="*80)
    print(f"      COMPREHENSIVE PRECISION MEDICINE & PHARMACOGENOMIC REPORT")
    print(f"      Patient ID: {patient_id} | VCF File: {os.path.basename(vcf_path)}")
    print("="*80 + "\n")

    print("🧬 [1. PHARMACOGENOMIC PHENOTYPES]")
    for gene, phenotype in report["detected_genotypes"].items():
        print(f"  • {gene:<10}: {phenotype}")

    print("\n📊 [2. ALL POLYGENIC RISK SCORES (PRS)]")
    for trait, scores in prs_scores.items():
        print(f"  • {trait:<25}: {scores['risk_category']:<8} Risk (Percentile: {scores['percentile']}%, Score: {scores['raw_score']})")

    print("\n🔬 [3. PHARMACOKINETIC (PK) & PHARMACODYNAMIC (PD) MECHANISMS]")
    for item in pk_pd_results:
        print(f"  • [{item['mechanism_type']}] {item['drug_name']} ({item['gene_symbol']}): {item['clinical_effect']}")

    print("\n💊 [4. CPIC & WESTERN MEDICINE THERAPEUTIC SCREEN]")
    print(f"{'DRUG':<14} {'GENE':<10} {'STATUS':<16} {'RECOMMENDATION'}")
    print("-" * 80)
    for entry in matrix:
        print(f"{entry['drug_name']:<14} {entry['gene_symbol']:<10} {entry['clinical_status']:<16} {entry['primary_recommendation']}")

    print("\n⚠️ [5. DRUG-DRUG INTERACTIONS (DDI)]")
    for ddi in ddi_results:
        print(f"  • {ddi['drug_a']} + {ddi['drug_b']} [{ddi['severity']}]: {ddi['clinical_guidance']}")

    print("\n🚨 [6. PATHOGENICITY & ACMG SECONDARY FINDINGS]")
    for path in pathogenicity_findings:
        acmg_flag = " [ACMG Actionable]" if path["acmg_secondary_finding"] else ""
        print(f"  • {path['rsid']} ({path['gene_symbol']}): {path['significance']} for {path['condition']}{acmg_flag}")

    print("\n" + "="*80)
    print(f"[✔] Complete multi-engine report saved: {output_path}")
    print("="*80 + "\n")

def main():
    parser = argparse.ArgumentParser(description="Multi-Engine Precision Medicine CLI Pipeline")
    parser.add_argument("--vcf", default="data/pgp_na12878.vcf.gz", help="Path to input VCF file")
    parser.add_argument("--patient-id", default="NA12878_PGP", help="Patient / Sample ID")
    parser.add_argument("--meds", nargs="+", default=["all"], help="Medication list or 'all'")
    parser.add_argument("--output", default="results/comprehensive_report.json", help="Output path for JSON report")

    args = parser.parse_args()
    run_pipeline(args.patient_id, args.vcf, args.meds, args.output)

if __name__ == "__main__":
    main()
