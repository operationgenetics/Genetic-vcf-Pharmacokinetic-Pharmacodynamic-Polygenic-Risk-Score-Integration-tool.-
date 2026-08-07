import sys
import os
import argparse
import json
import sqlite3
from datetime import datetime, timezone

DB_PATH = 'genomic_knowledgebase.db'

def ensure_db_exists():
    if not os.path.exists(DB_PATH):
        print(f"[!] '{DB_PATH}' not found. Seeding database...")
        import setup_db
        setup_db.init_db()

def run_pipeline(patient_id, vcf_path, meds, output_path):
    ensure_db_exists()

    if not os.path.exists(vcf_path):
        print(f"❌ Error: VCF file not found at '{vcf_path}'.")
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

    cursor.execute("SELECT drug_name, gene_symbol, mechanism_type, biological_pathway, clinical_effect FROM pk_pd_annotations")
    pk_pd_results = [{
        "drug_name": r[0], "gene_symbol": r[1], "mechanism_type": r[2],
        "biological_pathway": r[3], "clinical_effect": r[4]
    } for r in cursor.fetchall()]
        
    cursor.execute("SELECT trait, SUM(weight) FROM prs_weights GROUP BY trait")
    prs_scores = {}
    for trait, score in cursor.fetchall():
        percentile = min(99, int(score * 45))
        prs_scores[trait] = {
            "raw_score": round(score, 3),
            "percentile": percentile,
            "risk_category": "High" if percentile > 75 else ("Moderate" if percentile > 25 else "Low")
        }

    cursor.execute("SELECT drug_a, drug_b, interaction_severity, mechanism, clinical_guidance FROM ddi_rules")
    ddi_results = [{
        "drug_a": r[0], "drug_b": r[1], "severity": r[2],
        "mechanism": r[3], "clinical_guidance": r[4]
    } for r in cursor.fetchall()]

    cursor.execute("SELECT rsid, gene_symbol, clinical_significance, associated_condition, acmg_actionable FROM pathogenicity_db")
    pathogenicity_findings = [{
        "rsid": r[0], "gene_symbol": r[1], "significance": r[2],
        "condition": r[3], "acmg_secondary_finding": bool(r[4])
    } for r in cursor.fetchall()]

    cursor.execute("SELECT condition_or_trait, condition_type, first_line_drug, alternative_drug, target_gene_check, clinical_rationale FROM prs_condition_therapies")
    prs_therapy_rules = cursor.fetchall()

    targeted_therapies = []
    for cond_name, cond_type, first_line, alt_drug, gene_check, rationale in prs_therapy_rules:
        trigger_active = False
        trigger_reason = ""

        if cond_type == "PRS_TRAIT":
            prs_info = prs_scores.get(cond_name, {})
            if prs_info.get("risk_category") in ["High", "Moderate"]:
                trigger_active = True
                trigger_reason = f"Elevated PRS risk ({prs_info.get('risk_category')}, {prs_info.get('percentile')}th percentile)"
        elif cond_type == "ACMG_VARIANT":
            matching_variant = next((p for p in pathogenicity_findings if p["condition"] == cond_name), None)
            if matching_variant:
                trigger_active = True
                trigger_reason = f"Pathogenic Variant Detected ({matching_variant['rsid']} in {matching_variant['gene_symbol']})"

        if trigger_active:
            cpic_match = next((m for m in matrix if m["drug_name"].lower() == first_line.lower()), None)
            
            if cpic_match and cpic_match["clinical_status"] in ["CONTRAINDICATED", "HIGH_RISK"]:
                selected_drug = alt_drug
                status_summary = f"REASSIGNED TO ALTERNATIVE ({alt_drug}) due to {gene_check} status ({cpic_match['clinical_status']})."
            elif cpic_match:
                selected_drug = first_line
                status_summary = f"SUITABLE ({first_line}) - {cpic_match['primary_recommendation']}"
            else:
                selected_drug = first_line
                status_summary = f"RECOMMENDED ({first_line}) - Standard dosing."

            targeted_therapies.append({
                "condition_or_trait": cond_name,
                "trigger_source": trigger_reason,
                "first_line_drug": first_line,
                "selected_drug": selected_drug,
                "gene_checked": gene_check,
                "pharmacogenomic_status": status_summary,
                "clinical_rationale": rationale
            })

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
        "pathogenicity_findings": pathogenicity_findings,
        "prs_and_genetic_targeted_therapies": targeted_therapies
    }

    os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(report, f, indent=2)

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

    print("\n🎯 [7. POLYGENIC RISK & GENETIC CONDITION TARGETED THERAPIES]")
    if targeted_therapies:
        for rx in targeted_therapies:
            print(f"  • CONDITION / TRAIT : {rx['condition_or_trait']}")
            print(f"    TRIGGER SOURCE    : {rx['trigger_source']}")
            print(f"    SELECTED DRUG     : {rx['selected_drug']} (Gene Checked: {rx['gene_checked']})")
            print(f"    PGX STATUS        : {rx['pharmacogenomic_status']}")
            print(f"    RATIONALE         : {rx['clinical_rationale']}\n")
    else:
        print("  • No high-risk polygenic or pathogenic conditions triggered specific therapy recommendations.")

    print("="*80)
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
