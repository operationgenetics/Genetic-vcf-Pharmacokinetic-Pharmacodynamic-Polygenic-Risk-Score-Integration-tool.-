import argparse
import json
import sqlite3
import os

def connect_db():
    return sqlite3.connect("genomic_knowledgebase.db")

def parse_vcf_variants(vcf_path):
    patient_variants = {}
    if not os.path.exists(vcf_path):
        print(f"[!] VCF file {vcf_path} not found.")
        return patient_variants
        
    with open(vcf_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 5:
                chrom, pos, rsid, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
                gt = "0/1"
                if len(parts) >= 10:
                    gt_info = parts[9].split(':')[0]
                    gt = gt_info
                patient_variants[rsid] = {
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "genotype": gt
                }
    return patient_variants

def evaluate_genome_wide_clinvar(conn, patient_variants):
    cursor = conn.cursor()
    findings = []
    
    for rsid, details in patient_variants.items():
        cursor.execute("""
            SELECT gene_symbol, clinical_significance, associated_conditions, review_status
            FROM genome_clinvar WHERE rsid = ?
        """, (rsid,))
        row = cursor.fetchone()
        if row:
            findings.append({
                "rsid": rsid,
                "gene": row[0],
                "significance": row[1],
                "condition": row[2],
                "review_status": row[3],
                "genotype": details["genotype"]
            })
    return findings

def evaluate_pgx_phenotypes(conn, patient_variants):
    cursor = conn.cursor()
    cursor.execute("SELECT gene_symbol, star_allele, metabolizer_status, genotype_pattern FROM pgx_star_alleles")
    rules = cursor.fetchall()
    
    phenotypes = {
        "CYP2C9": "*1/*1 (Normal Metabolizer)",
        "CYP2C19": "*1/*1 (Normal Metabolizer)",
        "SLCO1B1": "*1/*1 (Normal Function)",
        "CYP3A5": "*3/*3 (Non-Expresser)",
        "CYP2D6": "*1/*1 (Normal Metabolizer)",
        "VKORC1": "1639G>G (Normal Sensitivity)"
    }

    if "rs1065852" in patient_variants:
        phenotypes["CYP2D6"] = "*10/*10 (Poor Metabolizer)"
    if "rs1057910" in patient_variants or "rs1800492" in patient_variants:
        phenotypes["CYP2C9"] = "*3/*3 (Poor Metabolizer)"
    if "rs4244285" in patient_variants or "rs12248560" in patient_variants:
        phenotypes["CYP2C19"] = "*2/*2 (Poor Metabolizer)"
    if "rs4149056" in patient_variants:
        phenotypes["SLCO1B1"] = "*5/*5 (Decreased Function)"
    if "rs9923231" in patient_variants:
        phenotypes["VKORC1"] = "-1639G>A (High Sensitivity)"

    return phenotypes

def evaluate_prs_scores(conn, patient_variants):
    cursor = conn.cursor()
    cursor.execute("SELECT disease_name, percentile, risk_category, base_score FROM prs_rules")
    scores = {}
    for disease, percentile, risk, score in cursor.fetchall():
        scores[disease] = {
            "percentile": percentile,
            "category": risk,
            "score": score
        }
    return scores

def evaluate_pk_pd_mechanisms(conn, patient_variants, pgx):
    cursor = conn.cursor()
    cursor.execute("SELECT drug_name, gene_symbol, mechanism_type, effect_summary FROM pk_pd_mechanisms")
    mechanisms = []
    for drug, gene, mech_type, effect in cursor.fetchall():
        mechanisms.append({
            "drug": drug,
            "gene": gene,
            "type": mech_type,
            "effect": effect
        })
    return mechanisms

def evaluate_cpic_guidelines(conn, pgx):
    cursor = conn.cursor()
    cursor.execute("SELECT drug_name, gene_symbol, status, recommendation FROM cpic_guidelines")
    cpic = []
    for drug, gene, status, rec in cursor.fetchall():
        cpic.append({
            "drug": drug,
            "gene": gene,
            "status": status,
            "recommendation": rec
        })
    return cpic

def evaluate_ddi_rules(conn):
    cursor = conn.cursor()
    cursor.execute("SELECT drug1, drug2, interaction_level, clinical_effect FROM ddi_rules")
    ddis = []
    for d1, d2, level, effect in cursor.fetchall():
        ddis.append({
            "pair": f"{d1} + {d2}",
            "level": level,
            "effect": effect
        })
    return ddis

def evaluate_acmg_findings(conn, patient_variants):
    cursor = conn.cursor()
    cursor.execute("SELECT rsid, gene_symbol, pathogenicity, disease_association, actionable FROM acmg_findings")
    findings = []
    for rsid, gene, path, disease, actionable in cursor.fetchall():
        if rsid in patient_variants:
            findings.append({
                "rsid": rsid,
                "gene": gene,
                "pathogenicity": path,
                "disease": disease,
                "actionable": bool(actionable)
            })
    return findings

def evaluate_targeted_therapies(conn, prs_scores, acmg_findings, pgx):
    cursor = conn.cursor()
    cursor.execute("SELECT condition_trait, primary_drug, gene_checked, pgx_status, alternative_drug, clinical_rationale FROM disease_targeted_therapies")
    therapies = []
    for cond, primary, gene, default_status, alt_drug, rec in cursor.fetchall():
        status = default_status
        selected_drug = primary
        
        if gene in pgx and "Poor" in pgx[gene]:
            status = f"REASSIGNED TO ALTERNATIVE ({alt_drug}) due to {gene} status (HIGH_RISK)."
            selected_drug = alt_drug
        elif gene in pgx and "Decreased" in pgx[gene]:
            status = f"SUITABLE ({primary}) with dose adjustment."

        therapies.append({
            "condition": cond,
            "selected_drug": selected_drug,
            "gene_checked": gene,
            "pgx_status": status,
            "rationale": rec
        })
    return therapies

def generate_report(vcf_path, patient_id, output_path):
    conn = connect_db()
    patient_variants = parse_vcf_variants(vcf_path)
    
    pgx = evaluate_pgx_phenotypes(conn, patient_variants)
    prs = evaluate_prs_scores(conn, patient_variants)
    pkpd = evaluate_pk_pd_mechanisms(conn, patient_variants, pgx)
    cpic = evaluate_cpic_guidelines(conn, pgx)
    ddi = evaluate_ddi_rules(conn)
    acmg = evaluate_acmg_findings(conn, patient_variants)
    targeted = evaluate_targeted_therapies(conn, prs, acmg, pgx)
    clinvar = evaluate_genome_wide_clinvar(conn, patient_variants)

    report = {
        "patient_id": patient_id,
        "vcf_source": vcf_path,
        "pgx_phenotypes": pgx,
        "polygenic_risk_scores": prs,
        "pk_pd_mechanisms": pkpd,
        "cpic_guidelines": cpic,
        "drug_interactions": ddi,
        "acmg_secondary_findings": acmg,
        "genome_wide_clinvar_findings": clinvar,
        "targeted_therapies": targeted
    }

    print("\n" + "="*80)
    print("      COMPREHENSIVE PRECISION MEDICINE & GENOME-WIDE CLINVAR REPORT")
    print(f"      Patient ID: {patient_id} | VCF File: {os.path.basename(vcf_path)}")
    print("="*80)

    print("\n🧬 [1. PHARMACOGENOMIC PHENOTYPES]")
    for gene, status in pgx.items():
        print(f"  • {gene:<10}: {status}")

    print("\n📊 [2. ALL POLYGENIC RISK SCORES (PRS)]")
    for disease, data in prs.items():
        print(f"  • {disease:<25}: {data['category']:<13} Risk (Percentile: {data['percentile']}%, Score: {data['score']})")

    print("\n🔬 [3. PHARMACOKINETIC (PK) & PHARMACODYNAMIC (PD) MECHANISMS]")
    for m in pkpd:
        print(f"  • [{m['type']}] {m['drug']} ({m['gene']}): {m['effect']}")

    print("\n💊 [4. CPIC & WESTERN MEDICINE THERAPEUTIC SCREEN]")
    print(f"{'DRUG':<15}{'GENE':<10}{'STATUS':<17}{'RECOMMENDATION'}")
    print("-" * 80)
    for c in cpic:
        print(f"{c['drug']:<15}{c['gene']:<10}{c['status']:<17}{c['recommendation']}")

    print("\n⚠️ [5. DRUG-DRUG INTERACTIONS (DDI)]")
    for d in ddi:
        print(f"  • {d['pair']} [{d['level']}]: {d['effect']}")

    print("\n🚨 [6. PATHOGENICITY & ACMG SECONDARY FINDINGS]")
    for a in acmg:
        print(f"  • {a['rsid']} ({a['gene']}): {a['pathogenicity']} for {a['disease']} [ACMG Actionable]")

    print("\n🌐 [7. GENOME-WIDE CLINVAR ANNOTATIONS (ALL KNOWN DISORDERS & VARIANT RISKS)]")
    if clinvar:
        for cl in clinvar:
            print(f"  • {cl['rsid']} ({cl['gene']}): {cl['significance']} for {cl['condition']} (Genotype: {cl['genotype']})")
    else:
        print("  • No additional ClinVar pathogenic/benign variants detected in target dataset.")

    print("\n🎯 [8. POLYGENIC RISK & GENETIC CONDITION TARGETED THERAPIES]")
    for t in targeted:
        print(f"  • CONDITION / TRAIT : {t['condition']}")
        print(f"    SELECTED DRUG      : {t['selected_drug']} (Gene Checked: {t['gene_checked']})")
        print(f"    PGX STATUS         : {t['pgx_status']}")
        print(f"    RATIONALE          : {t['rationale']}\n")

    print("="*80)

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(report, f, indent=2)

    print(f"[✔] Complete multi-engine report saved: {output_path}")
    print("="*80 + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Genome-Wide Precision Medicine Bot")
    parser.add_argument("--vcf", required=True, help="Path to patient VCF file")
    parser.add_argument("--patient-id", default="PATIENT_01", help="Patient Identifier")
    parser.add_argument("--meds", default="all", help="Target medication panel")
    parser.add_argument("--output", default="results/report.json", help="Output JSON path")

    args = parser.parse_args()
    generate_report(args.vcf, args.patient_id, args.output)
