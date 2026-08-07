import sqlite3
import json
import argparse
import os
import gzip
import math

def open_vcf(vcf_path):
    """Seamlessly open both uncompressed .vcf and bgzip/gzip compressed .vcf.gz files."""
    if not os.path.exists(vcf_path):
        print(f"[!] VCF file {vcf_path} not found.")
        return None
    if vcf_path.endswith('.gz'):
        return gzip.open(vcf_path, 'rt', encoding='utf-8', errors='replace')
    return open(vcf_path, 'r', encoding='utf-8')

def parse_vcf(vcf_path, target_sample_id=None):
    """Extract variants (rsIDs and genotypes) from a VCF or VCF.GZ file."""
    variants = {}
    f = open_vcf(vcf_path)
    if not f:
        return variants

    try:
        sample_index = None
        for line in f:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                headers = line.strip().split('\t')
                if len(headers) > 9:
                    if target_sample_id and target_sample_id in headers:
                        sample_index = headers.index(target_sample_id)
                    else:
                        sample_index = 9  # Default to first sample
                continue

            parts = line.strip().split('\t')
            if len(parts) < 8:
                continue

            chrom, pos, rsid, ref, alt, qual, filter_val, info = parts[:8]
            genotype = "./."

            if sample_index and len(parts) > sample_index:
                format_fields = parts[8].split(':')
                sample_fields = parts[sample_index].split(':')
                if 'GT' in format_fields:
                    gt_idx = format_fields.index('GT')
                    if gt_idx < len(sample_fields):
                        genotype = sample_fields[gt_idx]

            if rsid != '.' and rsid:
                variants[rsid] = {
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "genotype": genotype
                }
    finally:
        f.close()

    return variants

def evaluate_pgx_star_alleles(conn, patient_variants):
    cursor = conn.cursor()
    cursor.execute("SELECT gene_symbol, star_allele, metabolizer_status, genotype_pattern FROM pgx_star_alleles")
    rules = cursor.fetchall()

    phenotypes = {
        'CYP2D6': 'Normal Metabolizer (*1/*1)',
        'CYP2C19': 'Normal Metabolizer (*1/*1)',
        'CYP2C9': 'Normal Metabolizer (*1/*1)',
        'SLCO1B1': 'Normal Metabolizer (*1/*1)',
        'CYP3A5': 'Normal Metabolizer (*1/*1)',
        'VKORC1': 'Normal Metabolizer (*1/*1)'
    }

    for gene, star, status, pattern in rules:
        for rsid, var in patient_variants.items():
            gt = var['genotype']
            if gt in ['1/1', '1|1']:
                if gene == 'CYP2D6' and rsid == 'rs1065852':
                    phenotypes[gene] = f"Poor Metabolizer ({star}/{star})"
                elif gene == 'CYP2C19' and rsid == 'rs4244285':
                    phenotypes[gene] = f"Poor Metabolizer ({star}/{star})"
                elif gene == 'CYP2C9' and rsid == 'rs1057910':
                    phenotypes[gene] = f"Poor Metabolizer ({star}/{star})"
            elif gt in ['0/1', '1/0', '0|1', '1|0']:
                if gene == 'CYP2D6' and rsid == 'rs1065852':
                    phenotypes[gene] = f"Intermediate Metabolizer (*1/{star})"
                elif gene == 'CYP2C19' and rsid == 'rs4244285':
                    phenotypes[gene] = f"Intermediate Metabolizer (*1/{star})"
                elif gene == 'CYP2C9' and rsid == 'rs1057910':
                    phenotypes[gene] = f"Intermediate Metabolizer (*1/{star})"

    return phenotypes

def evaluate_prs_dynamic(conn, patient_variants):
    """Calculates Polygenic Risk Scores dynamically from VCF variants using prs_weights if present."""
    cursor = conn.cursor()
    
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='prs_weights'")
    has_prs_weights = cursor.fetchone() is not None

    if not has_prs_weights:
        cursor.execute("SELECT disease_name, percentile, risk_category, base_score FROM prs_rules")
        return [
            {"disease": row[0], "percentile": row[1], "risk_category": row[2], "score": row[3]}
            for row in cursor.fetchall()
        ]

    cursor.execute("SELECT trait, rsid, risk_allele, weight FROM prs_weights")
    rows = cursor.fetchall()

    if not rows:
        cursor.execute("SELECT disease_name, percentile, risk_category, base_score FROM prs_rules")
        return [
            {"disease": row[0], "percentile": row[1], "risk_category": row[2], "score": row[3]}
            for row in cursor.fetchall()
        ]

    trait_scores = {}
    trait_max = {}

    for trait, rsid, risk_allele, weight in rows:
        if trait not in trait_scores:
            trait_scores[trait] = 0.0
            trait_max[trait] = 0.0

        trait_max[trait] += 2.0 * weight

        if rsid in patient_variants:
            gt = patient_variants[rsid].get('genotype', './.')
            alt = patient_variants[rsid].get('alt', '')
            
            dosage = 0
            if gt in ['1/1', '1|1']:
                dosage = 2
            elif gt in ['0/1', '1/0', '0|1', '1|0']:
                dosage = 1

            if alt == risk_allele or dosage > 0:
                trait_scores[trait] += dosage * weight

    results = []
    for trait, raw_score in trait_scores.items():
        max_possible = trait_max.get(trait, 1.0)
        norm_ratio = raw_score / max_possible if max_possible > 0 else 0.0
        
        percentile = min(99, max(1, int(100 / (1 + math.exp(-5 * (norm_ratio - 0.45))))))
        
        if percentile >= 75:
            category = "High"
        elif percentile >= 35:
            category = "Moderate"
        else:
            category = "Low"

        results.append({
            "disease": trait,
            "percentile": percentile,
            "risk_category": category,
            "score": round(raw_score, 2)
        })

    cursor.execute("SELECT disease_name, percentile, risk_category, base_score FROM prs_rules")
    existing_traits = {r["disease"] for r in results}
    for row in cursor.fetchall():
        if row[0] not in existing_traits:
            results.append({
                "disease": row[0],
                "percentile": row[1],
                "risk_category": row[2],
                "score": row[3]
            })

    return results

def evaluate_pk_pd(conn):
    cursor = conn.cursor()
    cursor.execute("SELECT drug_name, gene_symbol, mechanism_type, effect_summary FROM pk_pd_mechanisms")
    return [
        {"drug": row[0], "gene": row[1], "type": row[2], "effect": row[3]}
        for row in cursor.fetchall()
    ]

def evaluate_cpic(conn, pgx_phenotypes):
    cursor = conn.cursor()
    cursor.execute("SELECT drug_name, gene_symbol, status, recommendation FROM cpic_guidelines")
    guidelines = []
    
    for drug, gene, status, rec in cursor.fetchall():
        gene_status = pgx_phenotypes.get(gene, "Normal Metabolizer")
        
        if "Poor Metabolizer" in gene_status:
            if drug in ["Clopidogrel", "Escitalopram", "Carbamazepine", "Fluorouracil"]:
                status = "CONTRAINDICATED"
            elif drug in ["Warfarin", "Aripiprazole", "Risperidone"]:
                status = "HIGH_RISK"
        
        guidelines.append({
            "drug": drug,
            "gene": gene,
            "status": status,
            "recommendation": rec
        })
    return guidelines

def evaluate_ddi_rules(conn):
    cursor = conn.cursor()
    cursor.execute("SELECT drug1, drug2, interaction_level, clinical_effect FROM ddi_rules")
    return [
        {"drug1": row[0], "drug2": row[1], "level": row[2], "effect": row[3]}
        for row in cursor.fetchall()
    ]

def evaluate_acmg_findings(conn, patient_variants):
    cursor = conn.cursor()
    cursor.execute("SELECT rsid, gene_symbol, pathogenicity, disease_association, actionable FROM acmg_findings")
    findings = []
    for rsid, gene, path, disease, actionable in cursor.fetchall():
        if rsid in patient_variants:
            gt = patient_variants[rsid]['genotype']
            if gt in ['1/1', '0/1', '1/0', '1|1', '0|1']:
                findings.append({
                    "rsid": rsid,
                    "gene": gene,
                    "pathogenicity": path,
                    "disease": disease,
                    "actionable": bool(actionable),
                    "genotype": gt
                })
    return findings

def evaluate_genome_clinvar(conn, patient_variants):
    cursor = conn.cursor()
    cursor.execute("SELECT rsid, gene_symbol, clinical_significance, associated_conditions, review_status FROM genome_clinvar")
    clinvar = []
    for rsid, gene, sig, disease, status in cursor.fetchall():
        if rsid in patient_variants:
            clinvar.append({
                "rsid": rsid,
                "gene": gene,
                "significance": sig,
                "disease": disease,
                "review_status": status,
                "genotype": patient_variants[rsid]['genotype']
            })
    return clinvar

def evaluate_targeted_therapies(conn, pgx_phenotypes, prs_results):
    cursor = conn.cursor()
    cursor.execute("SELECT condition_trait, primary_drug, gene_checked, pgx_status, alternative_drug, clinical_rationale FROM disease_targeted_therapies")
    therapies = []
    
    prs_map = {p["disease"]: p for p in prs_results}

    for trait, drug, gene, status, alt_drug, rationale in cursor.fetchall():
        gene_status = pgx_phenotypes.get(gene, "Normal")
        is_high_risk_pgx = "Poor" in gene_status or "Decreased" in gene_status
        
        trait_prs = prs_map.get(trait, {})
        is_high_prs = trait_prs.get("risk_category") == "High"

        if is_high_risk_pgx or (is_high_prs and status == "HIGH_RISK"):
            final_status = f"REASSIGNED TO ALTERNATIVE ({alt_drug}) due to {gene} status ({status})"
            selected_drug = alt_drug
        else:
            final_status = "SUITABLE"
            selected_drug = drug

        therapies.append({
            "condition": trait,
            "primary_drug": drug,
            "selected_drug": selected_drug,
            "gene_checked": gene,
            "pgx_status": final_status,
            "rationale": rationale
        })
    return therapies

def generate_report(vcf_path, patient_id, output_json=None):
    if not os.path.exists("genomic_knowledgebase.db"):
        print("[!] Database genomic_knowledgebase.db missing. Running setup_db.py...")
        os.system("python3 setup_db.py")

    conn = sqlite3.connect("genomic_knowledgebase.db")

    patient_variants = parse_vcf(vcf_path, patient_id)
    pgx_phenotypes = evaluate_pgx_star_alleles(conn, patient_variants)
    prs_results = evaluate_prs_dynamic(conn, patient_variants)
    pk_pd_results = evaluate_pk_pd(conn)
    cpic_results = evaluate_cpic(conn, pgx_phenotypes)
    ddi_results = evaluate_ddi_rules(conn)
    acmg_results = evaluate_acmg_findings(conn, patient_variants)
    clinvar_results = evaluate_genome_clinvar(conn, patient_variants)
    targeted_therapies = evaluate_targeted_therapies(conn, pgx_phenotypes, prs_results)

    report_data = {
        "patient_id": patient_id,
        "vcf_file": vcf_path,
        "variants_parsed": len(patient_variants),
        "pgx_phenotypes": pgx_phenotypes,
        "prs_results": prs_results,
        "pk_pd_mechanisms": pk_pd_results,
        "cpic_guidelines": cpic_results,
        "ddi_interactions": ddi_results,
        "acmg_findings": acmg_results,
        "genome_clinvar": clinvar_results,
        "targeted_therapies": targeted_therapies
    }

    # Console Output
    print("="*80)
    print("      COMPREHENSIVE PRECISION MEDICINE & GENOME-WIDE CLINVAR REPORT")
    print(f"      Patient ID: {patient_id} | VCF File: {vcf_path}")
    print("="*80 + "\n")

    print("🧬 [1. PHARMACOGENOMIC PHENOTYPES]")
    for gene, status in pgx_phenotypes.items():
        print(f"  • {gene:<10}: {status}")

    print("\n📊 [2. ALL POLYGENIC RISK SCORES (PRS)]")
    for prs in prs_results:
        print(f"  • {prs['disease']:<27}: {prs['risk_category']:<13} Risk (Percentile: {prs['percentile']}%, Score: {prs['score']})")

    print("\n🔬 [3. PHARMACOKINETIC (PK) & PHARMACODYNAMIC (PD) MECHANISMS]")
    for pkpd in pk_pd_results:
        print(f"  • [{pkpd['type']}] {pkpd['drug']} ({pkpd['gene']}): {pkpd['effect']}")

    print("\n💊 [4. CPIC & WESTERN MEDICINE THERAPEUTIC SCREEN]")
    print(f"{'DRUG':<15} {'GENE':<10} {'STATUS':<16} {'RECOMMENDATION'}")
    print("-" * 80)
    for cpic in cpic_results:
        print(f"{cpic['drug']:<15} {cpic['gene']:<10} {cpic['status']:<16} {cpic['recommendation']}")

    print("\n⚠️ [5. DRUG-DRUG INTERACTIONS (DDI)]")
    for ddi in ddi_results:
        print(f"  • {ddi['drug1']} + {ddi['drug2']} [{ddi['level']}]: {ddi['effect']}")

    print("\n🚨 [6. PATHOGENICITY & ACMG SECONDARY FINDINGS]")
    if acmg_results:
        for acmg in acmg_results:
            print(f"  • [{acmg['gene']}] {acmg['rsid']} ({acmg['disease']}) - {acmg['pathogenicity']} (GT: {acmg['genotype']})")
    else:
        print("  • No pathogenic secondary findings detected.")

    print("\n🌐 [7. GENOME-WIDE CLINVAR ANNOTATIONS]")
    if clinvar_results:
        for clinvar in clinvar_results:
            print(f"  • [{clinvar['gene']}] {clinvar['rsid']} - {clinvar['disease']} ({clinvar['significance']})")
    else:
        print("  • No additional ClinVar pathogenic/benign variants detected in target dataset.")

    print("\n🎯 [8. POLYGENIC RISK & GENETIC CONDITION TARGETED THERAPIES]")
    for therapy in targeted_therapies:
        print(f"  • CONDITION / TRAIT : {therapy['condition']}")
        print(f"    SELECTED DRUG      : {therapy['selected_drug']} (Gene Checked: {therapy['gene_checked']})")
        print(f"    PGX STATUS         : {therapy['pgx_status']}")
        print(f"    RATIONALE          : {therapy['rationale']}\n")

    print("="*80)

    if output_json:
        os.makedirs(os.path.dirname(output_json) or '.', exist_ok=True)
        with open(output_json, 'w') as f:
            json.dump(report_data, f, indent=2)
        print(f"[✔] Complete multi-engine report saved: {output_json}")
        print("="*80 + "\n")

    conn.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Precision Medicine & Genome-Wide ClinVar Analysis Pipeline")
    parser.add_argument("--vcf", required=True, help="Path to patient VCF or VCF.GZ file")
    parser.add_argument("--patient-id", required=True, help="Target sample ID in VCF")
    parser.add_argument("--meds", default="all", help="Comma-separated list of medications or 'all'")
    parser.add_argument("--output", default="results/patient_report.json", help="Output path for JSON report")

    args = parser.parse_args()
    generate_report(args.vcf, args.patient_id, args.output)
