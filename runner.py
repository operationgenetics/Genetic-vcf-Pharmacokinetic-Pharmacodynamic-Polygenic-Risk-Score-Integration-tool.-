import sqlite3
import json
import argparse
import os
import gzip

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

    phenotypes = {}
    for gene, star, status, pattern in rules:
        if gene not in phenotypes:
            phenotypes[gene] = "Normal Metabolizer (*1/*1)"
        for rsid, var in patient_variants.items():
            if var['genotype'] in ['1/1', '1|1', '0/1', '1/0']:
                if gene == 'CYP2D6' and rsid == 'rs1065852':
                    phenotypes[gene] = f"Poor Metabolizer ({star}/{star})"
                elif gene == 'CYP2C19' and rsid == 'rs4244285':
                    phenotypes[gene] = f"Poor Metabolizer ({star}/{star})"
                elif gene == 'CYP2C9' and rsid == 'rs1057910':
                    phenotypes[gene] = f"Poor Metabolizer ({star}/{star})"

    return phenotypes

def evaluate_prs(conn):
    cursor = conn.cursor()
    cursor.execute("SELECT disease_name, percentile, risk_category, base_score FROM prs_rules")
    return [
        {"disease": row[0], "percentile": row[1], "risk_category": row[2], "score": row[3]}
        for row in cursor.fetchall()
    ]

def evaluate_pk_pd(conn):
    cursor = conn.cursor()
    cursor.execute("SELECT drug_name, gene_symbol, mechanism_type, effect_summary FROM pk_pd_mechanisms")
    return [
        {"drug": row[0], "gene": row[1], "type": row[2], "effect": row[3]}
        for row in cursor.fetchall()
    ]

def evaluate_cpic(conn):
    cursor = conn.cursor()
    cursor.execute("SELECT drug_name, gene_symbol, status, recommendation FROM cpic_guidelines")
    return [
        {"drug": row[0], "gene": row[1], "status": row[2], "recommendation": row[3]}
        for row in cursor.fetchall()
    ]

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

def evaluate_targeted_therapies(conn, pgx_phenotypes):
    cursor = conn.cursor()
    cursor.execute("SELECT condition_trait, primary_drug, gene_checked, pgx_status, alternative_drug, clinical_rationale FROM disease_targeted_therapies")
    therapies = []
    for trait, drug, gene, status, alt_drug, rationale in cursor.fetchall():
        gene_status = pgx_phenotypes.get(gene, "Normal")
        is_high_risk = "Poor" in gene_status or "Decreased" in gene_status

        final_status = "REASSIGNED TO ALTERNATIVE (" + alt_drug + ") due to " + gene + " status (" + status + ")" if is_high_risk and status == "HIGH_RISK" else "SUITABLE"
        selected_drug = alt_drug if is_high_risk and status == "HIGH_RISK" else drug

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
    prs_results = evaluate_prs(conn)
    pk_pd_results = evaluate_pk_pd(conn)
    cpic_results = evaluate_cpic(conn)
    ddi_results = evaluate_ddi_rules(conn)
    acmg_results = evaluate_acmg_findings(conn, patient_variants)
    clinvar_results = evaluate_genome_clinvar(conn, patient_variants)
    targeted_therapies = evaluate_targeted_therapies(conn, pgx_phenotypes)

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
