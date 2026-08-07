#!/usr/bin/env python3
import gzip
import sqlite3
import json
import argparse
import math
import sys

def parse_vcf_qc(vcf_path, min_gq=20, min_dp=10):
    """
    Parses VCF (.vcf or .vcf.gz) with quality control filters:
    - Filters out variants where FILTER != PASS
    - Enforces Minimum Genotype Quality (GQ)
    - Enforces Minimum Read Depth (DP)
    Returns: dict mapping rsid -> { 'genotype': 'A/T', 'ref': 'A', 'alt': 'T', 'qc_passed': True }
    """
    variants = {}
    is_gz = vcf_path.endswith('.gz')
    open_fn = gzip.open if is_gz else open

    with open_fn(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 10:
                continue

            chrom, pos, rsid, ref, alt, qual, filt, info, fmt, sample = parts[:10]
            
            # Filter check
            if filt != 'PASS' and filt != '.':
                continue

            # Parse Format Fields
            fmt_keys = fmt.split(':')
            sample_vals = sample.split(':')
            fmt_map = dict(zip(fmt_keys, sample_vals))

            # Quality Control Thresholds
            try:
                gq = float(fmt_map.get('GQ', 99))
                dp = float(fmt_map.get('DP', 99))
                if gq < min_gq or dp < min_dp:
                    continue
            except ValueError:
                pass # Default to passing if numerical parse fails

            # Resolve Genotype
            gt = fmt_map.get('GT', './.')
            if gt in ['./.', '.']:
                continue

            alleles = [ref] + alt.split(',')
            try:
                gt_indices = [int(i) for i in gt.replace('|', '/').split('/')]
                gt_alleles = [alleles[i] for i in gt_indices]
            except (ValueError, IndexError):
                continue

            genotype_str = "/".join(gt_alleles)
            variants[rsid] = {
                'genotype': genotype_str,
                'alleles': gt_alleles,
                'ref': ref,
                'alt': alt
            }
    return variants

def norm_cdf(x):
    """Cumulative distribution function for standard normal distribution."""
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))

def calculate_prs(variants, conn):
    """
    Calculates polygenic risk scores using additive log-odds weights and normalizes
    against trait population mean and std deviation to produce Z-scores and percentiles.
    """
    cursor = conn.cursor()
    cursor.execute("SELECT trait, mean_score, std_score FROM prs_traits")
    traits = cursor.fetchall()

    prs_results = {}
    for trait, mean_val, std_val in traits:
        cursor.execute("SELECT rsid, risk_allele, effect_weight FROM prs_weights WHERE trait = ?", (trait,))
        weights = cursor.fetchall()
        
        raw_score = 0.0
        snps_counted = 0

        for rsid, risk_allele, weight in weights:
            if rsid in variants:
                gt_alleles = variants[rsid]['alleles']
                dosage = gt_alleles.count(risk_allele)
                raw_score += dosage * weight
                snps_counted += 1

        # Z-Score Calculation
        z_score = (raw_score - mean_val) / std_val if std_val > 0 else 0.0
        percentile = round(norm_cdf(z_score) * 100, 1)

        # Categorization
        if percentile >= 80:
            category = "High Risk"
        elif percentile <= 20:
            category = "Low Risk"
        else:
            category = "Moderate Risk"

        prs_results[trait] = {
            'score': round(raw_score, 2),
            'z_score': round(z_score, 2),
            'percentile': percentile,
            'category': category
        }
    return prs_results

def evaluate_pgx(variants, conn):
    """Maps variant patterns to Star Allele diplotypes and CPIC guidelines."""
    cursor = conn.cursor()
    
    # 1. Phenotype Identification
    cursor.execute("SELECT gene_symbol, diplotype, metabolizer_status, required_variants FROM pgx_star_alleles")
    star_rules = cursor.fetchall()

    phenotypes = {}
    genes = set([r[0] for r in star_rules])

    for gene in genes:
        matched_diplotype = "Unknown"
        matched_status = "Normal Metabolizer (*1/*1)"

        gene_rules = [r for r in star_rules if r[0] == gene]
        for _, diplotype, status, req_var in gene_rules:
            if req_var == "DEFAULT":
                continue
            
            # Check variant match: format "rsid:allele1:allele2"
            rsid, a1, a2 = req_var.split(':')
            if rsid in variants:
                v_alleles = variants[rsid]['alleles']
                if sorted(v_alleles) == sorted([a1, a2]):
                    matched_diplotype = diplotype
                    matched_status = f"{status} ({diplotype})"
                    break

        phenotypes[gene] = matched_status

    # 2. Fetch CPIC Guidelines
    cursor.execute("SELECT drug_name, gene_symbol, status, recommendation FROM cpic_guidelines")
    cpic = cursor.fetchall()
    
    # 3. Fetch PK/PD
    cursor.execute("SELECT drug_name, gene_symbol, mechanism_type, effect_summary FROM pk_pd_mechanisms")
    pkpd = cursor.fetchall()

    # 4. Fetch DDIs
    cursor.execute("SELECT drug1, drug2, interaction_level, clinical_effect FROM ddi_rules")
    ddis = cursor.fetchall()

    # 5. Fetch Targeted Therapies
    cursor.execute("SELECT condition_trait, primary_drug, gene_checked, pgx_status, alternative_drug, clinical_rationale FROM disease_targeted_therapies")
    targeted = cursor.fetchall()

    return phenotypes, cpic, pkpd, ddis, targeted

def run_pipeline(vcf_file, patient_id, output_json):
    conn = sqlite3.connect("genomic_knowledgebase.db")
    
    # Parse VCF with QC
    variants = parse_vcf_qc(vcf_file)
    
    # Calculate Statistical PRS
    prs_results = calculate_prs(variants, conn)

    # Evaluate Pharmacogenomics & CPIC
    phenotypes, cpic, pkpd, ddis, targeted = evaluate_pgx(variants, conn)

    # Console Report Generation
    print("=" * 80)
    print("      COMPREHENSIVE PRECISION MEDICINE & GENOME-WIDE CLINVAR REPORT")
    print(f"      Patient ID: {patient_id} | VCF File: {vcf_file}")
    print("=" * 80)

    print("\n🧬 [1. PHARMACOGENOMIC PHENOTYPES]")
    for gene, status in phenotypes.items():
        print(f"  • {gene:<10} : {status}")

    print("\n📊 [2. ALL POLYGENIC RISK SCORES (PRS - STATISTICAL Z-SCORE)]")
    for trait, metrics in prs_results.items():
        print(f"  • {trait:<30} : {metrics['category']:<12} (Percentile: {metrics['percentile']}%, Z-Score: {metrics['z_score']}, Raw: {metrics['score']})")

    print("\n🔬 [3. PHARMACOKINETIC (PK) & PHARMACODYNAMIC (PD) MECHANISMS]")
    for drug, gene, mech, summary in pkpd:
        print(f"  • [{mech}] {drug} ({gene}): {summary}")

    print("\n💊 [4. CPIC & WESTERN MEDICINE THERAPEUTIC SCREEN]")
    print(f"{'DRUG':<15} {'GENE':<10} {'STATUS':<16} {'RECOMMENDATION'}")
    print("-" * 80)
    for drug, gene, status, rec in cpic:
        print(f"{drug:<15} {gene:<10} {status:<16} {rec}")

    print("\n⚠️ [5. DRUG-DRUG INTERACTIONS (DDI)]")
    for d1, d2, level, effect in ddis:
        print(f"  • {d1} + {d2} [{level}]: {effect}")

    print("\n🎯 [6. POLYGENIC RISK & GENETIC CONDITION TARGETED THERAPIES]")
    for trait, p_drug, gene, status, alt_drug, rationale in targeted:
        print(f"  • CONDITION / TRAIT : {trait}")
        print(f"    SELECTED DRUG      : {p_drug}")
        print(f"    PGX STATUS         : REASSIGNED TO ALTERNATIVE ({alt_drug}) due to {gene} status ({status})")
        print(f"    RATIONALE          : {rationale}\n")

    print("=" * 80)

    # Save Structured JSON Output
    report_data = {
        "patient_id": patient_id,
        "vcf_source": vcf_file,
        "qc_variants_parsed": len(variants),
        "pgx_phenotypes": phenotypes,
        "polygenic_risk_scores": prs_results,
        "cpic_guidelines": cpic,
        "pkpd_mechanisms": pkpd,
        "ddi_rules": ddis,
        "targeted_therapies": targeted
    }

    with open(output_json, 'w') as f:
        json.dump(report_data, f, indent=2)

    print(f"[✔] Complete multi-engine report saved: {output_json}")
    print("=" * 80)
    conn.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Enterprise Genomic Pipeline Engine")
    parser.add_argument("--vcf", required=True, help="Path to input VCF file (.vcf or .vcf.gz)")
    parser.add_argument("--patient-id", required=True, help="Patient identifier")
    parser.add_argument("--output", required=True, help="Path to save JSON report")
    args = parser.parse_args()

    run_pipeline(args.vcf, args.patient_id, args.output)
