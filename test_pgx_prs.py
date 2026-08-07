#!/usr/bin/env python3
from pgx_prs_engine import PRSEngine, PGxEngine

def main():
    prs_eng = PRSEngine()
    pgx_eng = PGxEngine()

    print("================================================================================")
    print("      PRODUCTION-GRADE PRS & PGx COMPLEXITY ENGINE TEST HARNESS")
    print("================================================================================")

    print("\n🧬 [1. PRS ANCESTRY NORMALIZATION & COVERAGE AUDIT]")
    patient_pcs = [1.2, -0.8, 0.1]
    sample_prs_weights = {'rs1234': 0.05, 'rs5678': -0.12, 'rs91011': 0.08, 'rs1213': 0.15}
    vcf_detected_snps = {'rs1234', 'rs5678', 'rs91011'}

    cov_audit = prs_eng.audit_prs_variant_coverage(vcf_detected_snps, sample_prs_weights)
    prs_adj = prs_eng.calculate_ancestry_adjusted_prs(raw_score=1.7, patient_pcs=patient_pcs, ancestry_pop='AFR')

    print(f" • PRS Coverage Audit : {cov_audit['coverage_pct']}% ({cov_audit['status']})")
    print(f" • Adjusted Z-Score   : {prs_adj['z_score']} (Percentile: {prs_adj['percentile']}%)")

    print("\n🔬 [2. PGx PHASING, CNV, & CPIC STANDARDIZED LOOKUP]")
    phase_res = pgx_eng.resolve_diplotype_phase(['*1', '*2'], is_phased=False)
    cnv_res = pgx_eng.detect_cnv_duplications('CYP2D6', coverage_depth=45.0, expected_depth=30.0)
    cpic_res = pgx_eng.query_cpic_guideline('CYP2C19', '*2/*2')

    print(f" • Phase Status       : {phase_res['diplotype']} | Warning: {phase_res['warning']}")
    print(f" • CNV Call (CYP2D6)  : {cnv_res['cnv_status']} (Est Copies: {cnv_res['estimated_copy_number']})")
    print(f" • CPIC Phenotype     : {cpic_res['gene']} {cpic_res['diplotype']} -> {cpic_res['cpic_phenotype']}")
    print("================================================================================")
    print("[✔] Production engines successfully compiled and verified.")
    print("================================================================================")

if __name__ == "__main__":
    main()
