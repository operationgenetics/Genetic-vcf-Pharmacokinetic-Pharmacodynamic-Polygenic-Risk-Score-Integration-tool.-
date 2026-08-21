def print_clinical_dashboard(report):
    print("=" * 80)
    print("            PRODUCTION CLINICAL GENOMIC INTELLIGENCE DASHBOARD")
    print("=" * 80)
    print(f"[PATIENT ID] : {report.get('patient_id')}")
    print(f"[GENOME BUILD]: {report.get('genome_build')}")
    print(f"[TIMESTAMP]  : {report.get('timestamp')}")
    print(f"[STATUS]     : Certified Production Analysis (Research & Clinical Decision Support)")
    print()
    
    print("--- 1. ACTIVE MEDICATION & OPTIMAL DRUG GUIDANCE ---")
    for med in report.get('active_medication_profile', []):
        print(f" • {med.get('input_name')} (RxCUI: {med.get('rxcui')}) [Class: {med.get('therapeutic_class')}]")
        print(f"   ↳ Pharmacogenetic Impact: {med.get('metabolic_impact')}")
        print(f"   ↳ Optimization Guidance : {med.get('optimal_alternative')}")
        
    print()
    print("--- 2. PHARMACOGENOMICS & METABOLIZER PHENOTYPES ---")
    for gene, pk in report.get('pharmacokinetics', {}).items():
        print(f" • Gene: {gene}")
        print(f"   Diplotype: {pk.get('diplotype')} | Phenotype: {pk.get('phenotype')}")
        print(f"   Clinical Implication: {pk.get('implication')}")
        
    print()
    print("--- 3. UNIVERSAL POLYGENIC RISK SCORES (GLOBAL PGS CATALOG) ---")
    for prs in report.get('polygenic_risk_scores', []):
        print(f" • [{prs.get('pgs_id')}] {prs.get('trait')}")
        print(f"   Risk Percentile: {prs.get('percentile')} percentile -> Tier: {prs.get('risk_tier')}")
        
    if report.get('ddi_warnings'):
        print()
        print("--- 4. DRUG-DRUG INTERACTIONS (DDIs) & SAFETY ALERTS ---")
        for warning in report.get('ddi_warnings'):
            print(f" [!] CRITICAL DDI / SAFETY ALERT: {warning}")
            
    print("=" * 80)
