import re

with open('runner.py', 'r') as f:
    code = f.read()

# 1. Ensure clean, robust parsing of targeted loop and report data block
start_marker = "for trait, p_drug, gene, status, alt_drug, rationale in targeted:"
end_marker = "with open(args.output, 'w') as f:"

if start_marker in code and end_marker in code:
    start_idx = code.find(start_marker)
    end_idx = code.find(end_marker)

    enhanced_section = """for trait, p_drug, gene, status, alt_drug, rationale in targeted:
        prs_info = prs_results.get(trait, {})
        z_score = prs_info.get('z_score', 'N/A')
        percentile = prs_info.get('percentile', 'N/A')
        category = prs_info.get('category', 'N/A')
        
        # Format terminal card output
        print(f"  ┌── [ DISEASE / TRAIT ]: {trait.upper()}")
        print(f"  │   ├── Polygenic Risk   : {category} (Z-Score: {z_score}, Percentile: {percentile}%)")
        print(f"  │   ├── Targeted Drug    : {p_drug}")
        print(f"  │   ├── Gene Evaluated   : {gene} (Patient PGx Status: {status})")
        if any(keyword in str(status).upper() for keyword in ['CONTRAINDICATED', 'HIGH_RISK', 'REASSIGNED']):
            print(f"  │   ├── Action Required  : ⚠️ SWITCH / REASSIGN -> {alt_drug}")
        else:
            print(f"  │   ├── Action Required  : ✔ APPROVED -> Maintain {p_drug}")
        print(f"  │   └── Clinical Rationale: {rationale}")
        print("  └───────────────────────────────────────────────────────────────────\\n")

    report_data = {
        "patient_id": args.patient_id,
        "vcf_file": args.vcf,
        "pgx_phenotypes": pgx_phenotypes if 'pgx_phenotypes' in locals() else {},
        "prs_scores": prs_results if 'prs_results' in locals() else {},
        "targeted_therapies": [
            {
                "trait": trait,
                "polygenic_risk": prs_results.get(trait, {}),
                "primary_drug": p_drug,
                "gene": gene,
                "pgx_status": status,
                "alternative_drug": alt_drug,
                "rationale": rationale
            }
            for trait, p_drug, gene, status, alt_drug, rationale in targeted
        ]
    }

    """

    new_code = code[:start_idx] + enhanced_section + code[end_idx:]

    # Remove any lingering duplicate completion print banners at the bottom of the file
    duplicate_banner = """print("=" * 80)
print("[✔] Complete multi-engine report saved:", args.output)
print("=" * 80)"""

    # Keep exactly one instance of the output banner at the very end
    if new_code.count(duplicate_banner) > 1:
        parts = new_code.split(duplicate_banner)
        new_code = duplicate_banner.join(parts[:-1]) + duplicate_banner + parts[-1]

    with open('runner.py', 'w') as f:
        f.write(new_code)
    print("[✔] runner.py successfully updated and deduplicated.")
else:
    print("[!] Target markers not found in runner.py. Check file state.")
