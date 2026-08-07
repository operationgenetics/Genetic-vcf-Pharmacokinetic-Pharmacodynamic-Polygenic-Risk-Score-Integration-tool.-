with open('runner.py', 'r') as f:
    content = f.read()

# Locate the targeted loop block and replace it cleanly
target_str = "for trait, p_drug, gene, status, alt_drug, rationale in targeted:"

new_loop_code = """for trait, p_drug, gene, status, alt_drug, rationale in targeted:
        prs_info = prs_results.get(trait, {})
        z_score = prs_info.get('z_score', 'N/A')
        percentile = prs_info.get('percentile', 'N/A')
        category = prs_info.get('category', 'N/A')
        
        print(f"  ┌── [ DISEASE / TRAIT ]: {trait.upper()}")
        print(f"  │   ├── Polygenic Risk   : {category} (Z-Score: {z_score}, Percentile: {percentile}%)")
        print(f"  │   ├── Targeted Drug    : {p_drug}")
        print(f"  │   ├── Gene Evaluated   : {gene} (Patient PGx Status: {status})")
        if 'CONTRAINDICATED' in str(status) or 'HIGH_RISK' in str(status) or 'REASSIGNED' in str(status):
            print(f"  │   ├── Action Required  : ⚠️ SWITCH / REASSIGN -> {alt_drug}")
        else:
            print(f"  │   ├── Action Required  : ✔ APPROVED -> Maintain {p_drug}")
        print(f"  │   └── Clinical Rationale: {rationale}")
        print("  └───────────────────────────────────────────────────────────────────\\n")"""

# Find where the old loop starts and replace through the print statements
if target_str in content:
    start_idx = content.find(target_str)
    # Find end of old print loop before section divider or next block
    end_marker = "print(\"=\" * 80)"
    end_idx = content.find(end_marker, start_idx)
    
    if end_idx != -1:
        updated_content = content[:start_idx] + new_loop_code + "\n\n" + content[end_idx:]
        with open('runner.py', 'w') as f:
            f.write(updated_content)
        print("[✔] runner.py successfully updated with detailed clinical layout.")
    else:
        print("[!] Could not locate end marker in runner.py.")
else:
    print("[!] Could not locate target loop in runner.py.")
