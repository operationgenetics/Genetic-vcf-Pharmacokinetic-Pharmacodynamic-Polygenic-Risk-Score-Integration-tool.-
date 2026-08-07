with open('runner.py', 'r') as f:
    lines = f.readlines()

start_idx = None
end_idx = None

for i, line in enumerate(lines):
    if 'for trait, p_drug, gene, status, alt_drug, rationale in targeted:' in line:
        start_idx = i
    if start_idx is not None and i > start_idx and ('report_data =' in line or 'report_data' in line):
        end_idx = i
        break

if start_idx is not None and end_idx is not None:
    indent = "    "
    new_loop_lines = [
        f"{indent}for trait, p_drug, gene, status, alt_drug, rationale in targeted:\n",
        f"{indent}    prs_info = prs_results.get(trait, {{}})\n",
        f"{indent}    z_score = prs_info.get('z_score', 'N/A')\n",
        f"{indent}    percentile = prs_info.get('percentile', 'N/A')\n",
        f"{indent}    category = prs_info.get('category', 'N/A')\n",
        f"{indent}    \n",
        f"{indent}    print(f\"  ┌── [ DISEASE / TRAIT ]: {{trait.upper()}}\")\n",
        f"{indent}    print(f\"  │   ├── Polygenic Risk   : {{category}} (Z-Score: {{z_score}}, Percentile: {{percentile}}%)\")\n",
        f"{indent}    print(f\"  │   ├── Targeted Drug    : {{p_drug}}\")\n",
        f"{indent}    print(f\"  │   ├── Gene Evaluated   : {{gene}} (Patient PGx Status: {{status}})\")\n",
        f"{indent}    if 'CONTRAINDICATED' in str(status) or 'HIGH_RISK' in str(status) or 'REASSIGNED' in str(status):\n",
        f"{indent}        print(f\"  │   ├── Action Required  : ⚠️ SWITCH / REASSIGN -> {{alt_drug}}\")\n",
        f"{indent}    else:\n",
        f"{indent}        print(f\"  │   ├── Action Required  : ✔ APPROVED -> Maintain {{p_drug}}\")\n",
        f"{indent}    print(f\"  │   └── Clinical Rationale: {{rationale}}\")\n",
        f"{indent}    print(\"  └───────────────────────────────────────────────────────────────────\\n\")\n",
        "\n",
        f"{indent}print(\"=\" * 80)\n",
        f"{indent}print(\"[✔] Complete multi-engine report saved:\", args.output)\n",
        f"{indent}print(\"=\" * 80)\n",
        "\n"
    ]
    
    updated_lines = lines[:start_idx] + new_loop_lines + lines[end_idx:]
    
    with open('runner.py', 'w') as f:
        f.writelines(updated_lines)
    print("[✔] runner.py successfully repaired and formatted.")
else:
    print(f"[!] Could not locate target bounds (start: {start_idx}, end: {end_idx}).")
