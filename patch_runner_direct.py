with open("runner.py", "r") as f:
    lines = f.readlines()

start_line = None
end_line = None

for idx, line in enumerate(lines):
    if "for trait, p_drug, gene, status, alt_drug, rationale in targeted:" in line:
        start_line = idx
        break

if start_line is not None:
    # Find the next block that starts at outer indentation or opens/writes a file
    for idx in range(start_line + 1, len(lines)):
        line_str = lines[idx]
        if "report_data =" in line_str or "with open" in line_str or "json.dump" in line_str or "print(\"=\" * 80)" in line_str:
            end_line = idx
            break

if start_line is not None and end_line is not None:
    indent = "    "
    new_block = [
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
        f"{indent}    if any(k in str(status).upper() for k in ['CONTRAINDICATED', 'HIGH_RISK', 'REASSIGNED']):\n",
        f"{indent}        print(f\"  │   ├── Action Required  : ⚠️ SWITCH / REASSIGN -> {{alt_drug}}\")\n",
        f"{indent}    else:\n",
        f"{indent}        print(f\"  │   ├── Action Required  : ✔ APPROVED -> Maintain {{p_drug}}\")\n",
        f"{indent}    print(f\"  │   └── Clinical Rationale: {{rationale}}\")\n",
        f"{indent}    print(\"  └───────────────────────────────────────────────────────────────────\\n\")\n",
        "\n",
    ]
    
    updated_lines = lines[:start_line] + new_block + lines[end_line:]
    
    with open("runner.py", "w") as f:
        f.writelines(updated_lines)
    print(f"[✔] Successfully patched runner.py between lines {start_line} and {end_line}.")
else:
    print(f"[!] Target resolution failed: start_line={start_line}, end_line={end_line}")
