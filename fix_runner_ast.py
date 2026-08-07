import ast

with open("runner.py", "r") as f:
    lines = f.readlines()

# Locate line index for the target loop and the output block
start_line = None
end_line = None

for idx, line in enumerate(lines):
    if "for trait, p_drug, gene, status, alt_drug, rationale in targeted:" in line:
        start_line = idx
    if start_line is not None and idx > start_line and ("with open(args.output" in line or "with open(args" in line):
        end_line = idx
        break

if start_line is not None and end_line is not None:
    indent = "    "
    replacement_code = [
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
        f"{indent}report_data = {{\n",
        f"{indent}    'patient_id': args.patient_id,\n",
        f"{indent}    'vcf_file': args.vcf,\n",
        f"{indent}    'pgx_phenotypes': pgx_phenotypes if 'pgx_phenotypes' in locals() else {{}},\n",
        f"{indent}    'prs_scores': prs_results if 'prs_results' in locals() else {{}},\n",
        f"{indent}    'targeted_therapies': [\n",
        f"{indent}        {{\n",
        f"{indent}            'trait': trait,\n",
        f"{indent}            'polygenic_risk': prs_results.get(trait, {{}}),\n",
        f"{indent}            'primary_drug': p_drug,\n",
        f"{indent}            'gene': gene,\n",
        f"{indent}            'pgx_status': status,\n",
        f"{indent}            'alternative_drug': alt_drug,\n",
        f"{indent}            'rationale': rationale\n",
        f"{indent}        }}\n",
        f"{indent}        for trait, p_drug, gene, status, alt_drug, rationale in targeted\n",
        f"{indent}    ]\n",
        f"{indent}}}\n",
        "\n"
    ]

    new_lines = lines[:start_line] + replacement_code + lines[end_line:]
    
    with open("runner.py", "w") as f:
        f.writelines(new_lines)
    print("[✔] runner.py successfully re-built with AST line indices.")
else:
    print(f"[!] Unable to match markers (start_line: {start_line}, end_line: {end_line})")
