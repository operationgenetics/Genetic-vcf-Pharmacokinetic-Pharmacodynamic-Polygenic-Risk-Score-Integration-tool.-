
import subprocess
from pathlib import Path

print("=== STARTING COMPREHENSIVE MULTI-DISORDER GENOMIC VALIDATION ===")

# 1. Generate an expansive VCF containing representative risk alleles for top-tier clinical conditions
vcf_path = Path("genomic_bot_workspace/exhaustive_all_disorders.vcf")
vcf_path.parent.mkdir(parents=True, exist_ok=True)

vcf_content = """##fileformat=VCFv4.2
##contig=<ID=1,length=248956422>
##contig=<ID=3,length=198295559>
##contig=<ID=6,length=170340590>
##contig=<ID=9,length=141213431>
##contig=<ID=10,length=133797422>
##contig=<ID=13,length=114364328>
##contig=<ID=17,length=83257441>
##contig=<ID=19,length=58617616>
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene symbol">
##INFO=<ID=TRAIT,Number=1,Type=String,Description="Associated trait">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
19	44908684	rs429358	T	C	100	PASS	GENE=APOE
10	94942290	rs1057910	A	C	100	PASS	GENE=CYP2C9
9	22125503	rs1333049	T	C	100	PASS	TRAIT=Cardiovascular_Disease
10	114754029	rs7903146	C	T	100	PASS	TRAIT=Type_2_Diabetes
3	51234500	rs2522122	A	G	100	PASS	TRAIT=Major_Depressive_Disease
6	31234567	rs9271101	A	G	100	PASS	TRAIT=Autoimmune_Disorder
1	150123456	rs1234567	G	A	100	PASS	TRAIT=Cardiometabolic_Syndrome
17	43044295	rs80357906	C	T	100	PASS	GENE=BRCA1
13	32900000	rs80359550	G	A	100	PASS	GENE=BRCA2
"""
vcf_path.write_text(vcf_content)
print("[✔] Comprehensive multi-disorder VCF generated.")

# 2. Execute the production pipeline bot
cmd = ["python3", "bot.py", "--vcf", str(vcf_path), "--output", "production_clinical_report.json"]
print(f"[ℹ] Running command: \" ".join(cmd))
result = subprocess.run(cmd, capture_output=True, text=True)

print(result.stdout)
if result.stderr:
    print("[STDERR]:", result.stderr)

print("=== VALIDATION COMPLETE ===")
