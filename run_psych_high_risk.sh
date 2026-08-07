#!/bin/bash
set -e

echo "================================================================================"
echo "    GENERATING HIGH-RISK PSYCHIATRIC (ANXIETY, SCHIZO, BIPOLAR) VCF.GZ & TBI    "
echo "================================================================================"

# 1. Prepare directories
mkdir -p data results

UNCOMPRESSED_VCF="data/psych_high_risk_triple.vcf"
COMPRESSED_VCF="data/psych_high_risk_triple.vcf.gz"

# 2. Build synthetic VCF with homozygous high-risk effect alleles for Anxiety, Schizoaffective, and Bipolar
echo "[1/5] Building sample VCF with high-risk PRS variant alleles..."
cat << 'VCF' > $UNCOMPRESSED_VCF
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr1	1000	rs1065852	C	T	99	PASS	.	GT:GQ:DP	1/1:99:45
chr2	2000	rs4244285	G	A	99	PASS	.	GT:GQ:DP	1/1:99:50
chr3	3000	rs1024611	C	A	99	PASS	.	GT:GQ:DP	1/1:99:55
chr4	4000	rs9272219	C	T	99	PASS	.	GT:GQ:DP	1/1:99:60
chr5	5000	rs1799853	C	T	99	PASS	.	GT:GQ:DP	1/1:99:40
chr6	6000	rs28399433	C	T	99	PASS	.	GT:GQ:DP	1/1:99:42
chr7	7000	rs11178997	A	G	99	PASS	.	GT:GQ:DP	1/1:99:50
VCF

# 3. Compress VCF with gzip
echo "[2/5] Compressing VCF to .vcf.gz format..."
gzip -c $UNCOMPRESSED_VCF > $COMPRESSED_VCF
rm -f $UNCOMPRESSED_VCF

# 4. Generate tabix index (.tbi) if bgzip/tabix tools are present
echo "[3/5] Indexing VCF file..."
if command -v tabix &> /dev/null && command -v bgzip &> /dev/null; then
    gunzip -c $COMPRESSED_VCF | bgzip -c > "${COMPRESSED_VCF}.tmp"
    mv "${COMPRESSED_VCF}.tmp" $COMPRESSED_VCF
    tabix -p vcf $COMPRESSED_VCF
    echo "[✔] Tabix index generated successfully: ${COMPRESSED_VCF}.tbi"
else
    echo "[i] Proceeding with Python stream parsing."
fi

# 5. Initialize Database Schema
echo "[4/5] Initializing genomic database..."
python3 setup_db.py

# 6. Execute Unit Tests
echo "[5/5] Executing Pytest suite..."
pytest -v

# 7. Execute Runner
echo ""
echo "================================================================================"
echo "          RUNNING PIPELINE ON HIGH-RISK PSYCHIATRIC SAMPLE VCF                  "
echo "================================================================================"
python3 runner.py \
  --vcf $COMPRESSED_VCF \
  --patient-id PSYCH_HIGH_RISK_TRIPLE_01 \
  --output results/psych_high_risk_report.json

