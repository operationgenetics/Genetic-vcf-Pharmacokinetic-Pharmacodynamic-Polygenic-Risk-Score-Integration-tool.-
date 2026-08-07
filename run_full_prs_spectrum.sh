#!/bin/bash
set -e

echo "================================================================================"
echo "      EXPANDING DB SCHEMA & RUNNING FULL-SPECTRUM PRS TEST (ALL TRAITS)        "
echo "================================================================================"

mkdir -p data results

# 1. Update/Seed Database with Comprehensive Psychiatric & Complex Trait PRS Metadata
python3 -c "
import sqlite3

conn = sqlite3.connect('genomics.db')
cursor = conn.cursor()

# Ensure tables exist
cursor.execute('''
CREATE TABLE IF NOT EXISTS prs_traits (
    trait_id INTEGER PRIMARY KEY AUTOINCREMENT,
    trait_name TEXT UNIQUE,
    mean_raw_score REAL,
    std_raw_score REAL
)''')

cursor.execute('''
CREATE TABLE IF NOT EXISTS prs_variants (
    variant_id INTEGER PRIMARY KEY AUTOINCREMENT,
    rsid TEXT,
    trait_name TEXT,
    effect_allele TEXT,
    weight REAL
)''')

# Extended list of PRS traits and baseline population parameters
traits = [
    ('Anxiety Disorder', 0.5, 0.3),
    ('Generalized Anxiety Disorder', 0.4, 0.25),
    ('Panic Disorder', 0.3, 0.2),
    ('Post-Traumatic Stress Disorder', 0.5, 0.3),
    ('Schizoaffective Disorder', 0.5, 0.4),
    ('Bipolar Disorder', 0.4, 0.35),
    ('Major Depressive Disorder', 0.6, 0.3),
    ('Coronary Artery Disease', 0.5, 0.3),
    ('Type 2 Diabetes', 0.5, 0.3)
]

for name, mean, std in traits:
    cursor.execute('''
        INSERT INTO prs_traits (trait_name, mean_raw_score, std_raw_score)
        VALUES (?, ?, ?)
        ON CONFLICT(trait_name) DO UPDATE SET mean_raw_score=excluded.mean_raw_score, std_raw_score=excluded.std_raw_score
    ''', (name, mean, std))

# Map high-risk variant weights for all traits
variants = [
    ('rs11178997', 'Anxiety Disorder', 'G', 0.8),
    ('rs28399433', 'Generalized Anxiety Disorder', 'T', 0.75),
    ('rs1024611',  'Schizoaffective Disorder', 'A', 0.85),
    ('rs9272219',  'Bipolar Disorder', 'T', 0.90),
    ('rs1799853',  'Panic Disorder', 'T', 0.70),
    ('rs4244285',  'Post-Traumatic Stress Disorder', 'A', 0.80),
    ('rs1065852',  'Major Depressive Disorder', 'T', 0.65)
]

for rsid, trait, allele, weight in variants:
    cursor.execute('''
        INSERT OR REPLACE INTO prs_variants (rsid, trait_name, effect_allele, weight)
        VALUES (?, ?, ?, ?)
    ''', (rsid, trait, allele, weight))

conn.commit()
conn.close()
print('[✔] SQLite Database successfully upgraded with full-spectrum PRS traits.')
"

# 2. Build synthetic VCF with risk alleles across all expanded PRS traits
UNCOMPRESSED_VCF="data/full_prs_spectrum.vcf"
COMPRESSED_VCF="data/full_prs_spectrum.vcf.gz"

echo "[1/4] Building full-spectrum PRS sample VCF..."
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

# 3. Compress and index VCF
echo "[2/4] Compressing and indexing VCF..."
if command -v bgzip &> /dev/null && command -v tabix &> /dev/null; then
    bgzip -cf $UNCOMPRESSED_VCF > $COMPRESSED_VCF
    tabix -f -p vcf $COMPRESSED_VCF
    rm -f $UNCOMPRESSED_VCF
else
    gzip -cf $UNCOMPRESSED_VCF > $COMPRESSED_VCF
    rm -f $UNCOMPRESSED_VCF
fi

# 4. Run Pytest Suite
echo "[3/4] Running test suite..."
pytest -v

# 5. Execute Runner
echo ""
echo "================================================================================"
echo "          RUNNING PIPELINE FOR ALL SPECTRUM PRS DISORDERS & TRAITS              "
echo "================================================================================"
python3 runner.py \
  --vcf $COMPRESSED_VCF \
  --patient-id FULL_SPECTRUM_PRS_01 \
  --output results/full_prs_spectrum_report.json

