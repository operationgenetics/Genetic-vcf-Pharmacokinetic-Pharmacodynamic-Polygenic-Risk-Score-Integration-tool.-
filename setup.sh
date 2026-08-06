#!/usr/bin/env bash
set -e

echo "🚀 Setting up Genetic-vcf pipeline workspace..."
mkdir -p data results

if [ ! -f "genomic_knowledgebase.db" ]; then
    python3 setup_db.py
fi

if [ ! -f "data/pgp_na12878.vcf.gz" ]; then
    echo "📥 Downloading sample PGP benchmark VCF (NA12878)..."
    wget -q -O data/pgp_na12878.vcf.gz https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
fi

echo "✅ Setup complete! You can now run:"
echo "   python runner.py --vcf data/pgp_na12878.vcf.gz --patient-id Sample01"
