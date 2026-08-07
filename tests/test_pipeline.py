import pytest
import sqlite3
import os
import json
import math
from runner import parse_vcf_qc, calculate_prs, evaluate_pgx, norm_cdf

DB_PATH = "genomic_knowledgebase.db"

@pytest.fixture(scope="module")
def db_conn():
    """Ensure database exists and yields connection."""
    assert os.path.exists(DB_PATH), "Database file missing. Run setup_db.py first."
    conn = sqlite3.connect(DB_PATH)
    yield conn
    conn.close()

def test_norm_cdf_values():
    """Verify standard normal CDF math."""
    assert round(norm_cdf(0.0), 2) == 0.50
    assert round(norm_cdf(1.96), 3) == 0.975
    assert round(norm_cdf(-1.96), 3) == 0.025

def test_vcf_qc_parser(tmp_path):
    """Test VCF QC filtering with pass and fail records."""
    vcf_content = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chr1\t100\trs1065852\tC\tT\t99\tPASS\t.\tGT:GQ:DP\t1/1:99:30
chr1\t200\trs4244285\tG\tA\t10\tLowQual\t.\tGT:GQ:DP\t0/1:10:5
"""
    vcf_file = tmp_path / "test.vcf"
    vcf_file.write_text(vcf_content)

    variants = parse_vcf_qc(str(vcf_file), min_gq=20, min_dp=10)
    
    # rs1065852 should pass QC
    assert "rs1065852" in variants
    assert variants["rs1065852"]["genotype"] == "T/T"
    
    # rs4244285 should be filtered out due to LowQual / low GQ
    assert "rs4244285" not in variants

def test_prs_zscore_calculation(db_conn):
    """Validate Z-score calculation against mock variant input."""
    mock_variants = {
        "rs9272219": {"alleles": ["T", "T"]}, # dosage = 2 -> 2 * 0.42 = 0.84
        "rs1024611": {"alleles": ["A", "A"]}, # dosage = 2 -> 2 * 0.38 = 0.76
        "rs1065852": {"alleles": ["T", "T"]}  # dosage = 2 -> 2 * 0.25 = 0.50
    }
    # Raw score total = 2.10
    prs = calculate_prs(mock_variants, db_conn)
    
    schiz = prs["Schizoaffective Disorder"]
    assert schiz["score"] == 2.10
    # mean = 1.20, std = 0.45 -> Z = (2.10 - 1.20) / 0.45 = 2.0
    assert schiz["z_score"] == 2.0
    assert schiz["percentile"] == 97.7
    assert schiz["category"] == "High Risk"

def test_pgx_star_allele_matching(db_conn):
    """Validate CYP2D6 *10/*10 diplotype resolution."""
    mock_variants = {
        "rs1065852": {"alleles": ["T", "T"]}
    }
    phenotypes, cpic, _, _, _ = evaluate_pgx(mock_variants, db_conn)
    
    assert "Poor Metabolizer (*10/*10)" in phenotypes["CYP2D6"]

