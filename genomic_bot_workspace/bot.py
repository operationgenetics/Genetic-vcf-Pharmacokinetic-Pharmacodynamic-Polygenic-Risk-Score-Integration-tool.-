#!/usr/bin/env python3
"""
bot.py - CLI tool for ultimate-genomic-bot.
"""

import argparse
import json
import sys
from pathlib import Path

from init_db import init_database
from precision_medicine_pipeline import PrecisionMedicinePipeline, DB_PATH


def main():
    parser = argparse.ArgumentParser(prog="genomic-bot", description="Automated Genomic & Precision Medicine CLI Tool")
    subparsers = parser.add_subparsers(dest="command", help="Available subcommands")

    subparsers.add_parser("initdb", help="Initialize local genomic SQLite database")

    run_parser = subparsers.add_parser("run", help="Run full pipeline analysis")
    run_parser.add_argument("--patient-id", required=True, help="Patient Unique Identifier")
    run_parser.add_argument("--vcf", required=True, help="Path to input VCF file")
    run_parser.add_argument("--meds", nargs="+", required=True, help="List of active medication names")
    run_parser.add_argument("--output", default="dashboard.json", help="Path for JSON output")

    sample_parser = subparsers.add_parser("samplevcf", help="Generate a mock patient VCF file")
    sample_parser.add_argument("--out", default="sample_patient.vcf", help="Output VCF path")

    args = parser.parse_args()

    if args.command == "initdb":
        init_database()
        print("[✔] Database initialized with comprehensive PRS, ClinVar, and PGx tables.")
    elif args.command == "samplevcf":
        out_path = Path(args.out)
        with open(out_path, "w") as f:
            f.write("##fileformat=VCFv4.2\n")
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            f.write("19\t44908684\trs429358\tT\tC\t100\tPASS\tGENE=APOE\n")
            f.write("10\t94942290\trs1057910\tA\tC\t100\tPASS\tGENE=CYP2C9\n")
            f.write("9\t22125503\trs1333049\tT\tC\t100\tPASS\tGENE=CDKN2B\n")
            f.write("10\t114754029\trs7903146\tC\tT\t100\tPASS\tGENE=TCF7L2\n")
        print(f"[✔] Comprehensive sample VCF created at '{out_path}'.")
    elif args.command == "run":
        if not DB_PATH.exists():
            print("[!] Knowledgebase missing. Generating database automatically...")
            init_database()

        pipeline = PrecisionMedicinePipeline()
        results = pipeline.execute_pipeline(
            patient_id=args.patient_id,
            raw_vcf_path=Path(args.vcf),
            active_meds_list=args.meds
        )

        with open(args.output, "w") as f:
            json.dump(results, f, indent=2)

        print(f"[✔] Pipeline complete! Comprehensive clinical dashboard written to '{args.output}'.")
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
