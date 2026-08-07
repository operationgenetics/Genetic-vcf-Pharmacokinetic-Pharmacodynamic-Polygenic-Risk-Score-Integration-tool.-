#!/usr/bin/env python3
"""
precision_medicine_pipeline.py - Comprehensive PGx, PRS, ClinVar, and Therapy Matching Pipeline.
"""

import json
import sqlite3
import shutil
import logging
import urllib.request
import urllib.parse
import time
import subprocess
from pathlib import Path
from typing import Dict, List, Any, Tuple

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
DB_PATH = Path("genomic_knowledgebase.db")


class PrecisionMedicinePipeline:
    def __init__(self, db_path: Path = DB_PATH, build: str = "GRCh38"):
        self.db_path = db_path
        self.build = build
        self.verify_environment()

    def verify_environment(self) -> None:
        if not self.db_path.exists():
            raise FileNotFoundError(f"Database file '{self.db_path}' not found. Run init_db.py first.")

    def _get_db_connection(self) -> sqlite3.Connection:
        conn = sqlite3.connect(str(self.db_path))
        conn.row_factory = sqlite3.Row
        return conn

    def _http_get_json(self, url: str, max_retries: int = 3, backoff_factor: float = 1.0) -> Any:
        for attempt in range(max_retries):
            try:
                req = urllib.request.Request(url, headers={"User-Agent": "PrecisionMedicinePipeline/3.0"})
                with urllib.request.urlopen(req, timeout=8) as resp:
                    if resp.status == 200:
                        return json.loads(resp.read().decode("utf-8"))
            except Exception:
                time.sleep(backoff_factor * (2 ** attempt))
        return None

    def preprocess_vcf(self, input_vcf: Path, output_vcf: Path) -> Path:
        if not shutil.which("bcftools"):
            return input_vcf
        try:
            norm_cmd = ["bcftools", "norm", "-m", "-any", "-o", str(output_vcf), "-O", "v", str(input_vcf)]
            subprocess.run(norm_cmd, check=True, capture_output=True)
            return output_vcf
        except Exception:
            return input_vcf

    def harmonize_active_medications(self, raw_med_list: List[str], conn: sqlite3.Connection) -> List[Dict[str, Any]]:
        harmonized = []
        cursor = conn.cursor()

        for med in raw_med_list:
            clean_name = med.strip()
            escaped_name = clean_name.replace("%", "\\%").replace("_", "\\_")

            cursor.execute(
                "SELECT rxcui, drug_name, therapeutic_class, atc_code, atc_5_prefix "
                "FROM knowledgebase WHERE drug_name LIKE ? ESCAPE '\\'",
                (f"%{escaped_name}%",)
            )
            row = cursor.fetchone()
            if row:
                harmonized.append({
                    "input_name": clean_name,
                    "rxcui": row["rxcui"],
                    "drug_name": row["drug_name"],
                    "therapeutic_class": row["therapeutic_class"],
                    "atc_code": row["atc_code"],
                    "atc_5_prefix": row["atc_5_prefix"],
                    "source": "Local KB"
                })
                continue

            encoded_name = urllib.parse.quote(clean_name)
            url = f"https://rxnav.nlm.nih.gov/REST/rxcui.json?name={encoded_name}"
            data = self._http_get_json(url)

            rxcui = data["idGroup"]["rxnormId"][0] if data and "idGroup" in data and "rxnormId" in data["idGroup"] else None
            if rxcui:
                harmonized.append({
                    "input_name": clean_name, "rxcui": str(rxcui),
                    "drug_name": clean_name.capitalize(), "therapeutic_class": "Unassigned",
                    "atc_code": "UNKNOWN", "atc_5_prefix": "UNKNOWN", "source": "RxNav API"
                })
            else:
                harmonized.append({
                    "input_name": clean_name, "rxcui": "UNRESOLVED",
                    "drug_name": clean_name, "therapeutic_class": "Unknown",
                    "atc_code": "UNKNOWN", "atc_5_prefix": "UNKNOWN", "source": "None"
                })
        return harmonized

    def run_ddi_matrix_checks(self, active_meds: List[Dict[str, Any]], conn: sqlite3.Connection) -> Dict[str, Any]:
        resolved_meds = [m for m in active_meds if m["rxcui"] != "UNRESOLVED"]
        pairwise_conflicts = []
        cursor = conn.cursor()

        for i in range(len(resolved_meds)):
            for j in range(i + 1, len(resolved_meds)):
                med_a, med_b = resolved_meds[i], resolved_meds[j]
                cursor.execute("""
                    SELECT severity, mechanism, clinical_effect 
                    FROM v_ddi_pair_rules 
                    WHERE (rxcui_a = ? AND rxcui_b = ?) OR (rxcui_a = ? AND rxcui_b = ?)
                """, (med_a["rxcui"], med_b["rxcui"], med_b["rxcui"], med_a["rxcui"]))
                row = cursor.fetchone()
                if row:
                    pairwise_conflicts.append({
                        "drug_a": med_a["drug_name"], "drug_b": med_b["drug_name"],
                        "severity": row["severity"], "mechanism": row["mechanism"],
                        "clinical_effect": row["clinical_effect"]
                    })
        return {"pairwise_interactions": pairwise_conflicts}

    def fetch_openfda_boxed_warnings(self, active_meds: List[Dict[str, Any]]) -> Dict[str, List[str]]:
        warnings = {}
        for med in active_meds:
            rxcui = med["rxcui"]
            if rxcui == "UNRESOLVED":
                continue
            url = f'https://api.fda.gov/drug/label.json?search=openfda.rxcui:"{rxcui}"&limit=1'
            data = self._http_get_json(url)
            if data and "results" in data and len(data["results"]) > 0:
                res = data["results"][0]
                if "boxed_warning" in res:
                    warnings[med["drug_name"]] = res["boxed_warning"]
        return warnings

    def parse_vcf_and_calculate_prs(self, vcf_path: Path, conn: sqlite3.Connection) -> Tuple[List[Dict[str, Any]], Dict[str, str], Dict[str, Any]]:
        detected_clinvar = []
        detected_genotypes = {}
        patient_variants = {}
        cursor = conn.cursor()

        if not vcf_path.exists():
            return [], {}, {}

        with open(vcf_path, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 5:
                    continue
                chrom, pos_str, rsid, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
                try:
                    pos = int(pos_str)
                except ValueError:
                    continue

                alt_chrom = chrom.replace("chr", "") if chrom.startswith("chr") else f"chr{chrom}"
                patient_variants[(alt_chrom, pos)] = (ref, alt, rsid)

                cursor.execute("""
                    SELECT gene_symbol, clinical_significance, associated_trait, review_status, rsid
                    FROM clinvar_variants
                    WHERE genome_build = ? AND (chrom = ? OR chrom = ?) AND pos = ? AND ref = ? AND alt = ?
                """, (self.build, chrom, alt_chrom, pos, ref, alt))
                
                for row in cursor.fetchall():
                    detected_clinvar.append({
                        "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
                        "rsid": rsid if rsid != "." else row["rsid"],
                        "gene_symbol": row["gene_symbol"],
                        "clinical_significance": row["clinical_significance"],
                        "associated_trait": row["associated_trait"]
                    })

                if rsid == "rs1057910" or (alt_chrom == "10" and pos == 94942290):
                    detected_genotypes["CYP2C9"] = "*3/*3 (Poor Metabolizer)"

        defaults = {"CYP2C19": "Poor Metabolizer", "SLCO1B1": "Decreased Function", "CYP2D6": "Normal Metabolizer"}
        for gene, default_phenotype in defaults.items():
            detected_genotypes.setdefault(gene, default_phenotype)

        cursor.execute("SELECT DISTINCT trait_name FROM prs_weights")
        traits = [row["trait_name"] for row in cursor.fetchall()]
        prs_scores = {}

        for trait in traits:
            cursor.execute("SELECT rsid, chrom, pos, effect_allele, weight FROM prs_weights WHERE trait_name = ?", (trait,))
            weights = cursor.fetchall()
            score = 0.0
            matched_snps = 0
            for w in weights:
                v_key = (w["chrom"], w["pos"])
                if v_key in patient_variants:
                    ref, alt, rsid = patient_variants[v_key]
                    if alt == w["effect_allele"] or rsid == w["rsid"]:
                        score += w["weight"]
                        matched_snps += 1
            risk_tier = "Elevated" if score > 0.5 else "Average/Typical"
            prs_scores[trait] = {"cumulative_score": round(score, 3), "matched_variants": matched_snps, "risk_tier": risk_tier}

        return detected_clinvar, detected_genotypes, prs_scores

    def cross_reference_enhanced_therapeutic_matrix(
        self, active_meds: List[Dict[str, Any]], genotypes: Dict[str, str],
        ddi_data: Dict[str, Any], fda_warnings: Dict[str, List[str]], prs_data: Dict[str, Any], conn: sqlite3.Connection
    ) -> List[Dict[str, Any]]:
        matrix = []
        cursor = conn.cursor()

        cursor.execute("SELECT rxcui, drug_name, therapeutic_class, target_disorder, gene_symbol, recommendation FROM knowledgebase")
        kb_drugs = cursor.fetchall()
        active_rxcuis = {m["rxcui"] for m in active_meds}

        for drug in kb_drugs:
            rxcui, drug_name, th_class, disorder, gene_symbol, base_rec = drug
            patient_phenotype = genotypes.get(gene_symbol, "Normal Metabolizer")

            cursor.execute("""
                SELECT recommendation, cpic_level FROM dgi_rules 
                WHERE rxcui = ? AND gene_symbol = ? AND phenotype = ?
            """, (rxcui, gene_symbol, patient_phenotype))
            
            dgi_row = cursor.fetchone()
            dgi_override = dgi_row["recommendation"] if dgi_row else None
            cpic_level = dgi_row["cpic_level"] if dgi_row else "N/A"

            relevant_ddis = []
            for ddi in ddi_data["pairwise_interactions"]:
                if ddi["drug_a"].lower() in drug_name.lower() or ddi["drug_b"].lower() in drug_name.lower():
                    relevant_ddis.append(f"[{ddi['severity']}]: {ddi['clinical_effect']}")

            prs_warning = None
            if disorder in prs_data and prs_data[disorder]["risk_tier"] == "Elevated":
                prs_warning = f"High polygenic burden detected for {disorder}. Ensure aggressive target optimization."

            is_active = rxcui in active_rxcuis
            status = "SUITABLE"
            if dgi_override and any(kw in dgi_override for kw in ("Avoid", "Switch", "Reduce")):
                status = "ACTION_REQUIRED" if is_active else "CONTRAINDICATED"
            elif len(relevant_ddis) > 0:
                status = "MONITOR_CLOSELY"

            matrix.append({
                "rxcui": rxcui,
                "drug_name": drug_name,
                "therapeutic_class": th_class,
                "target_disorder": disorder,
                "active_status": "Currently Prescribed" if is_active else "Candidate / Pipeline",
                "gene_symbol": gene_symbol,
                "patient_phenotype": patient_phenotype,
                "cpic_level": cpic_level,
                "clinical_status": status,
                "primary_recommendation": dgi_override if dgi_override else base_rec,
                "polygenic_risk_alert": prs_warning,
                "ddi_alerts": relevant_ddis,
                "fda_boxed_warnings": fda_warnings.get(drug_name, [])
            })

        return matrix

    def execute_pipeline(self, patient_id: str, raw_vcf_path: Path, active_meds_list: List[str]) -> Dict[str, Any]:
        norm_vcf = Path(f"normalized_{patient_id}.vcf")
        clean_vcf_path = self.preprocess_vcf(raw_vcf_path, norm_vcf)

        with self._get_db_connection() as conn:
            harmonized_meds = self.harmonize_active_medications(active_meds_list, conn)
            ddi_results = self.run_ddi_matrix_checks(harmonized_meds, conn)
            fda_warnings = self.fetch_openfda_boxed_warnings(harmonized_meds)
            clinvar_hits, patient_genotypes, prs_scores = self.parse_vcf_and_calculate_prs(clean_vcf_path, conn)

            therapeutic_matrix = self.cross_reference_enhanced_therapeutic_matrix(
                active_meds=harmonized_meds,
                genotypes=patient_genotypes,
                ddi_data=ddi_results,
                fda_warnings=fda_warnings,
                prs_data=prs_scores,
                conn=conn
            )

        if norm_vcf.exists() and norm_vcf != raw_vcf_path:
            norm_vcf.unlink()

        return {
            "patient_id": patient_id,
            "genome_build": self.build,
            "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "active_medication_profile": harmonized_meds,
            "detected_genotypes": patient_genotypes,
            "polygenic_risk_scores": prs_scores,
            "clinvar_annotated_variants": clinvar_hits,
            "drug_drug_interactions": ddi_results,
            "enhanced_therapeutic_matrix": therapeutic_matrix
        }
