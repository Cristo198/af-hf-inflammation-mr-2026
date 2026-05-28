"""Run single-variant Wald-ratio MR for UKB-PPP inflammation instruments."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path
from statistics import NormalDist


PROJECT_ROOT = Path(__file__).resolve().parents[1]
EXPOSURE = PROJECT_ROOT / "data" / "processed" / "ukbppp_inflammation_cis_pqtl_instruments_lead_no_mhc.csv"
OUTCOME_DIR = PROJECT_ROOT / "data" / "outcomes"
RESULT_DIR = PROJECT_ROOT / "results" / "mr"

RESULT_FIELDS = [
    "protein",
    "protein_name",
    "panel",
    "outcome_id",
    "outcome_name",
    "SNP",
    "method",
    "beta_exposure",
    "se_exposure",
    "beta_outcome",
    "se_outcome",
    "beta_mr",
    "se_mr",
    "or",
    "or_lci95",
    "or_uci95",
    "pval_mr",
    "fdr",
    "harmonise_action",
    "effect_allele_exposure",
    "other_allele_exposure",
    "effect_allele_outcome",
    "other_allele_outcome",
]


COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def to_float(value: str) -> float:
    if value is None or value == "":
        return float("nan")
    return float(value)


def norm_allele(value: str) -> str:
    return (value or "").upper()


def complement(value: str) -> str:
    value = norm_allele(value)
    if len(value) != 1:
        return value
    return COMPLEMENT.get(value, value)


def harmonise(exp: dict[str, str], out: dict[str, str]) -> tuple[int, str]:
    ea_exp = norm_allele(exp["effect_allele"])
    oa_exp = norm_allele(exp["other_allele"])
    ea_out = norm_allele(out["effect_allele"])
    oa_out = norm_allele(out["other_allele"])

    if ea_exp == ea_out and oa_exp == oa_out:
        return 1, "aligned"
    if ea_exp == oa_out and oa_exp == ea_out:
        return -1, "flipped"
    if complement(ea_exp) == ea_out and complement(oa_exp) == oa_out:
        return 1, "aligned_complement"
    if complement(ea_exp) == oa_out and complement(oa_exp) == ea_out:
        return -1, "flipped_complement"
    return 0, "allele_mismatch"


def bh_fdr(rows: list[dict[str, str]]) -> None:
    valid = [(i, to_float(row["pval_mr"])) for i, row in enumerate(rows) if row["pval_mr"]]
    valid = [(i, p) for i, p in valid if p == p]
    n = len(valid)
    ranked = sorted(valid, key=lambda x: x[1])
    adjusted = [1.0] * n
    min_seen = 1.0
    for rank_idx in range(n - 1, -1, -1):
        original_index, p = ranked[rank_idx]
        rank = rank_idx + 1
        value = min(min_seen, p * n / rank)
        min_seen = value
        adjusted[rank_idx] = value
    for (original_index, _), fdr in zip(ranked, adjusted):
        rows[original_index]["fdr"] = f"{min(fdr, 1.0):.6g}"


def run_for_outcome(outcome_file: Path) -> list[dict[str, str]]:
    exposure_rows = {row["SNP"]: row for row in read_csv(EXPOSURE)}
    outcome_rows = read_csv(outcome_file)
    results = []
    normal = NormalDist()
    for out in outcome_rows:
        snp = out["SNP"]
        if snp not in exposure_rows:
            continue
        exp = exposure_rows[snp]
        flip, harmonise_action = harmonise(exp, out)
        if flip == 0:
            continue
        beta_exp = to_float(exp["beta"])
        se_exp = to_float(exp["se"])
        beta_out = to_float(out["beta"]) * flip
        se_out = to_float(out["se"])
        if not all(value == value for value in [beta_exp, se_exp, beta_out, se_out]) or beta_exp == 0:
            continue

        beta_mr = beta_out / beta_exp
        se_mr = abs(se_out / beta_exp)
        z = beta_mr / se_mr if se_mr else float("nan")
        pval = 2 * (1 - normal.cdf(abs(z))) if z == z else float("nan")
        ci_low = beta_mr - 1.96 * se_mr
        ci_high = beta_mr + 1.96 * se_mr
        results.append(
            {
                "protein": exp["protein"],
                "protein_name": exp["protein_name"],
                "panel": exp["panel"],
                "outcome_id": out["outcome_id"],
                "outcome_name": out["outcome_name"],
                "SNP": snp,
                "method": "Wald ratio",
                "beta_exposure": exp["beta"],
                "se_exposure": exp["se"],
                "beta_outcome": f"{beta_out:.8g}",
                "se_outcome": out["se"],
                "beta_mr": f"{beta_mr:.8g}",
                "se_mr": f"{se_mr:.8g}",
                "or": f"{math.exp(beta_mr):.8g}",
                "or_lci95": f"{math.exp(ci_low):.8g}",
                "or_uci95": f"{math.exp(ci_high):.8g}",
                "pval_mr": f"{pval:.8g}" if pval == pval else "",
                "fdr": "",
                "harmonise_action": harmonise_action,
                "effect_allele_exposure": exp["effect_allele"],
                "other_allele_exposure": exp["other_allele"],
                "effect_allele_outcome": out["effect_allele"],
                "other_allele_outcome": out["other_allele"],
            }
        )
    bh_fdr(results)
    return results


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outcome-files", nargs="*", default=[])
    args = parser.parse_args()

    outcome_files = [Path(path) for path in args.outcome_files]
    if not outcome_files:
        outcome_files = sorted(OUTCOME_DIR.glob("*_outcome_effects_*.csv"))
    if not outcome_files:
        print("No outcome effect files found in data/outcomes. Run script 04 first.")
        sys.exit(0)

    all_results = []
    for outcome_file in outcome_files:
        rows = run_for_outcome(outcome_file)
        all_results.extend(rows)
        out_path = RESULT_DIR / f"{outcome_file.stem.replace('_outcome_effects_local','').replace('_outcome_effects_opengwas','')}_wald_mr.csv"
        write_csv(out_path, rows, RESULT_FIELDS)
        print(f"{outcome_file.name}: MR rows {len(rows)} -> {out_path}")

    if all_results:
        write_csv(RESULT_DIR / "wald_mr_all_outcomes.csv", all_results, RESULT_FIELDS)
        print(f"Combined MR rows {len(all_results)} -> {RESULT_DIR / 'wald_mr_all_outcomes.csv'}")


if __name__ == "__main__":
    main()
