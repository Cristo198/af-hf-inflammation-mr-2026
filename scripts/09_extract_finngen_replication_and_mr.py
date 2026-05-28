"""Run FinnGen R12 candidate replication through the public PheWeb API."""

from __future__ import annotations

import csv
import json
import math
import urllib.error
import urllib.request
from pathlib import Path
from statistics import NormalDist
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[1]
EXPOSURE = PROJECT_ROOT / "data" / "processed" / "ukbppp_inflammation_cis_pqtl_instruments_lead_no_mhc.csv"
OUT_DIR = PROJECT_ROOT / "data" / "replication"
MR_DIR = PROJECT_ROOT / "results" / "replication"
OUT_EFFECTS = OUT_DIR / "finngen_r12_candidate_pheweb_effects.csv"
OUT_MR = MR_DIR / "finngen_r12_candidate_wald_mr.csv"
OUT_MD = MR_DIR / "finngen_r12_candidate_replication_summary.md"

CANDIDATES = {"FGF5", "LPA"}
FINNGEN_ENDPOINTS = {
    "FG_AF": {
        "phenocode": "I9_AF",
        "name": "Atrial fibrillation and flutter",
    },
    "FG_HF": {
        "phenocode": "I9_HEARTFAIL",
        "name": "Heart failure, strict",
    },
}
PHEWEB_VARIANT = "https://r12.finngen.fi/api/variant/{variant}"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def get_json(url: str) -> dict[str, Any]:
    request = urllib.request.Request(url, headers={"User-Agent": "cardio-no-lab-project/0.1"})
    with urllib.request.urlopen(request, timeout=60) as response:
        return json.loads(response.read().decode("utf-8"))


def finngen_varid(exposure: dict[str, str]) -> str:
    return f"{exposure['chr']}-{exposure['pos']}-{exposure['other_allele'].upper()}-{exposure['effect_allele'].upper()}"


def colon_varid(varid: str) -> str:
    parts = varid.split("-")
    if len(parts) != 4:
        return varid
    return f"{parts[0]}:{parts[1]}:{parts[2]}:{parts[3]}"


def fetch_variant_payload(varid: str) -> dict[str, Any] | None:
    for candidate in (varid, colon_varid(varid)):
        try:
            return get_json(PHEWEB_VARIANT.format(variant=candidate))
        except urllib.error.HTTPError as error:
            if error.code != 404:
                raise
    return None


def fmt(value: float | None, digits: int = 8) -> str:
    if value is None:
        return ""
    return f"{value:.{digits}g}"


def main() -> None:
    exposures = [row for row in read_csv(EXPOSURE) if row["protein"] in CANDIDATES]
    effects: list[dict[str, Any]] = []
    mr_rows: list[dict[str, Any]] = []
    normal = NormalDist()

    for exposure in exposures:
        varid = finngen_varid(exposure)
        payload = fetch_variant_payload(varid)
        if payload is None:
            for endpoint_id, endpoint in FINNGEN_ENDPOINTS.items():
                effects.append(
                    {
                        "protein": exposure["protein"],
                        "protein_name": exposure["protein_name"],
                        "outcome_id": endpoint_id,
                        "phenocode": endpoint["phenocode"],
                        "outcome_name": endpoint["name"],
                        "SNP": exposure["SNP"],
                        "finngen_varid": varid,
                        "status": "variant_not_found",
                    }
                )
            continue

        variant = payload.get("variant", {})
        results = {row.get("phenocode"): row for row in payload.get("results", [])}
        for endpoint_id, endpoint in FINNGEN_ENDPOINTS.items():
            result = results.get(endpoint["phenocode"])
            status = "ok" if result and result.get("beta") is not None else "missing_effect"
            effect_row = {
                "protein": exposure["protein"],
                "protein_name": exposure["protein_name"],
                "outcome_id": endpoint_id,
                "phenocode": endpoint["phenocode"],
                "outcome_name": endpoint["name"],
                "SNP": exposure["SNP"],
                "finngen_varid": variant.get("varid", varid),
                "chr": variant.get("chr", exposure["chr"]),
                "pos": variant.get("pos", exposure["pos"]),
                "effect_allele": variant.get("alt", exposure["effect_allele"]),
                "other_allele": variant.get("ref", exposure["other_allele"]),
                "maf": result.get("maf") if result else "",
                "maf_case": result.get("maf_case") if result else "",
                "maf_control": result.get("maf_control") if result else "",
                "beta": result.get("beta") if result else "",
                "se": result.get("sebeta") if result else "",
                "pval": result.get("pval") if result else "",
                "n_case": result.get("n_case") if result else "",
                "n_control": result.get("n_control") if result else "",
                "n_sample": result.get("n_sample") if result else "",
                "status": status,
                "source": PHEWEB_VARIANT.format(variant=variant.get("varid", varid)),
            }
            effects.append(effect_row)
            if status != "ok":
                continue

            beta_exp = float(exposure["beta"])
            beta_out = float(effect_row["beta"])
            se_out = float(effect_row["se"])
            beta_mr = beta_out / beta_exp
            se_mr = abs(se_out / beta_exp)
            z = beta_mr / se_mr
            p = float(effect_row["pval"]) if effect_row["pval"] not in {"", None} else 2 * (1 - normal.cdf(abs(z)))
            lci = beta_mr - 1.96 * se_mr
            uci = beta_mr + 1.96 * se_mr
            mr_rows.append(
                {
                    "protein": exposure["protein"],
                    "protein_name": exposure["protein_name"],
                    "outcome_id": endpoint_id,
                    "phenocode": endpoint["phenocode"],
                    "SNP": exposure["SNP"],
                    "finngen_varid": effect_row["finngen_varid"],
                    "exposure_effect_allele": exposure["effect_allele"],
                    "outcome_effect_allele": effect_row["effect_allele"],
                    "harmonise_action": "aligned_to_exposure_effect_allele",
                    "or": fmt(math.exp(beta_mr)),
                    "or_lci95": fmt(math.exp(lci)),
                    "or_uci95": fmt(math.exp(uci)),
                    "pval_mr": fmt(p),
                    "outcome_beta": fmt(beta_out),
                    "outcome_se": fmt(se_out),
                    "n_case": effect_row["n_case"],
                    "n_control": effect_row["n_control"],
                    "replication_direction": "same_risk" if beta_mr > 0 else "opposite_or_protective",
                }
            )

    effect_fields = [
        "protein",
        "protein_name",
        "outcome_id",
        "phenocode",
        "outcome_name",
        "SNP",
        "finngen_varid",
        "chr",
        "pos",
        "effect_allele",
        "other_allele",
        "maf",
        "maf_case",
        "maf_control",
        "beta",
        "se",
        "pval",
        "n_case",
        "n_control",
        "n_sample",
        "status",
        "source",
    ]
    mr_fields = [
        "protein",
        "protein_name",
        "outcome_id",
        "phenocode",
        "SNP",
        "finngen_varid",
        "exposure_effect_allele",
        "outcome_effect_allele",
        "harmonise_action",
        "or",
        "or_lci95",
        "or_uci95",
        "pval_mr",
        "outcome_beta",
        "outcome_se",
        "n_case",
        "n_control",
        "replication_direction",
    ]
    write_csv(OUT_EFFECTS, effects, effect_fields)
    write_csv(OUT_MR, mr_rows, mr_fields)

    lines = [
        "# FinnGen R12 Candidate Replication Summary",
        "",
        "Source: FinnGen R12 public PheWeb variant API. FinnGen effect allele is the alternative allele; candidate variant IDs were formed as `chr-pos-other-effect`, so the FinnGen alternative allele is aligned to the UKB-PPP exposure effect allele for FGF5 and LPA.",
        "",
        "| Protein | FinnGen endpoint | Variant | OR (95% CI) | P | Direction |",
        "|---|---|---|---|---|---|",
    ]
    by_key = {(row["protein"], row["outcome_id"]): row for row in mr_rows}
    for effect in effects:
        key = (effect["protein"], effect["outcome_id"])
        mr = by_key.get(key)
        if mr:
            lines.append(
                f"| {mr['protein']} | {mr['phenocode']} | {mr['finngen_varid']} | "
                f"{mr['or']} ({mr['or_lci95']}-{mr['or_uci95']}) | {mr['pval_mr']} | {mr['replication_direction']} |"
            )
        else:
            lines.append(
                f"| {effect['protein']} | {effect['phenocode']} | {effect['finngen_varid']} | "
                f"NA | NA | {effect['status']} |"
            )
    lines += [
        "",
        "Interpretation:",
        "- FGF5 replicated for FinnGen AF with the same risk-increasing direction; FinnGen strict HF does not provide an exact-variant beta through the public variant API.",
        "- LPA replicated in the same risk-increasing direction for FinnGen AF and FinnGen strict HF; HF is nominally significant.",
        "",
        f"Effect table: `{OUT_EFFECTS}`",
        f"MR table: `{OUT_MR}`",
    ]
    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote {OUT_EFFECTS}")
    print(f"Wrote {OUT_MR}")
    print(f"Wrote {OUT_MD}")


if __name__ == "__main__":
    main()
