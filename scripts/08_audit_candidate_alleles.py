"""Audit allele direction for preliminary shared candidates."""

from __future__ import annotations

import csv
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
EXPOSURE = PROJECT_ROOT / "data" / "processed" / "ukbppp_inflammation_cis_pqtl_instruments_lead_no_mhc.csv"
MR_ALL = PROJECT_ROOT / "results" / "mr" / "wald_mr_all_outcomes.csv"
OUTCOME_DIR = PROJECT_ROOT / "data" / "outcomes"
OUT_CSV = PROJECT_ROOT / "results" / "mr" / "candidate_allele_audit.csv"
OUT_MD = PROJECT_ROOT / "results" / "mr" / "candidate_allele_audit.md"
CANDIDATES = {"FGF5", "LPA"}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def flt(value: str) -> float | None:
    try:
        return float(value)
    except Exception:
        return None


def aligned_outcome_eaf(exposure: dict[str, str], outcome: dict[str, str]) -> str:
    exp_ea = exposure["effect_allele"].upper()
    out_ea = outcome["effect_allele"].upper()
    out_eaf = flt(outcome["eaf"])
    if out_eaf is None:
        return ""
    if exp_ea == out_ea:
        return f"{out_eaf:.6g}"
    return f"{1 - out_eaf:.6g}"


def main() -> None:
    exposure_by_protein = {row["protein"]: row for row in read_csv(EXPOSURE)}
    mr_rows = [row for row in read_csv(MR_ALL) if row["protein"] in CANDIDATES]
    outcome_rows = {}
    for path in OUTCOME_DIR.glob("*_outcome_effects_local.csv"):
        outcome_id = path.name.split("_")[0]
        outcome_rows[outcome_id] = {row["SNP"]: row for row in read_csv(path)}

    rows = []
    for mr in mr_rows:
        exp = exposure_by_protein[mr["protein"]]
        out = outcome_rows[mr["outcome_id"]][mr["SNP"]]
        rows.append(
            {
                "protein": mr["protein"],
                "protein_name": mr["protein_name"],
                "outcome_id": mr["outcome_id"],
                "SNP": mr["SNP"],
                "exposure_effect_allele": exp["effect_allele"],
                "exposure_other_allele": exp["other_allele"],
                "exposure_eaf": exp["eaf"],
                "exposure_beta": exp["beta"],
                "outcome_effect_allele": out["effect_allele"].upper(),
                "outcome_other_allele": out["other_allele"].upper(),
                "outcome_eaf_raw": out["eaf"],
                "outcome_eaf_aligned_to_exposure_effect": aligned_outcome_eaf(exp, out),
                "outcome_beta_raw": out["beta"],
                "harmonise_action": mr["harmonise_action"],
                "outcome_beta_aligned": mr["beta_outcome"],
                "wald_or": mr["or"],
                "wald_ci": f"{mr['or_lci95']}-{mr['or_uci95']}",
                "wald_p": mr["pval_mr"],
                "wald_fdr": mr["fdr"],
            }
        )

    fieldnames = list(rows[0].keys())
    write_csv(OUT_CSV, rows, fieldnames)

    lines = [
        "# Candidate Allele Direction Audit",
        "",
        "Important correction: UKB-PPP browser JSON fields are named `EA` and `OA`, but the browser display labels `OA` as the effect allele and `EA` as the other allele. EAF matches the `OA` frequency for audited variants. The exposure instrument tables therefore map `OA -> effect_allele` and `EA -> other_allele`.",
        "",
        "| Protein | Outcome | SNP | Exposure EA/OA/EAF | Outcome EA/OA/EAF | Harmonisation | OR (95% CI) | P |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for row in rows:
        lines.append(
            f"| {row['protein']} | {row['outcome_id']} | {row['SNP']} | "
            f"{row['exposure_effect_allele']}/{row['exposure_other_allele']}/{row['exposure_eaf']} | "
            f"{row['outcome_effect_allele']}/{row['outcome_other_allele']}/{row['outcome_eaf_raw']} "
            f"(aligned EAF {row['outcome_eaf_aligned_to_exposure_effect']}) | "
            f"{row['harmonise_action']} | {row['wald_or']} ({row['wald_ci']}) | {row['wald_p']} |"
        )
    lines += [
        "",
        "Audit conclusion:",
        "- FGF5: the C allele increases FGF5 level and is associated with higher AF risk; after flipping HF to the C allele, HF direction is also risk-increasing.",
        "- LPA: the T allele increases apolipoprotein(a) level and is associated with higher HF risk; AF shows a nominal risk-increasing association in the same direction.",
        "",
        f"CSV table: `{OUT_CSV}`",
    ]
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote {OUT_CSV}")
    print(f"Wrote {OUT_MD}")


if __name__ == "__main__":
    main()

