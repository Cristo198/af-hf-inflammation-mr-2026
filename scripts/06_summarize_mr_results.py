"""Summarise preliminary Wald-ratio MR results."""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
MR_ALL = PROJECT_ROOT / "results" / "mr" / "wald_mr_all_outcomes.csv"
OUT_DIR = PROJECT_ROOT / "results" / "mr"
SUMMARY_MD = OUT_DIR / "mr_summary_2026-05-26.md"
SHARED_OUT = OUT_DIR / "shared_candidate_preliminary.csv"
FDR_OUT = OUT_DIR / "fdr_significant_preliminary.csv"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def flt(value: str) -> float:
    return float(value) if value not in {"", "NA", "nan"} else float("nan")


def top_rows(rows: list[dict[str, str]], n: int = 10) -> list[dict[str, str]]:
    return sorted(rows, key=lambda row: flt(row["pval_mr"]))[:n]


def fmt_result(row: dict[str, str]) -> str:
    return (
        f"{row['protein']} ({row['protein_name']}), OR {row['or']} "
        f"[{row['or_lci95']}, {row['or_uci95']}], "
        f"P={row['pval_mr']}, FDR={row['fdr']}"
    )


def main() -> None:
    rows = read_csv(MR_ALL)
    by_outcome: dict[str, list[dict[str, str]]] = defaultdict(list)
    by_protein: dict[str, dict[str, dict[str, str]]] = defaultdict(dict)
    for row in rows:
        by_outcome[row["outcome_id"]].append(row)
        by_protein[row["protein"]][row["outcome_id"]] = row

    fdr_rows = [row for row in rows if row["fdr"] and flt(row["fdr"]) < 0.05]
    nominal_rows = [row for row in rows if row["pval_mr"] and flt(row["pval_mr"]) < 0.05]

    shared = []
    for protein, outcome_rows in by_protein.items():
        if "AF" not in outcome_rows or "HF" not in outcome_rows:
            continue
        af = outcome_rows["AF"]
        hf = outcome_rows["HF"]
        af_beta = flt(af["beta_mr"])
        hf_beta = flt(hf["beta_mr"])
        same_direction = af_beta * hf_beta > 0
        af_fdr = flt(af["fdr"]) < 0.05
        hf_fdr = flt(hf["fdr"]) < 0.05
        af_nom = flt(af["pval_mr"]) < 0.05
        hf_nom = flt(hf["pval_mr"]) < 0.05
        if same_direction and ((af_fdr and hf_nom) or (hf_fdr and af_nom)):
            shared.append(
                {
                    "protein": protein,
                    "protein_name": af["protein_name"],
                    "panel": af["panel"],
                    "SNP": af["SNP"],
                    "AF_OR": af["or"],
                    "AF_LCI95": af["or_lci95"],
                    "AF_UCI95": af["or_uci95"],
                    "AF_P": af["pval_mr"],
                    "AF_FDR": af["fdr"],
                    "HF_OR": hf["or"],
                    "HF_LCI95": hf["or_lci95"],
                    "HF_UCI95": hf["or_uci95"],
                    "HF_P": hf["pval_mr"],
                    "HF_FDR": hf["fdr"],
                    "direction": "same_risk" if af_beta > 0 else "same_protective",
                    "shared_rule": "AF_FDR+HF_nominal" if af_fdr else "HF_FDR+AF_nominal",
                }
            )
    shared = sorted(shared, key=lambda row: min(flt(row["AF_P"]), flt(row["HF_P"])))

    fdr_fieldnames = [
        "protein",
        "protein_name",
        "panel",
        "outcome_id",
        "outcome_name",
        "SNP",
        "or",
        "or_lci95",
        "or_uci95",
        "pval_mr",
        "fdr",
        "harmonise_action",
    ]
    shared_fieldnames = [
        "protein",
        "protein_name",
        "panel",
        "SNP",
        "AF_OR",
        "AF_LCI95",
        "AF_UCI95",
        "AF_P",
        "AF_FDR",
        "HF_OR",
        "HF_LCI95",
        "HF_UCI95",
        "HF_P",
        "HF_FDR",
        "direction",
        "shared_rule",
    ]
    write_csv(FDR_OUT, fdr_rows, fdr_fieldnames)
    write_csv(SHARED_OUT, shared, shared_fieldnames)

    lines = [
        "# Preliminary Wald-ratio MR Summary",
        "",
        "Date: 2026-05-26",
        "",
        "Input:",
        "- Exposure: 529 UKB-PPP inflammation-panel lead cis-pQTL instruments after MHC exclusion.",
        "- Outcomes: Nielsen 2018 AF and HERMES 2019 HF local summary statistics.",
        "- Method: single-variant Wald ratio; FDR by Benjamini-Hochberg within outcome.",
        "",
        "Coverage:",
        f"- AF outcome effects extracted: {len(by_outcome.get('AF', []))} harmonised MR rows.",
        f"- HF outcome effects extracted: {len(by_outcome.get('HF', []))} harmonised MR rows.",
        "",
        "Significance counts:",
    ]
    for outcome_id in ["AF", "HF"]:
        subset = by_outcome.get(outcome_id, [])
        nominal = [row for row in subset if flt(row["pval_mr"]) < 0.05]
        fdr = [row for row in subset if flt(row["fdr"]) < 0.05]
        lines.append(f"- {outcome_id}: nominal P<0.05 = {len(nominal)}, FDR<0.05 = {len(fdr)}")

    lines += [
        "",
        "FDR-significant AF results:",
        *(f"- {fmt_result(row)}" for row in top_rows([row for row in fdr_rows if row["outcome_id"] == "AF"], 20)),
        "",
        "FDR-significant HF results:",
        *(f"- {fmt_result(row)}" for row in top_rows([row for row in fdr_rows if row["outcome_id"] == "HF"], 20)),
        "",
        "Preliminary shared candidates:",
    ]
    if shared:
        for row in shared:
            lines.append(
                f"- {row['protein']} ({row['protein_name']}): "
                f"AF OR {row['AF_OR']} [{row['AF_LCI95']}, {row['AF_UCI95']}], "
                f"P={row['AF_P']}, FDR={row['AF_FDR']}; "
                f"HF OR {row['HF_OR']} [{row['HF_LCI95']}, {row['HF_UCI95']}], "
                f"P={row['HF_P']}, FDR={row['HF_FDR']}; {row['shared_rule']}, {row['direction']}."
            )
    else:
        lines.append("- None by the current rule.")

    lines += [
        "",
        "Important interpretation note:",
        "These are first-pass Wald-ratio results. Shared candidates require replication, allele-coding audit, LD/proxy review, and colocalization before being interpreted as causal targets.",
        "",
        f"FDR table: `{FDR_OUT}`",
        f"Shared candidate table: `{SHARED_OUT}`",
    ]
    SUMMARY_MD.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote summary to {SUMMARY_MD}")
    print(f"FDR rows: {len(fdr_rows)} -> {FDR_OUT}")
    print(f"Shared candidates: {len(shared)} -> {SHARED_OUT}")


if __name__ == "__main__":
    main()

