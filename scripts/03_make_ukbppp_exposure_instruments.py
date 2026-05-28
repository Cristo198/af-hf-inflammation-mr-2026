"""Prepare UKB-PPP cis-pQTL exposure instrument tables for MR."""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[1]
CIS_IN = PROJECT_ROOT / "data" / "processed" / "ukbppp_inflammation_panel_cis_pqtl_p5e-8_1mb.csv"
ALL_OUT = PROJECT_ROOT / "data" / "processed" / "ukbppp_inflammation_cis_pqtl_instruments_all.csv"
LEAD_OUT = PROJECT_ROOT / "data" / "processed" / "ukbppp_inflammation_cis_pqtl_instruments_lead.csv"
LEAD_NO_MHC_OUT = PROJECT_ROOT / "data" / "processed" / "ukbppp_inflammation_cis_pqtl_instruments_lead_no_mhc.csv"
SUMMARY_OUT = PROJECT_ROOT / "data" / "processed" / "ukbppp_inflammation_cis_pqtl_instruments_summary.md"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def numeric(value: str) -> float:
    try:
        return float(value)
    except Exception:
        return float("nan")


def p_sort_key(row: dict[str, str]) -> tuple[float, float]:
    p = numeric(row["p"])
    log10 = numeric(row["log10"])
    if p == 0:
        return (0.0, -log10)
    return (p, -log10)


def format_row(row: dict[str, str]) -> dict[str, Any]:
    beta = numeric(row["beta"])
    se = numeric(row["SE"])
    f_stat = (beta / se) ** 2 if se and se == se else ""
    return {
        "protein": row["target_gene"],
        "protein_name": row["olink_target_name"],
        "uniprot": row["target_uniprot"],
        "panel": row["target_panel"],
        "SNP": row["rsid"],
        "chr": row["chr"],
        "pos": row["pos"],
        # UKB-PPP JSON exposes columns named EA/OA, but the browser labels
        # them oppositely and EAF matches OA. Treat OA as the effect allele.
        "effect_allele": row["OA"],
        "other_allele": row["EA"],
        "eaf": row["EAF"],
        "beta": row["beta"],
        "se": row["SE"],
        "pval": row["p"],
        "minus_log10_p": row["log10"],
        "f_stat": f"{f_stat:.3f}" if f_stat != "" else "",
        "olink_id": row["olink"],
        "ukbppp_gene_symbol": row["gene_symbol"],
        "ensembl_id": row["gene_ensembl"],
        "cis_gene_match": row["cis_gene_match"],
        "cis_window_bp": row["cis_window_bp"],
        "is_mhc_region": row["is_mhc_region"],
        "locus_annotation": row["locus_annotation"],
    }


def main() -> None:
    cis_rows = read_csv(CIS_IN)
    instruments = [format_row(row) for row in cis_rows]

    by_protein: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in cis_rows:
        by_protein[row["target_gene"]].append(row)

    lead_rows = [format_row(sorted(rows, key=p_sort_key)[0]) for rows in by_protein.values()]
    lead_rows = sorted(lead_rows, key=lambda row: (row["panel"], row["protein"]))
    lead_no_mhc = [row for row in lead_rows if row["is_mhc_region"] != "TRUE"]

    fieldnames = [
        "protein",
        "protein_name",
        "uniprot",
        "panel",
        "SNP",
        "chr",
        "pos",
        "effect_allele",
        "other_allele",
        "eaf",
        "beta",
        "se",
        "pval",
        "minus_log10_p",
        "f_stat",
        "olink_id",
        "ukbppp_gene_symbol",
        "ensembl_id",
        "cis_gene_match",
        "cis_window_bp",
        "is_mhc_region",
        "locus_annotation",
    ]

    write_csv(ALL_OUT, instruments, fieldnames)
    write_csv(LEAD_OUT, lead_rows, fieldnames)
    write_csv(LEAD_NO_MHC_OUT, lead_no_mhc, fieldnames)

    weak = [row for row in lead_rows if row["f_stat"] and float(row["f_stat"]) <= 10]
    mhc = [row for row in lead_rows if row["is_mhc_region"] == "TRUE"]
    SUMMARY_OUT.write_text(
        "\n".join(
            [
                "# UKB-PPP Inflammation cis-pQTL Instrument Summary",
                "",
                f"All cis-pQTL instrument rows: {len(instruments)}",
                f"Proteins with at least one cis-pQTL: {len(lead_rows)}",
                f"Lead instruments after MHC exclusion: {len(lead_no_mhc)}",
                f"Lead instruments in MHC region: {len(mhc)}",
                f"Lead instruments with F <= 10: {len(weak)}",
                "",
                "Output files:",
                f"- All cis-pQTL instruments: `{ALL_OUT}`",
                f"- Lead cis-pQTL instruments: `{LEAD_OUT}`",
                f"- Lead cis-pQTL instruments excluding MHC: `{LEAD_NO_MHC_OUT}`",
                "",
                "Allele note:",
                "The UKB-PPP browser JSON reports fields named `EA` and `OA`, but the table display labels `OA` as the effect allele and `EA` as the other allele. EAF also matches the `OA` frequency for audited variants. Therefore this project maps JSON `OA` to `effect_allele` and JSON `EA` to `other_allele`.",
            ]
        ),
        encoding="utf-8",
    )

    print(f"Wrote {len(instruments)} all instruments to {ALL_OUT}")
    print(f"Wrote {len(lead_rows)} lead instruments to {LEAD_OUT}")
    print(f"Wrote {len(lead_no_mhc)} lead instruments excluding MHC to {LEAD_NO_MHC_OUT}")


if __name__ == "__main__":
    main()
