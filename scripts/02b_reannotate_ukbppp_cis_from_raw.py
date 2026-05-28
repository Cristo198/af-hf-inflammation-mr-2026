"""Re-annotate downloaded UKB-PPP raw pQTL rows for cis status.

Use this when the raw UKB-PPP pQTL table has already been downloaded and only
gene-coordinate annotation logic changed.
"""

from __future__ import annotations

import csv
import importlib.util
from collections import defaultdict
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCREEN_SCRIPT = PROJECT_ROOT / "scripts" / "02_screen_ukbppp_cis_pqtl.py"
EXPOSURE_LIST = PROJECT_ROOT / "data" / "exposure" / "ukbppp_olink_inflammation_panel_proteins.csv"
RAW_IN = PROJECT_ROOT / "data" / "interim" / "ukbppp_inflammation_panel_pqtl_p5e-8_raw.csv"
CIS_OUT = PROJECT_ROOT / "data" / "processed" / "ukbppp_inflammation_panel_cis_pqtl_p5e-8_1mb.csv"
SUMMARY_OUT = PROJECT_ROOT / "data" / "processed" / "ukbppp_inflammation_panel_cis_pqtl_summary.md"


def load_screen_module():
    spec = importlib.util.spec_from_file_location("screen_ukbppp", SCREEN_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load {SCREEN_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    screen = load_screen_module()
    proteins = read_csv(EXPOSURE_LIST)
    raw_rows = read_csv(RAW_IN)
    cache = screen.load_gene_cache()

    protein_meta = {row["Gene name"]: row for row in proteins}
    raw_by_gene: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in raw_rows:
        raw_by_gene[row["target_gene"]].append(row)

    updated_raw = []
    cis_rows = []
    for target_gene, rows in raw_by_gene.items():
        meta = protein_meta[target_gene]
        genes = []
        for symbol in screen.gene_symbols(target_gene):
            hit = screen.lookup_gene(symbol, cache, delay=0.05)
            if hit is not None:
                genes.append(hit)
        for row in rows:
            is_cis, cis_gene = screen.annotate_cis(row, genes, window=1_000_000)
            row["target_gene_symbols_checked"] = ";".join(gene["symbol"] for gene in genes)
            row["is_cis_1mb"] = "TRUE" if is_cis else "FALSE"
            row["cis_gene_match"] = cis_gene
            row["is_mhc_region"] = "TRUE" if screen.in_mhc(row) else "FALSE"
            row["target_panel"] = meta["Explore 384 panel"]
            row["target_uniprot"] = meta["UniProt ID"]
            updated_raw.append(row)
            if is_cis:
                cis_rows.append(row)

    screen.save_gene_cache(cache)

    fieldnames = [
        "chr",
        "pos",
        "rsid",
        "EA",
        "OA",
        "EAF",
        "olink",
        "olink_target_name",
        "uniprot",
        "gene_symbol",
        "locus_annotation",
        "gene_ensembl",
        "beta",
        "SE",
        "p",
        "log10",
        "target_uniprot",
        "target_gene",
        "target_panel",
        "target_gene_symbols_checked",
        "cis_window_bp",
        "is_cis_1mb",
        "cis_gene_match",
        "is_mhc_region",
    ]
    for row in updated_raw:
        row["cis_window_bp"] = "1000000"

    write_csv(RAW_IN, updated_raw, fieldnames)
    write_csv(CIS_OUT, cis_rows, fieldnames)

    proteins_with_assoc = set(raw_by_gene)
    no_assoc_count = len(proteins) - len(proteins_with_assoc)
    no_coordinate_count = sum(
        1
        for gene in proteins_with_assoc
        if not any(symbol in cache and cache[symbol] for symbol in screen.gene_symbols(gene))
    )
    SUMMARY_OUT.write_text(
        "\n".join(
            [
                "# UKB-PPP Inflammation Panel cis-pQTL Screening Summary",
                "",
                "Population: EUR",
                "p-value threshold: 5e-8",
                "cis window: 1000000 bp",
                f"Targets queried: {len(proteins)}",
                f"Raw significant pQTL rows: {len(updated_raw)}",
                f"proteins with at least one significant pQTL: {len(proteins_with_assoc)}",
                f"proteins without significant pQTL at this threshold: {no_assoc_count}",
                f"cis-pQTL rows: {len(cis_rows)}",
                f"proteins with at least one cis-pQTL: {len({row['target_gene'] for row in cis_rows})}",
                f"proteins with significant pQTL but no mapped Ensembl coordinates: {no_coordinate_count}",
                "failed target queries: 0",
                "",
                f"Raw pQTL table: `{RAW_IN}`",
                f"cis-pQTL table: `{CIS_OUT}`",
            ]
        ),
        encoding="utf-8",
    )
    print(f"Reannotated {len(updated_raw)} raw rows")
    print(f"Wrote {len(cis_rows)} cis rows to {CIS_OUT}")


if __name__ == "__main__":
    main()

