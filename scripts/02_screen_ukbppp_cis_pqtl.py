"""Fetch UKB-PPP pQTL associations and screen cis-pQTLs.

This script queries the public UKB-PPP browser endpoint by Olink protein target,
then annotates genome-wide significant associations as cis if the variant lies
within 1 Mb of the target gene coordinates from Ensembl GRCh38.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import ssl
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[1]
EXPOSURE_LIST = PROJECT_ROOT / "data" / "exposure" / "ukbppp_olink_inflammation_panel_proteins.csv"
GENE_CACHE = PROJECT_ROOT / "data" / "metadata" / "ensembl_gene_coordinates_cache.json"
RAW_OUT = PROJECT_ROOT / "data" / "interim" / "ukbppp_inflammation_panel_pqtl_p5e-8_raw.csv"
CIS_OUT = PROJECT_ROOT / "data" / "processed" / "ukbppp_inflammation_panel_cis_pqtl_p5e-8_1mb.csv"
SUMMARY_OUT = PROJECT_ROOT / "data" / "processed" / "ukbppp_inflammation_panel_cis_pqtl_summary.md"

UKBPPP_ENDPOINT = "https://metabolomics.helmholtz-munich.de/ukbbpgwas/processing/retr.php"
ENSEMBL_LOOKUP = "https://rest.ensembl.org/lookup/symbol/homo_sapiens/{symbol}"
ENSEMBL_XREF = "https://rest.ensembl.org/xrefs/symbol/homo_sapiens/{symbol}"
ENSEMBL_LOOKUP_ID = "https://rest.ensembl.org/lookup/id/{ensembl_id}"

SSL_CONTEXT = ssl._create_unverified_context()
MHC_CHR = "6"
MHC_START = 25_500_000
MHC_END = 34_000_000


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def load_gene_cache() -> dict[str, Any]:
    if not GENE_CACHE.exists():
        return {}
    return json.loads(GENE_CACHE.read_text(encoding="utf-8"))


def save_gene_cache(cache: dict[str, Any]) -> None:
    GENE_CACHE.parent.mkdir(parents=True, exist_ok=True)
    GENE_CACHE.write_text(json.dumps(cache, indent=2, ensure_ascii=False), encoding="utf-8")


def fetch_json(url: str, retries: int = 4, pause: float = 1.0) -> dict[str, Any]:
    last_error: Exception | None = None
    for attempt in range(retries):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "cardio-no-lab-project/0.1"})
            with urllib.request.urlopen(request, context=SSL_CONTEXT, timeout=60) as response:
                return json.loads(response.read().decode("utf-8"))
        except urllib.error.HTTPError as error:
            last_error = error
            if error.code == 429:
                retry_after = error.headers.get("retry-after")
                wait = float(retry_after) if retry_after else pause * (attempt + 2)
                time.sleep(wait)
                continue
            if error.code in {500, 502, 503, 504}:
                time.sleep(pause * (attempt + 2))
                continue
            raise
        except Exception as error:  # noqa: BLE001
            last_error = error
            time.sleep(pause * (attempt + 2))
    raise RuntimeError(f"Failed after {retries} attempts: {url}") from last_error


def gene_symbols(gene_name: str) -> list[str]:
    parts = re.split(r"[_;/, ]+", gene_name.strip())
    return [part for part in parts if part and not part.lower().startswith("isoform")]


def lookup_gene(symbol: str, cache: dict[str, Any], delay: float) -> dict[str, Any] | None:
    if symbol in cache:
        result = cache[symbol]
        return None if result is None else result

    url = ENSEMBL_LOOKUP.format(symbol=urllib.parse.quote(symbol)) + "?content-type=application/json;expand=0"
    try:
        result = fetch_json(url, retries=3, pause=1.0)
        cache[symbol] = {
            "symbol": symbol,
            "ensembl_id": result.get("id"),
            "chr": str(result.get("seq_region_name")),
            "start": int(result.get("start")),
            "end": int(result.get("end")),
            "assembly": result.get("assembly_name"),
            "biotype": result.get("biotype"),
        }
    except Exception:
        try:
            xref_url = (
                ENSEMBL_XREF.format(symbol=urllib.parse.quote(symbol))
                + "?external_db=HGNC;content-type=application/json"
            )
            xrefs = fetch_json(xref_url, retries=3, pause=1.0)
            ensembl_id = xrefs[0]["id"]
            id_url = (
                ENSEMBL_LOOKUP_ID.format(ensembl_id=urllib.parse.quote(ensembl_id))
                + "?content-type=application/json;expand=0"
            )
            result = fetch_json(id_url, retries=3, pause=1.0)
            cache[symbol] = {
                "symbol": symbol,
                "ensembl_id": result.get("id"),
                "chr": str(result.get("seq_region_name")),
                "start": int(result.get("start")),
                "end": int(result.get("end")),
                "assembly": result.get("assembly_name"),
                "biotype": result.get("biotype"),
                "display_name": result.get("display_name"),
            }
        except Exception:
            cache[symbol] = None
    time.sleep(delay)
    return None if cache[symbol] is None else cache[symbol]


def fetch_target_associations(target_name: str, pop: str, p_threshold: str) -> list[dict[str, Any]]:
    params = {
        "m": "target",
        "t": target_name,
        "pop": pop,
        "trait": "",
        "pma": p_threshold,
    }
    url = UKBPPP_ENDPOINT + "?" + urllib.parse.urlencode(params)
    payload = fetch_json(url)
    header = payload.get("header", [])
    rows = payload.get("results", [])
    return [dict(zip(header, row)) for row in rows]


def in_mhc(row: dict[str, Any]) -> bool:
    try:
        return str(row.get("chr")) == MHC_CHR and MHC_START <= int(row.get("pos")) <= MHC_END
    except Exception:
        return False


def annotate_cis(row: dict[str, Any], genes: list[dict[str, Any]], window: int) -> tuple[bool, str]:
    for gene in genes:
        try:
            if str(row.get("chr")) != str(gene["chr"]):
                continue
            pos = int(row.get("pos"))
            if int(gene["start"]) - window <= pos <= int(gene["end"]) + window:
                return True, gene["symbol"]
        except Exception:
            continue
    return False, ""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pop", default="EUR", help="UKB-PPP population code, default EUR")
    parser.add_argument("--p-threshold", default="5e-8", help="p-value threshold sent to UKB-PPP endpoint")
    parser.add_argument("--window", type=int, default=1_000_000, help="cis window around gene start/end")
    parser.add_argument("--delay", type=float, default=0.15, help="delay between web requests")
    parser.add_argument("--max-targets", type=int, default=0, help="limit targets for testing")
    args = parser.parse_args()

    proteins = read_csv(EXPOSURE_LIST)
    if args.max_targets:
        proteins = proteins[: args.max_targets]

    cache = load_gene_cache()
    raw_rows: list[dict[str, Any]] = []
    cis_rows: list[dict[str, Any]] = []
    failures: list[str] = []

    for index, protein in enumerate(proteins, start=1):
        target_name = protein["Protein name"]
        target_gene = protein["Gene name"]
        panel = protein["Explore 384 panel"]
        print(f"[{index}/{len(proteins)}] {target_gene} | {target_name}")

        try:
            associations = fetch_target_associations(target_name, args.pop, args.p_threshold)
        except Exception as error:  # noqa: BLE001
            failures.append(f"{target_gene}\t{target_name}\t{type(error).__name__}: {error}")
            continue

        genes = []
        if associations:
            for symbol in gene_symbols(target_gene):
                hit = lookup_gene(symbol, cache, args.delay)
                if hit is not None:
                    genes.append(hit)

        for row in associations:
            is_cis, cis_gene = annotate_cis(row, genes, args.window)
            enriched = {
                **row,
                "target_uniprot": protein["UniProt ID"],
                "target_gene": target_gene,
                "target_panel": panel,
                "target_gene_symbols_checked": ";".join(gene["symbol"] for gene in genes),
                "cis_window_bp": args.window,
                "is_cis_1mb": "TRUE" if is_cis else "FALSE",
                "cis_gene_match": cis_gene,
                "is_mhc_region": "TRUE" if in_mhc(row) else "FALSE",
            }
            raw_rows.append(enriched)
            if is_cis:
                cis_rows.append(enriched)

        time.sleep(args.delay)
        if index % 25 == 0:
            save_gene_cache(cache)

    save_gene_cache(cache)

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
    write_csv(RAW_OUT, raw_rows, fieldnames)
    write_csv(CIS_OUT, cis_rows, fieldnames)

    protein_count = len({row["target_gene"] for row in cis_rows})
    proteins_with_assoc = {row["target_gene"] for row in raw_rows}
    no_coordinate_count = sum(
        1
        for protein in proteins
        if protein["Gene name"] in proteins_with_assoc
        and not any(symbol in cache and cache[symbol] for symbol in gene_symbols(protein["Gene name"]))
    )
    no_assoc_count = len(proteins) - len(proteins_with_assoc)
    SUMMARY_OUT.parent.mkdir(parents=True, exist_ok=True)
    SUMMARY_OUT.write_text(
        "\n".join(
            [
                "# UKB-PPP Inflammation Panel cis-pQTL Screening Summary",
                "",
                f"Population: {args.pop}",
                f"p-value threshold: {args.p_threshold}",
                f"cis window: {args.window} bp",
                f"Targets queried: {len(proteins)}",
                f"Raw significant pQTL rows: {len(raw_rows)}",
                f"proteins with at least one significant pQTL: {len(proteins_with_assoc)}",
                f"proteins without significant pQTL at this threshold: {no_assoc_count}",
                f"cis-pQTL rows: {len(cis_rows)}",
                f"proteins with at least one cis-pQTL: {protein_count}",
                f"proteins with significant pQTL but no mapped Ensembl coordinates: {no_coordinate_count}",
                f"failed target queries: {len(failures)}",
                "",
                f"Raw pQTL table: `{RAW_OUT}`",
                f"cis-pQTL table: `{CIS_OUT}`",
            ]
        ),
        encoding="utf-8",
    )
    if failures:
        fail_path = PROJECT_ROOT / "data" / "interim" / "ukbppp_inflammation_panel_pqtl_failures.tsv"
        fail_path.parent.mkdir(parents=True, exist_ok=True)
        fail_path.write_text("\n".join(failures), encoding="utf-8")

    print(f"Wrote {len(raw_rows)} raw rows to {RAW_OUT}")
    print(f"Wrote {len(cis_rows)} cis rows to {CIS_OUT}")


if __name__ == "__main__":
    main()
