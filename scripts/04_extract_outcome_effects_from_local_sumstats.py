"""Extract exposure SNPs from local AF/HF GWAS summary-statistics files.

Place downloaded outcome files in `data/raw/outcomes/` and run this script.
It supports plain text, .gz, and .zip files and tries to infer common GWAS
summary-statistics column names.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import zipfile
from pathlib import Path
from typing import Iterable


PROJECT_ROOT = Path(__file__).resolve().parents[1]
EXPOSURE = PROJECT_ROOT / "data" / "processed" / "ukbppp_inflammation_cis_pqtl_instruments_lead_no_mhc.csv"
MANIFEST = PROJECT_ROOT / "config" / "outcomes.tsv"
RAW_DIR = PROJECT_ROOT / "data" / "raw" / "outcomes"
OUT_DIR = PROJECT_ROOT / "data" / "outcomes"


COLUMN_ALIASES = {
    "snp": ["snp", "rsid", "rs_id", "markername", "marker", "variant", "variant_id", "id"],
    "chr": ["chr", "chrom", "chromosome"],
    "pos": ["pos", "position", "base_pair_location", "bp"],
    "effect_allele": ["effect_allele", "ea", "a1", "allele1", "alt", "tested_allele", "coded_allele"],
    "other_allele": ["other_allele", "nea", "a2", "allele2", "ref", "non_effect_allele", "noncoded_allele"],
    "eaf": ["eaf", "effect_allele_frequency", "af", "freq", "freq1"],
    "beta": ["beta", "effect", "b", "log_odds", "logor", "estimate"],
    "se": ["se", "stderr", "standard_error", "sebeta"],
    "pval": ["p", "pval", "p_value", "p-value", "pvalue", "p.value"],
}

OUTCOME_COLUMN_MAPS = {
    "AF": {
        "snp": "rs_dbSNP147",
        "marker": "MarkerName",
        "chr": "CHR",
        "pos": "POS_GRCh37",
        "effect_allele": "A2",
        "other_allele": "A1",
        "eaf": "Freq_A2",
        "beta": "Effect_A2",
        "se": "StdErr",
        "pval": "Pvalue",
    },
    "HF": {
        "snp": "SNP",
        "chr": "CHR",
        "pos": "BP",
        "effect_allele": "A1",
        "other_allele": "A2",
        "eaf": "freq",
        "beta": "b",
        "se": "se",
        "pval": "p",
    },
}


def normalise(name: str) -> str:
    return name.strip().lower().replace("-", "_").replace(".", "_")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def open_text(path: Path) -> Iterable[str]:
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            yield from handle
    elif path.suffix == ".zip":
        with zipfile.ZipFile(path) as archive:
            members = [name for name in archive.namelist() if not name.endswith("/")]
            if not members:
                raise FileNotFoundError(f"No file inside {path}")
            member = members[0]
            with archive.open(member) as raw:
                text = io.TextIOWrapper(raw, encoding="utf-8", errors="replace")
                yield from text
    else:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            yield from handle


def detect_delimiter(header_line: str) -> str:
    if "\t" in header_line:
        return "\t"
    if "," in header_line:
        return ","
    return None  # type: ignore[return-value]


def map_columns(fieldnames: list[str], outcome_id: str) -> dict[str, str]:
    if outcome_id in OUTCOME_COLUMN_MAPS:
        configured = OUTCOME_COLUMN_MAPS[outcome_id]
        field_set = set(fieldnames)
        missing = [value for key, value in configured.items() if key != "marker" and value not in field_set]
        if missing:
            raise ValueError(f"Configured columns missing for {outcome_id}: {missing}. Found columns: {fieldnames[:30]}")
        return configured

    norm_to_original = {normalise(name): name for name in fieldnames}
    mapping = {}
    for target, aliases in COLUMN_ALIASES.items():
        for alias in aliases:
            if normalise(alias) in norm_to_original:
                mapping[target] = norm_to_original[normalise(alias)]
                break
    required = ["snp", "beta", "se", "effect_allele", "other_allele"]
    missing = [name for name in required if name not in mapping]
    if missing:
        raise ValueError(f"Missing required columns {missing}. Found columns: {fieldnames[:30]}")
    return mapping


def exposure_keys() -> tuple[set[str], dict[str, str], dict[str, dict[str, str]]]:
    rows = read_csv(EXPOSURE)
    snps = {row["SNP"] for row in rows}
    chrpos_to_snp = {f"{row['chr']}:{row['pos']}": row["SNP"] for row in rows}
    by_snp = {row["SNP"]: row for row in rows}
    return snps, chrpos_to_snp, by_snp


def norm_allele(value: str) -> str:
    return (value or "").upper()


def allele_score(row: dict[str, str], exposure: dict[str, str] | None) -> int:
    if exposure is None:
        return 0
    ea_exp = norm_allele(exposure["effect_allele"])
    oa_exp = norm_allele(exposure["other_allele"])
    ea_out = norm_allele(row["effect_allele"])
    oa_out = norm_allele(row["other_allele"])
    if ea_exp == ea_out and oa_exp == oa_out:
        return 2
    if ea_exp == oa_out and oa_exp == ea_out:
        return 2
    if ea_exp in {ea_out, oa_out} or oa_exp in {ea_out, oa_out}:
        return 1
    return 0


def deduplicate_rows(rows: list[dict[str, str]], exposure_by_snp: dict[str, dict[str, str]]) -> list[dict[str, str]]:
    grouped: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        grouped.setdefault(row["SNP"], []).append(row)

    deduped = []
    for snp, group in grouped.items():
        exposure = exposure_by_snp.get(snp)
        exposure_chrpos = f"{exposure['chr']}:{exposure['pos']}" if exposure else ""

        def key(row: dict[str, str]) -> tuple[int, int, float]:
            row_chrpos = f"{row['chr']}:{row['pos']}"
            pos_match = 1 if exposure_chrpos and row_chrpos == exposure_chrpos else 0
            score = allele_score(row, exposure)
            try:
                p_value = float(row["pval"])
            except Exception:
                p_value = 1.0
            return (pos_match, score, -p_value)

        deduped.append(max(group, key=key))
    return deduped


def extract_file(path: Path, outcome: dict[str, str]) -> list[dict[str, str]]:
    snps, chrpos_to_snp, exposure_by_snp = exposure_keys()
    iterator = iter(open_text(path))
    header = ""
    for line in iterator:
        if line.strip() and not line.startswith("#"):
            header = line
            break
    if not header:
        raise ValueError(f"No header found in {path}")
    delimiter = detect_delimiter(header)
    fieldnames = next(csv.reader([header], delimiter=delimiter))
    mapping = map_columns(fieldnames, outcome["outcome_id"])

    rows = []
    reader = csv.DictReader(iterator, fieldnames=fieldnames, delimiter=delimiter)
    for row in reader:
        raw_snp = row.get(mapping["snp"], "")
        candidate_snps = [value.strip() for value in raw_snp.replace(",", ";").split(";") if value.strip()]
        snp = next((value for value in candidate_snps if value in snps), raw_snp)
        marker = row.get(mapping.get("marker", ""), "")
        pos_key = ""
        if "chr" in mapping and "pos" in mapping:
            pos_key = f"{row.get(mapping['chr'], '')}:{row.get(mapping['pos'], '')}"
        matched_snp = snp if snp in snps else chrpos_to_snp.get(pos_key, "")
        if not matched_snp and marker in snps:
            matched_snp = marker
        if not matched_snp:
            continue
        rows.append(
            {
                "outcome_id": outcome["outcome_id"],
                "outcome_name": outcome["outcome_name"],
                "opengwas_id": outcome["opengwas_id"],
                "SNP": matched_snp,
                "chr": row.get(mapping.get("chr", ""), ""),
                "pos": row.get(mapping.get("pos", ""), ""),
                "effect_allele": row.get(mapping["effect_allele"], ""),
                "other_allele": row.get(mapping["other_allele"], ""),
                "eaf": row.get(mapping.get("eaf", ""), ""),
                "beta": row.get(mapping["beta"], ""),
                "se": row.get(mapping["se"], ""),
                "pval": row.get(mapping.get("pval", ""), ""),
                "source": str(path),
            }
        )
    return deduplicate_rows(rows, exposure_by_snp)


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    fieldnames = [
        "outcome_id",
        "outcome_name",
        "opengwas_id",
        "SNP",
        "chr",
        "pos",
        "effect_allele",
        "other_allele",
        "eaf",
        "beta",
        "se",
        "pval",
        "source",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def locate_file(filename: str) -> Path | None:
    candidates = [
        RAW_DIR / filename,
        RAW_DIR / f"{filename}.gz",
        RAW_DIR / f"{filename}.zip",
        RAW_DIR / filename.replace(".txt", ".tbl.gz"),
        RAW_DIR / filename.replace(".txt", ".gz"),
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outcomes", default="AF,HF")
    args = parser.parse_args()

    wanted = set(args.outcomes.split(","))
    outcomes = [row for row in read_tsv(MANIFEST) if row["outcome_id"] in wanted]
    missing = []
    for outcome in outcomes:
        path = locate_file(outcome["local_filename"])
        if path is None:
            missing.append(outcome)
            continue
        rows = extract_file(path, outcome)
        out_path = OUT_DIR / f"{outcome['outcome_id']}_outcome_effects_local.csv"
        write_csv(out_path, rows)
        print(f"{outcome['outcome_id']}: extracted {len(rows)} rows from {path} -> {out_path}")

    if missing:
        print("Missing local outcome files:")
        for outcome in missing:
            print(f"- {outcome['outcome_id']}: place `{outcome['local_filename']}` in {RAW_DIR}")
            print(f"  Source: {outcome['source_url']}")


if __name__ == "__main__":
    main()
