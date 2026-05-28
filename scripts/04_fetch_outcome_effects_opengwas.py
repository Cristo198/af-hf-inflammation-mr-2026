"""Fetch AF/HF outcome effects for exposure SNPs from OpenGWAS.

Requires internet access to https://api.opengwas.io and usually an OpenGWAS JWT.
Set OPENGWAS_JWT in the environment before running when required by the API.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[1]
EXPOSURE = PROJECT_ROOT / "data" / "processed" / "ukbppp_inflammation_cis_pqtl_instruments_lead_no_mhc.csv"
OUT_DIR = PROJECT_ROOT / "data" / "outcomes"
MANIFEST = PROJECT_ROOT / "config" / "outcomes.tsv"
API_URL = "https://api.opengwas.io/api/associations"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def chunks(values: list[str], size: int) -> list[list[str]]:
    return [values[i : i + size] for i in range(0, len(values), size)]


def post_query(url: str, params: dict[str, Any], jwt: str | None, retries: int = 4) -> Any:
    query = urllib.parse.urlencode(params, doseq=True)
    target = f"{url}?{query}"
    headers = {
        "Accept": "application/json",
        "User-Agent": "cardio-no-lab-project/0.1",
    }
    if jwt:
        headers["Authorization"] = f"Bearer {jwt}"

    last_error: Exception | None = None
    for attempt in range(retries):
        try:
            request = urllib.request.Request(target, data=b"", headers=headers, method="POST")
            with urllib.request.urlopen(request, timeout=120) as response:
                return json.loads(response.read().decode("utf-8"))
        except urllib.error.HTTPError as error:
            last_error = error
            if error.code == 429:
                retry_after = error.headers.get("retry-after")
                wait = float(retry_after) if retry_after else 5 * (attempt + 1)
                time.sleep(wait)
                continue
            raise
        except Exception as error:  # noqa: BLE001
            last_error = error
            time.sleep(3 * (attempt + 1))
    raise RuntimeError(f"OpenGWAS request failed after {retries} attempts") from last_error


def standardise(row: dict[str, Any], outcome: dict[str, str]) -> dict[str, Any]:
    return {
        "outcome_id": outcome["outcome_id"],
        "outcome_name": outcome["outcome_name"],
        "opengwas_id": outcome["opengwas_id"],
        "SNP": row.get("rsid") or row.get("variant") or row.get("SNP"),
        "chr": row.get("chr") or row.get("chromosome"),
        "pos": row.get("position") or row.get("pos"),
        "effect_allele": row.get("ea") or row.get("effect_allele"),
        "other_allele": row.get("nea") or row.get("other_allele"),
        "eaf": row.get("eaf"),
        "beta": row.get("beta"),
        "se": row.get("se"),
        "pval": row.get("p") or row.get("pval"),
        "source": "OpenGWAS API",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--batch-size", type=int, default=100)
    parser.add_argument("--outcomes", default="AF,HF", help="comma-separated outcome_id list from config/outcomes.tsv")
    args = parser.parse_args()

    exposure = read_csv(EXPOSURE)
    snps = sorted({row["SNP"] for row in exposure})
    wanted = set(args.outcomes.split(","))
    outcomes = [row for row in read_tsv(MANIFEST) if row["outcome_id"] in wanted]
    jwt = os.environ.get("OPENGWAS_JWT")

    for outcome in outcomes:
        rows: list[dict[str, Any]] = []
        for batch in chunks(snps, args.batch_size):
            params = {
                "variant": batch,
                "id": [outcome["opengwas_id"]],
                "proxies": 0,
            }
            data = post_query(API_URL, params, jwt)
            if isinstance(data, dict) and "data" in data:
                data_rows = data["data"]
            else:
                data_rows = data
            rows.extend(standardise(row, outcome) for row in data_rows)
            time.sleep(0.25)

        out_path = OUT_DIR / f"{outcome['outcome_id']}_outcome_effects_opengwas.csv"
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
        write_csv(out_path, rows, fieldnames)
        print(f"Wrote {len(rows)} rows to {out_path}")


if __name__ == "__main__":
    main()
