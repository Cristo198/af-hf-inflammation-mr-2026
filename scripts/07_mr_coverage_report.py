"""Create extraction and harmonisation coverage report for MR."""

from __future__ import annotations

import csv
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
EXPOSURE = PROJECT_ROOT / "data" / "processed" / "ukbppp_inflammation_cis_pqtl_instruments_lead_no_mhc.csv"
OUTCOME_DIR = PROJECT_ROOT / "data" / "outcomes"
MR_DIR = PROJECT_ROOT / "results" / "mr"
OUT = MR_DIR / "mr_coverage_report.csv"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def main() -> None:
    exposure = read_csv(EXPOSURE)
    exposure_snps = {row["SNP"] for row in exposure}
    rows = []
    for outcome_file in sorted(OUTCOME_DIR.glob("*_outcome_effects_*.csv")):
        outcome_id = outcome_file.name.split("_")[0]
        outcome = read_csv(outcome_file)
        outcome_snps = {row["SNP"] for row in outcome}
        mr_file = MR_DIR / f"{outcome_id}_wald_mr.csv"
        mr_rows = read_csv(mr_file) if mr_file.exists() else []
        mr_snps = {row["SNP"] for row in mr_rows}

        rows.append(
            {
                "outcome_id": outcome_id,
                "exposure_instruments": len(exposure_snps),
                "outcome_effects_extracted": len(outcome_snps),
                "missing_from_outcome": len(exposure_snps - outcome_snps),
                "harmonised_mr_rows": len(mr_snps),
                "dropped_after_extraction": len(outcome_snps - mr_snps),
                "extraction_rate": f"{len(outcome_snps) / len(exposure_snps):.4f}",
                "harmonisation_rate_among_extracted": f"{len(mr_snps) / len(outcome_snps):.4f}" if outcome_snps else "",
            }
        )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="", encoding="utf-8") as handle:
        fieldnames = [
            "outcome_id",
            "exposure_instruments",
            "outcome_effects_extracted",
            "missing_from_outcome",
            "harmonised_mr_rows",
            "dropped_after_extraction",
            "extraction_rate",
            "harmonisation_rate_among_extracted",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()

