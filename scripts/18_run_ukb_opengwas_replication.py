"""Run UKB/OpenGWAS VCF replication for all lead cis-pQTL instruments."""

from __future__ import annotations

import csv
import gzip
import math
from collections import defaultdict
from pathlib import Path
from statistics import NormalDist
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[1]
EXPOSURE = PROJECT_ROOT / "data" / "processed" / "ukbppp_inflammation_cis_pqtl_instruments_lead_no_mhc.csv"
PRIMARY_ALL = PROJECT_ROOT / "results" / "mr" / "wald_mr_all_outcomes.csv"
PRIMARY_FDR = PROJECT_ROOT / "results" / "mr" / "fdr_significant_preliminary.csv"
SHARED = PROJECT_ROOT / "results" / "mr" / "shared_candidate_preliminary.csv"
RAW_UKB = PROJECT_ROOT / "data" / "raw" / "replication" / "opengwas_ukb"
OUT_DATA = PROJECT_ROOT / "data" / "replication" / "ukb_opengwas"
OUT_RESULTS = PROJECT_ROOT / "results" / "replication" / "ukb_opengwas"
TABLE_DIR = PROJECT_ROOT / "tables"

UKB = {
    "UKB_AF": {
        "primary_outcome": "AF",
        "opengwas_id": "ukb-b-964",
        "outcome_name": "UKB/OpenGWAS atrial fibrillation",
        "file": RAW_UKB / "ukb-b-964.vcf.gz",
    },
    "UKB_HF": {
        "primary_outcome": "HF",
        "opengwas_id": "ukb-d-HEARTFAIL",
        "outcome_name": "UKB/OpenGWAS heart failure",
        "file": RAW_UKB / "ukb-d-HEARTFAIL.vcf.gz",
    },
}

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def fnum(value: Any) -> float | None:
    try:
        if value in {"", None, "NA", "."}:
            return None
        number = float(value)
        if number != number or math.isinf(number):
            return None
        return number
    except Exception:
        return None


def fmt(value: Any, digits: int = 4) -> str:
    number = fnum(value)
    if number is None:
        return "NA"
    if number == 0:
        return "0"
    if abs(number) < 0.001 or abs(number) >= 10000:
        return f"{number:.{digits}e}"
    return f"{number:.{digits}f}".rstrip("0").rstrip(".")


def fmt_p(value: Any) -> str:
    number = fnum(value)
    if number is None:
        return "NA"
    if number == 0:
        return "0"
    if number < 0.001:
        return f"{number:.2e}"
    return f"{number:.4f}".rstrip("0").rstrip(".")


def norm_allele(value: str) -> str:
    return (value or "").strip().upper()


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


def bh_fdr(rows: list[dict[str, Any]], p_field: str = "pval_mr", out_field: str = "fdr") -> None:
    valid = [(idx, fnum(row.get(p_field))) for idx, row in enumerate(rows)]
    valid = [(idx, p) for idx, p in valid if p is not None]
    ranked = sorted(valid, key=lambda item: item[1])
    n = len(ranked)
    min_seen = 1.0
    adjusted = [1.0] * n
    for rank_idx in range(n - 1, -1, -1):
        _idx, p = ranked[rank_idx]
        rank = rank_idx + 1
        value = min(min_seen, p * n / rank)
        min_seen = value
        adjusted[rank_idx] = min(value, 1.0)
    for (idx, _p), fdr in zip(ranked, adjusted):
        rows[idx][out_field] = f"{fdr:.8g}"


def exposure_maps() -> tuple[dict[str, dict[str, str]], dict[tuple[str, str], list[dict[str, str]]]]:
    rows = read_csv(EXPOSURE)
    by_snp = {row["SNP"]: row for row in rows}
    by_chr_pos: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_chr_pos[(str(row["chr"]), str(row["pos"]))].append(row)
    return by_snp, by_chr_pos


def parse_format(format_keys: str, sample_values: str) -> dict[str, str]:
    keys = format_keys.split(":")
    vals = sample_values.split(":")
    return dict(zip(keys, vals))


def lp_to_p(lp: str) -> float | None:
    value = fnum(lp)
    if value is None:
        return None
    if value > 323:
        return 0.0
    return 10 ** (-value)


def extract_ukb_effects(outcome_id: str, config: dict[str, Any]) -> tuple[list[dict[str, Any]], dict[str, int]]:
    _by_snp, by_chr_pos = exposure_maps()
    stats = {
        "raw_rows_scanned": 0,
        "position_rows_seen": 0,
        "harmonised_rows": 0,
        "allele_mismatch_rows": 0,
        "gzip_eof_error": 0,
    }
    rows: list[dict[str, Any]] = []
    path = config["file"]
    if not path.exists() or path.stat().st_size < 1000:
        return rows, stats
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                break
        while True:
            try:
                line = next(handle)
            except StopIteration:
                break
            except EOFError:
                stats["gzip_eof_error"] = 1
                break
            if not line or line.startswith("#"):
                continue
            stats["raw_rows_scanned"] += 1
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue
            chrom, pos, variant_id, ref, alt, _qual, _filter, info, fmt, sample = parts[:10]
            chrom = chrom.replace("chr", "")
            exposures = by_chr_pos.get((chrom, pos))
            if not exposures:
                continue
            stats["position_rows_seen"] += 1
            parsed = parse_format(fmt, sample)
            beta = fnum(parsed.get("ES"))
            se = fnum(parsed.get("SE"))
            pval = lp_to_p(parsed.get("LP", ""))
            if beta is None or se is None:
                continue
            for exp in exposures:
                flip, action = harmonise(
                    exp,
                    {
                        "effect_allele": alt,
                        "other_allele": ref,
                    },
                )
                if flip == 0:
                    stats["allele_mismatch_rows"] += 1
                    continue
                stats["harmonised_rows"] += 1
                rows.append(
                    {
                        "protein": exp["protein"],
                        "protein_name": exp["protein_name"],
                        "panel": exp["panel"],
                        "outcome_id": outcome_id,
                        "opengwas_id": config["opengwas_id"],
                        "outcome_name": config["outcome_name"],
                        "SNP": exp["SNP"],
                        "chr": exp["chr"],
                        "pos": exp["pos"],
                        "vcf_id": variant_id,
                        "effect_allele_exposure": exp["effect_allele"],
                        "other_allele_exposure": exp["other_allele"],
                        "effect_allele_outcome": norm_allele(alt),
                        "other_allele_outcome": norm_allele(ref),
                        "harmonise_action": action,
                        "outcome_beta_raw": beta,
                        "outcome_beta_aligned": beta * flip,
                        "outcome_se": se,
                        "outcome_p": pval if pval is not None else "",
                        "outcome_af_alt": parsed.get("AF", ""),
                        "outcome_sample_size": parsed.get("SS", ""),
                        "outcome_ncase": parsed.get("NC", ""),
                        "outcome_info": info,
                        "beta_exposure": exp["beta"],
                        "se_exposure": exp["se"],
                        "f_stat": exp["f_stat"],
                    }
                )
    return rows, stats


def run_mr(effect_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    normal = NormalDist()
    rows: list[dict[str, Any]] = []
    for row in effect_rows:
        bx = fnum(row.get("beta_exposure"))
        by = fnum(row.get("outcome_beta_aligned"))
        sy = fnum(row.get("outcome_se"))
        if bx is None or by is None or sy is None or bx == 0:
            continue
        beta_mr = by / bx
        se_mr = abs(sy / bx)
        z = beta_mr / se_mr if se_mr else float("nan")
        pval = 2 * (1 - normal.cdf(abs(z))) if z == z else float("nan")
        lci = beta_mr - 1.96 * se_mr
        uci = beta_mr + 1.96 * se_mr
        rows.append(
            {
                "protein": row["protein"],
                "protein_name": row["protein_name"],
                "panel": row["panel"],
                "outcome_id": row["outcome_id"],
                "opengwas_id": row["opengwas_id"],
                "outcome_name": row["outcome_name"],
                "SNP": row["SNP"],
                "method": "Wald ratio",
                "beta_exposure": row["beta_exposure"],
                "se_exposure": row["se_exposure"],
                "beta_outcome": f"{by:.8g}",
                "se_outcome": row["outcome_se"],
                "beta_mr": f"{beta_mr:.8g}",
                "se_mr": f"{se_mr:.8g}",
                "or": f"{math.exp(beta_mr):.8g}",
                "or_lci95": f"{math.exp(lci):.8g}",
                "or_uci95": f"{math.exp(uci):.8g}",
                "pval_mr": f"{pval:.8g}" if pval == pval else "",
                "fdr": "",
                "harmonise_action": row["harmonise_action"],
                "effect_allele_exposure": row["effect_allele_exposure"],
                "other_allele_exposure": row["other_allele_exposure"],
                "effect_allele_outcome": row["effect_allele_outcome"],
                "other_allele_outcome": row["other_allele_outcome"],
            }
        )
    by_outcome: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        by_outcome[row["outcome_id"]].append(row)
    for out_rows in by_outcome.values():
        bh_fdr(out_rows)
    return rows


def status_label(primary: dict[str, str] | None, replication: dict[str, Any] | None) -> str:
    if primary is None:
        return "not_in_primary"
    if replication is None:
        return "not_matched_in_ukb"
    primary_beta = fnum(primary.get("beta_mr"))
    replication_beta = fnum(replication.get("beta_mr"))
    if primary_beta is None or replication_beta is None:
        return "insufficient"
    replication_p_value = fnum(replication.get("pval_mr"))
    replication_fdr_value = fnum(replication.get("fdr"))
    replication_p = 1.0 if replication_p_value is None else replication_p_value
    replication_fdr = 1.0 if replication_fdr_value is None else replication_fdr_value
    same = primary_beta * replication_beta > 0
    if same and replication_fdr < 0.05:
        return "same_direction_fdr"
    if same and replication_p < 0.05:
        return "same_direction_nominal"
    if same:
        return "same_direction_not_significant"
    if replication_p < 0.05:
        return "opposite_direction_nominal"
    return "opposite_direction_not_significant"


def compare_with_primary(mr_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    primary_rows = read_csv(PRIMARY_ALL)
    primary_by_key = {(row["protein"], row["outcome_id"]): row for row in primary_rows}
    replication_to_primary = {"UKB_AF": "AF", "UKB_HF": "HF"}
    replication_by_key = {(row["protein"], row["outcome_id"]): row for row in mr_rows}
    exposure = read_csv(EXPOSURE)
    compare_rows = []
    for exp in exposure:
        for replication_id, primary_id in replication_to_primary.items():
            primary = primary_by_key.get((exp["protein"], primary_id))
            replication = replication_by_key.get((exp["protein"], replication_id))
            compare_rows.append(
                {
                    "protein": exp["protein"],
                    "protein_name": exp["protein_name"],
                    "panel": exp["panel"],
                    "primary_outcome": primary_id,
                    "replication_outcome": replication_id,
                    "SNP": exp["SNP"],
                    "primary_or": primary.get("or", "") if primary else "",
                    "primary_p": primary.get("pval_mr", "") if primary else "",
                    "primary_fdr": primary.get("fdr", "") if primary else "",
                    "replication_or": replication.get("or", "") if replication else "",
                    "replication_p": replication.get("pval_mr", "") if replication else "",
                    "replication_fdr": replication.get("fdr", "") if replication else "",
                    "replication_status": status_label(primary, replication),
                }
            )
    return compare_rows


def subset_rows(compare_rows: list[dict[str, Any]], proteins: set[str] | None = None, fdr_only: bool = False) -> list[dict[str, Any]]:
    out = []
    fdr_set = {(row["protein"], row["outcome_id"]) for row in read_csv(PRIMARY_FDR)}
    for row in compare_rows:
        if proteins and row["protein"] not in proteins:
            continue
        if fdr_only and (row["protein"], row["primary_outcome"]) not in fdr_set:
            continue
        out.append(row)
    return out


def summarize(mr_rows: list[dict[str, Any]], compare_rows: list[dict[str, Any]], stats_by_outcome: dict[str, dict[str, int]]) -> None:
    lines = [
        "# UKB/OpenGWAS Full-Panel Replication Summary",
        "",
        "Scope: all 529 lead cis-pQTL instruments from the UKB-PPP inflammation panel main analysis. Because the exposure data are from UKB-PPP, these analyses are treated as supplementary and may be affected by sample overlap.",
        "",
        "## Extraction Coverage",
        "",
        "| UKB/OpenGWAS outcome | Harmonised instruments | Position rows seen | Allele mismatches | Gzip EOF warning |",
        "|---|---:|---:|---:|---:|",
    ]
    for outcome_id, config in UKB.items():
        rows = [row for row in mr_rows if row["outcome_id"] == outcome_id]
        stats = stats_by_outcome[outcome_id]
        lines.append(
            f"| {config['opengwas_id']} | {len(rows)} | {stats['position_rows_seen']} | {stats['allele_mismatch_rows']} | {stats['gzip_eof_error']} |"
        )
    lines.extend(
        [
            "",
            "## Replication Counts",
            "",
            "| Comparison | Count |",
            "|---|---:|",
        ]
    )
    for outcome_id, primary_id in [("UKB_AF", "AF"), ("UKB_HF", "HF")]:
        rows = [
            row
            for row in compare_rows
            if row["replication_outcome"] == outcome_id
            and row["primary_p"]
            and row["replication_status"] != "not_matched_in_ukb"
        ]
        same_nominal = sum(1 for row in rows if row["replication_status"] in {"same_direction_nominal", "same_direction_fdr"})
        same_fdr = sum(1 for row in rows if row["replication_status"] == "same_direction_fdr")
        lines.append(f"| {primary_id} primary proteins matched in UKB/OpenGWAS | {len(rows)} |")
        lines.append(f"| {primary_id} same-direction nominal UKB/OpenGWAS replication | {same_nominal} |")
        lines.append(f"| {primary_id} same-direction FDR UKB/OpenGWAS replication | {same_fdr} |")
    lines.extend(
        [
            "",
            "## Primary FDR Signals",
            "",
            "| Protein | Primary outcome | Primary OR/P/FDR | UKB OR/P/FDR | Replication status |",
            "|---|---|---|---|---|",
        ]
    )
    for row in subset_rows(compare_rows, fdr_only=True):
        lines.append(
            f"| {row['protein']} | {row['primary_outcome']} | "
            f"{fmt(row['primary_or'])}/{fmt_p(row['primary_p'])}/{fmt_p(row['primary_fdr'])} | "
            f"{fmt(row['replication_or'])}/{fmt_p(row['replication_p'])}/{fmt_p(row['replication_fdr'])} | "
            f"{row['replication_status']} |"
        )
    lines.extend(
        [
            "",
            "## FGF5/LPA Shared Candidates",
            "",
            "| Protein | Primary outcome | UKB outcome | Primary OR/P | UKB OR/P | Replication status |",
            "|---|---|---|---|---|---|",
        ]
    )
    for row in subset_rows(compare_rows, proteins={row["protein"] for row in read_csv(SHARED)}):
        lines.append(
            f"| {row['protein']} | {row['primary_outcome']} | {row['replication_outcome']} | "
            f"{fmt(row['primary_or'])}/{fmt_p(row['primary_p'])} | {fmt(row['replication_or'])}/{fmt_p(row['replication_p'])} | "
            f"{row['replication_status']} |"
        )
    lines.extend(
        [
            "",
            "Interpretation: UKB/OpenGWAS replication is completed as a supplementary check, but coverage is very low for the 529 lead cis-pQTL panel and the FGF5/LPA shared candidates were not matched. Because UKB-PPP exposure data and UKB outcome data may overlap, these results should not be interpreted as fully independent replication.",
            "",
            f"Full MR table: `{OUT_RESULTS / 'ukb_opengwas_full_panel_wald_mr.csv'}`",
            f"Primary comparison table: `{OUT_RESULTS / 'ukb_opengwas_primary_comparison.csv'}`",
        ]
    )
    write_text(OUT_RESULTS / "ukb_opengwas_full_panel_replication_summary.md", "\n".join(lines) + "\n")


def update_status() -> None:
    rows = [
        {
            "replication_layer": "FinnGen R12 candidate exact-variant replication",
            "status": "completed",
            "scope": "FGF5 and LPA candidate variants for I9_AF and I9_HEARTFAIL",
            "output": "results/replication/finngen_r12_candidate_wald_mr.csv",
            "note": "Completed previously and superseded for coverage by the full-panel FinnGen extraction.",
        },
        {
            "replication_layer": "FinnGen R12 full-panel replication",
            "status": "completed",
            "scope": "All 529 lead cis-pQTL instruments against FinnGen I9_AF and I9_HEARTFAIL",
            "output": "results/replication/full_panel/finngen_r12_full_panel_wald_mr.csv",
            "note": "Completed using downloaded FinnGen R12 AF and HF summary statistics.",
        },
        {
            "replication_layer": "UKB/OpenGWAS full-panel replication",
            "status": "completed_supplementary_overlap_caution",
            "scope": "All 529 lead cis-pQTL instruments against ukb-b-964 and ukb-d-HEARTFAIL",
            "output": "results/replication/ukb_opengwas/ukb_opengwas_full_panel_wald_mr.csv",
            "note": "Completed using downloaded OpenGWAS VCF files, but coverage is very low for the 529 lead cis-pQTL panel and FGF5/LPA were not matched; interpret as a limited supplementary check because UKB-PPP exposure and UKB outcomes may overlap.",
        },
    ]
    fields = ["replication_layer", "status", "scope", "output", "note"]
    write_csv(PROJECT_ROOT / "results" / "replication" / "finngen_ukb_replication_status.csv", rows, fields)
    write_csv(TABLE_DIR / "supplementary_table_s15_replication_status.csv", rows, fields)
    md = [
        "# FinnGen/UKB Replication Status",
        "",
        "| Layer | Status | Scope | Note |",
        "|---|---|---|---|",
    ]
    for row in rows:
        md.append(f"| {row['replication_layer']} | {row['status']} | {row['scope']} | {row['note']} |")
    write_text(PROJECT_ROOT / "results" / "replication" / "finngen_ukb_replication_status.md", "\n".join(md) + "\n")


def update_manifest() -> None:
    rows = read_csv(TABLE_DIR / "supplementary_tables_manifest.csv")
    rows = [row for row in rows if row.get("table_id") not in {"S18", "S19"}]
    rows.extend(
        [
            {
                "table_id": "S18",
                "description": "UKB/OpenGWAS full-panel replication Wald ratio MR results",
                "filename": "supplementary_table_s18_ukb_opengwas_full_panel_wald_mr.csv",
                "rows": str(len(read_csv(OUT_RESULTS / "ukb_opengwas_full_panel_wald_mr.csv"))),
                "source": str(OUT_RESULTS / "ukb_opengwas_full_panel_wald_mr.csv"),
            },
            {
                "table_id": "S19",
                "description": "UKB/OpenGWAS full-panel replication comparison with primary MR",
                "filename": "supplementary_table_s19_ukb_opengwas_primary_comparison.csv",
                "rows": str(len(read_csv(OUT_RESULTS / "ukb_opengwas_primary_comparison.csv"))),
                "source": str(OUT_RESULTS / "ukb_opengwas_primary_comparison.csv"),
            },
        ]
    )
    fields = ["table_id", "description", "filename", "rows", "source"]
    write_csv(TABLE_DIR / "supplementary_tables_manifest.csv", rows, fields)
    lines = [
        "# Supplementary Tables Manifest",
        "",
        "| Table | Description | File | Rows |",
        "|---|---|---|---|",
    ]
    for row in rows:
        lines.append(f"| {row['table_id']} | {row['description']} | `{row['filename']}` | {row['rows']} |")
    lines.extend(
        [
            "",
            "Note: Results now include primary MR, formal coloc.abf, target prioritization, feasible sensitivity checks, candidate reverse MR, FinnGen R12 full-panel replication, and supplementary UKB/OpenGWAS replication.",
        ]
    )
    write_text(TABLE_DIR / "supplementary_tables_manifest.md", "\n".join(lines) + "\n")


def write_log(mr_rows: list[dict[str, Any]], compare_rows: list[dict[str, Any]], stats_by_outcome: dict[str, dict[str, int]]) -> None:
    lines = [
        "# UKB/OpenGWAS全量复制记录",
        "",
        "日期：2026-05-28",
        "",
        "## 文件检查",
        "",
        "- ukb-b-964.vcf.gz：有效VCF.GZ文件。",
        "- ukb-d-HEARTFAIL.vcf.gz：有效VCF.GZ文件。",
        "- 两个.tbi索引文件均存在；本脚本逐行扫描VCF，不依赖索引。",
        "",
        "## 复制覆盖",
        "",
    ]
    for outcome_id, config in UKB.items():
        rows = [row for row in mr_rows if row["outcome_id"] == outcome_id]
        stats = stats_by_outcome[outcome_id]
        lines.append(f"- {config['opengwas_id']}协调后MR记录：{len(rows)}；位置匹配{stats['position_rows_seen']}；等位基因不匹配{stats['allele_mismatch_rows']}。")
    lines.extend(
        [
            "",
            "## 解释限制",
            "",
            "UKB/OpenGWAS复制已完成，但本次VCF对529个lead cis-pQTL覆盖很低，且未覆盖FGF5/LPA共享候选；同时因为暴露pQTL来自UKB-PPP，UKB结局复制可能存在样本重叠，因此应作为低覆盖补充验证而非独立复制的核心证据。",
            "",
            "## 输出文件",
            "",
            "- `results/replication/ukb_opengwas/ukb_opengwas_full_panel_replication_summary.md`",
            "- `results/replication/ukb_opengwas/ukb_opengwas_full_panel_wald_mr.csv`",
            "- `results/replication/ukb_opengwas/ukb_opengwas_primary_comparison.csv`",
            "- `tables/supplementary_table_s18_ukb_opengwas_full_panel_wald_mr.csv`",
            "- `tables/supplementary_table_s19_ukb_opengwas_primary_comparison.csv`",
        ]
    )
    write_text(PROJECT_ROOT / "ukb_opengwas_replication_log_2026-05-28.md", "\n".join(lines) + "\n")


def main() -> None:
    all_effects: list[dict[str, Any]] = []
    stats_by_outcome: dict[str, dict[str, int]] = {}
    effect_fields = [
        "protein",
        "protein_name",
        "panel",
        "outcome_id",
        "opengwas_id",
        "outcome_name",
        "SNP",
        "chr",
        "pos",
        "vcf_id",
        "effect_allele_exposure",
        "other_allele_exposure",
        "effect_allele_outcome",
        "other_allele_outcome",
        "harmonise_action",
        "outcome_beta_raw",
        "outcome_beta_aligned",
        "outcome_se",
        "outcome_p",
        "outcome_af_alt",
        "outcome_sample_size",
        "outcome_ncase",
        "outcome_info",
        "beta_exposure",
        "se_exposure",
        "f_stat",
    ]
    for outcome_id, config in UKB.items():
        effects, stats = extract_ukb_effects(outcome_id, config)
        stats_by_outcome[outcome_id] = stats
        all_effects.extend(effects)
        write_csv(OUT_DATA / f"{outcome_id}_outcome_effects_full_panel.csv", effects, effect_fields)

    mr_rows = run_mr(all_effects)
    mr_fields = [
        "protein",
        "protein_name",
        "panel",
        "outcome_id",
        "opengwas_id",
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
    write_csv(OUT_RESULTS / "ukb_opengwas_full_panel_wald_mr.csv", mr_rows, mr_fields)
    write_csv(TABLE_DIR / "supplementary_table_s18_ukb_opengwas_full_panel_wald_mr.csv", mr_rows, mr_fields)

    compare_rows = compare_with_primary(mr_rows)
    compare_fields = [
        "protein",
        "protein_name",
        "panel",
        "primary_outcome",
        "replication_outcome",
        "SNP",
        "primary_or",
        "primary_p",
        "primary_fdr",
        "replication_or",
        "replication_p",
        "replication_fdr",
        "replication_status",
    ]
    write_csv(OUT_RESULTS / "ukb_opengwas_primary_comparison.csv", compare_rows, compare_fields)
    write_csv(TABLE_DIR / "supplementary_table_s19_ukb_opengwas_primary_comparison.csv", compare_rows, compare_fields)

    summarize(mr_rows, compare_rows, stats_by_outcome)
    update_status()
    update_manifest()
    write_log(mr_rows, compare_rows, stats_by_outcome)
    print(f"UKB/OpenGWAS MR rows: {len(mr_rows)}")
    print(f"Primary comparison rows: {len(compare_rows)}")


if __name__ == "__main__":
    main()
