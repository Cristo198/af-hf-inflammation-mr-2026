"""Run full-panel FinnGen R12 replication for all lead cis-pQTL instruments."""

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
RAW_REPLICATION = PROJECT_ROOT / "data" / "raw" / "replication"
OUT_DATA = PROJECT_ROOT / "data" / "replication" / "full_panel"
OUT_RESULTS = PROJECT_ROOT / "results" / "replication" / "full_panel"
TABLE_DIR = PROJECT_ROOT / "tables"

FINNGEN = {
    "FG_AF": {
        "phenocode": "I9_AF",
        "outcome_name": "FinnGen R12 atrial fibrillation and flutter",
        "file": RAW_REPLICATION / "finngen_r12" / "finngen_R12_I9_AF.gz",
        "primary_outcome_id": "AF",
    },
    "FG_HF": {
        "phenocode": "I9_HEARTFAIL",
        "outcome_name": "FinnGen R12 heart failure strict",
        "file": RAW_REPLICATION / "finngen_r12" / "finngen_R12_I9_HEARTFAIL.gz",
        "primary_outcome_id": "HF",
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
        if value in {"", None, "NA"}:
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


def candidate_score(row: dict[str, str], exposures: list[dict[str, str]]) -> tuple[dict[str, str] | None, int, str]:
    best: tuple[dict[str, str] | None, int, str] = (None, 0, "no_match")
    for exp in exposures:
        flip, action = harmonise(
            exp,
            {
                "effect_allele": row["alt"],
                "other_allele": row["ref"],
            },
        )
        if flip:
            return exp, flip, action
        rsids = {value.strip() for value in row.get("rsids", "").split(",") if value.strip()}
        if exp["SNP"] in rsids and best[0] is None:
            best = (exp, 0, "rsid_position_match_allele_mismatch")
    return best


def extract_finngen_effects(outcome_id: str, config: dict[str, Any]) -> tuple[list[dict[str, Any]], dict[str, int]]:
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
        reader = csv.DictReader(handle, delimiter="\t")
        while True:
            try:
                row = next(reader)
            except StopIteration:
                break
            except EOFError:
                stats["gzip_eof_error"] = 1
                break
            stats["raw_rows_scanned"] += 1
            chrom = str(row.get("#chrom", row.get("chrom", ""))).replace("chr", "")
            pos = str(row.get("pos", ""))
            exposures = by_chr_pos.get((chrom, pos))
            if not exposures:
                continue
            stats["position_rows_seen"] += 1
            exp, flip, action = candidate_score(row, exposures)
            if exp is None:
                continue
            if flip == 0:
                stats["allele_mismatch_rows"] += 1
                continue
            beta = fnum(row.get("beta"))
            se = fnum(row.get("sebeta"))
            pval = fnum(row.get("pval"))
            if beta is None or se is None:
                continue
            stats["harmonised_rows"] += 1
            rows.append(
                {
                    "protein": exp["protein"],
                    "protein_name": exp["protein_name"],
                    "panel": exp["panel"],
                    "outcome_id": outcome_id,
                    "phenocode": config["phenocode"],
                    "outcome_name": config["outcome_name"],
                    "SNP": exp["SNP"],
                    "chr": exp["chr"],
                    "pos": exp["pos"],
                    "finngen_chr": chrom,
                    "finngen_pos": pos,
                    "finngen_rsids": row.get("rsids", ""),
                    "effect_allele_exposure": exp["effect_allele"],
                    "other_allele_exposure": exp["other_allele"],
                    "effect_allele_outcome": norm_allele(row.get("alt", "")),
                    "other_allele_outcome": norm_allele(row.get("ref", "")),
                    "harmonise_action": action,
                    "outcome_beta_raw": beta,
                    "outcome_beta_aligned": beta * flip,
                    "outcome_se": se,
                    "outcome_p": pval if pval is not None else "",
                    "outcome_af_alt": row.get("af_alt", ""),
                    "outcome_af_alt_cases": row.get("af_alt_cases", ""),
                    "outcome_af_alt_controls": row.get("af_alt_controls", ""),
                    "nearest_genes": row.get("nearest_genes", ""),
                    "beta_exposure": exp["beta"],
                    "se_exposure": exp["se"],
                    "f_stat": exp["f_stat"],
                }
            )
    return rows, stats


def run_mr(effect_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    normal = NormalDist()
    mr_rows: list[dict[str, Any]] = []
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
        mr_rows.append(
            {
                "protein": row["protein"],
                "protein_name": row["protein_name"],
                "panel": row["panel"],
                "outcome_id": row["outcome_id"],
                "phenocode": row["phenocode"],
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
    for row in mr_rows:
        by_outcome[row["outcome_id"]].append(row)
    for rows in by_outcome.values():
        bh_fdr(rows)
    return mr_rows


def status_label(primary: dict[str, str] | None, replication: dict[str, Any] | None) -> str:
    if primary is None:
        return "not_in_primary"
    if replication is None:
        return "not_matched_in_finngen"
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
    replication_to_primary = {"FG_AF": "AF", "FG_HF": "HF"}
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
        "# FinnGen R12 Full-Panel Replication Summary",
        "",
        "Scope: all 529 lead cis-pQTL instruments from the UKB-PPP inflammation panel main analysis.",
        "",
        "## Extraction Coverage",
        "",
        "| FinnGen outcome | Harmonised instruments | Position rows seen | Allele mismatches | Gzip EOF warning |",
        "|---|---:|---:|---:|---:|",
    ]
    for outcome_id, config in FINNGEN.items():
        rows = [row for row in mr_rows if row["outcome_id"] == outcome_id]
        stats = stats_by_outcome[outcome_id]
        lines.append(
            f"| {config['phenocode']} | {len(rows)} | {stats['position_rows_seen']} | {stats['allele_mismatch_rows']} | {stats.get('gzip_eof_error', 0)} |"
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
    for outcome_id, primary_id in [("FG_AF", "AF"), ("FG_HF", "HF")]:
        rows = [
            row
            for row in compare_rows
            if row["replication_outcome"] == outcome_id
            and row["primary_p"]
            and row["replication_status"] != "not_matched_in_finngen"
        ]
        same_nominal = sum(1 for row in rows if row["replication_status"] in {"same_direction_nominal", "same_direction_fdr"})
        same_fdr = sum(1 for row in rows if row["replication_status"] == "same_direction_fdr")
        lines.append(f"| {primary_id} primary proteins matched in FinnGen | {len(rows)} |")
        lines.append(f"| {primary_id} same-direction nominal FinnGen replication | {same_nominal} |")
        lines.append(f"| {primary_id} same-direction FDR FinnGen replication | {same_fdr} |")

    fdr_compare = subset_rows(compare_rows, fdr_only=True)
    shared_compare = subset_rows(compare_rows, proteins={row["protein"] for row in read_csv(SHARED)})
    lines.extend(
        [
            "",
            "## Primary FDR Signals",
            "",
            "| Protein | Primary outcome | Primary OR/P/FDR | FinnGen OR/P/FDR | Replication status |",
            "|---|---|---|---|---|",
        ]
    )
    for row in fdr_compare:
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
            "| Protein | Primary outcome | FinnGen outcome | Primary OR/P | FinnGen OR/P | Replication status |",
            "|---|---|---|---|---|---|",
        ]
    )
    for row in shared_compare:
        lines.append(
            f"| {row['protein']} | {row['primary_outcome']} | {row['replication_outcome']} | "
            f"{fmt(row['primary_or'])}/{fmt_p(row['primary_p'])} | {fmt(row['replication_or'])}/{fmt_p(row['replication_p'])} | "
            f"{row['replication_status']} |"
        )

    af_eof = stats_by_outcome.get("FG_AF", {}).get("gzip_eof_error", 0)
    hf_eof = stats_by_outcome.get("FG_HF", {}).get("gzip_eof_error", 0)
    if af_eof or hf_eof:
        interpretation = (
            "FinnGen replication was completed for the readable parts of the downloaded files. "
            "At least one gzip EOF warning remains, so affected outcomes should be treated as partial until re-downloaded cleanly. "
            "UKB/OpenGWAS replication remains unavailable because the local UKB files are HTML error pages rather than valid VCF.GZ files."
        )
    else:
        interpretation = (
            "FinnGen R12 AF and HF full-panel replication is now complete for the 529 lead cis-pQTL instruments. "
            "UKB/OpenGWAS replication remains unavailable because the local UKB files are HTML error pages rather than valid VCF.GZ files."
        )
    lines.extend(
        [
            "",
            f"Interpretation: {interpretation}",
            "",
            f"Full MR table: `{OUT_RESULTS / 'finngen_r12_full_panel_wald_mr.csv'}`",
            f"Primary comparison table: `{OUT_RESULTS / 'finngen_r12_primary_comparison.csv'}`",
        ]
    )
    write_text(OUT_RESULTS / "finngen_r12_full_panel_replication_summary.md", "\n".join(lines) + "\n")


def update_status(stats_by_outcome: dict[str, dict[str, int]]) -> None:
    af_eof = stats_by_outcome.get("FG_AF", {}).get("gzip_eof_error", 0)
    hf_eof = stats_by_outcome.get("FG_HF", {}).get("gzip_eof_error", 0)
    finngen_status = "completed" if not af_eof and not hf_eof else "partial_due_to_gzip_eof"
    finngen_note = (
        "Completed using downloaded FinnGen R12 AF and HF summary statistics."
        if finngen_status == "completed"
        else "At least one FinnGen summary-statistics file raised a gzip EOF warning; affected results should be treated as partial until re-downloaded cleanly."
    )
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
            "status": finngen_status,
            "scope": "All 529 lead cis-pQTL instruments against FinnGen I9_AF and I9_HEARTFAIL",
            "output": "results/replication/full_panel/finngen_r12_full_panel_wald_mr.csv",
            "note": finngen_note,
        },
        {
            "replication_layer": "UKB outcome replication",
            "status": "blocked_invalid_download",
            "scope": "AF/HF replication in UKB/OpenGWAS outcome files",
            "output": "",
            "note": "Downloaded UKB/OpenGWAS files are HTML pages, not valid VCF.GZ files. Also consider UKB-PPP exposure overlap before interpreting UKB as independent replication.",
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
    manifest = TABLE_DIR / "supplementary_tables_manifest.csv"
    rows = read_csv(manifest)
    rows = [row for row in rows if row.get("table_id") not in {"S16", "S17"}]
    rows.extend(
        [
            {
                "table_id": "S16",
                "description": "FinnGen R12 full-panel replication Wald ratio MR results",
                "filename": "supplementary_table_s16_finngen_r12_full_panel_wald_mr.csv",
                "rows": str(len(read_csv(OUT_RESULTS / "finngen_r12_full_panel_wald_mr.csv"))),
                "source": str(OUT_RESULTS / "finngen_r12_full_panel_wald_mr.csv"),
            },
            {
                "table_id": "S17",
                "description": "FinnGen R12 full-panel replication comparison with primary MR",
                "filename": "supplementary_table_s17_finngen_r12_primary_comparison.csv",
                "rows": str(len(read_csv(OUT_RESULTS / "finngen_r12_primary_comparison.csv"))),
                "source": str(OUT_RESULTS / "finngen_r12_primary_comparison.csv"),
            },
        ]
    )
    fields = ["table_id", "description", "filename", "rows", "source"]
    write_csv(manifest, rows, fields)
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
            "Note: Results now include primary MR, formal coloc.abf, target prioritization, feasible sensitivity checks, candidate reverse MR, and FinnGen R12 full-panel replication.",
        ]
    )
    write_text(TABLE_DIR / "supplementary_tables_manifest.md", "\n".join(lines) + "\n")


def write_log(mr_rows: list[dict[str, Any]], compare_rows: list[dict[str, Any]], stats_by_outcome: dict[str, dict[str, int]]) -> None:
    af_rows = [row for row in mr_rows if row["outcome_id"] == "FG_AF"]
    hf_rows = [row for row in mr_rows if row["outcome_id"] == "FG_HF"]
    fdr_compare = subset_rows(compare_rows, fdr_only=True)
    lines = [
        "# FinnGen R12全量复制记录",
        "",
        "日期：2026-05-27",
        "",
        "## 文件检查",
        "",
        "- FinnGen R12 I9_AF summary statistics：" + ("有效gzip文件，已用于全量复制。" if stats_by_outcome.get("FG_AF", {}).get("gzip_eof_error", 0) == 0 else "可读但存在gzip结尾异常，已输出当前可恢复的部分复制结果。"),
        "- FinnGen R12 I9_HEARTFAIL summary statistics：" + ("有效gzip文件，已用于全量复制。" if stats_by_outcome.get("FG_HF", {}).get("gzip_eof_error", 0) == 0 else "可读但存在gzip结尾异常，已输出当前可恢复的部分复制结果。"),
        "- UKB/OpenGWAS AF和HF文件：当前下载结果为HTML页面，不是有效VCF.GZ文件，UKB复制暂不执行。",
        f"- FinnGen gzip结尾异常警告：AF={stats_by_outcome.get('FG_AF', {}).get('gzip_eof_error', 0)}；HF={stats_by_outcome.get('FG_HF', {}).get('gzip_eof_error', 0)}。",
        "",
        "## 复制覆盖",
        "",
        f"- FinnGen AF协调后MR记录：{len(af_rows)}。",
        f"- FinnGen HF协调后MR记录：{len(hf_rows)}。",
        "",
        "## 主FDR信号复制状态",
        "",
        "| Protein | Outcome | FinnGen status |",
        "|---|---|---|",
    ]
    for row in fdr_compare:
        lines.append(f"| {row['protein']} | {row['primary_outcome']} | {row['replication_status']} |")
    lines.extend(
        [
            "",
            "## 输出文件",
            "",
            "- `results/replication/full_panel/finngen_r12_full_panel_replication_summary.md`",
            "- `results/replication/full_panel/finngen_r12_full_panel_wald_mr.csv`",
            "- `results/replication/full_panel/finngen_r12_primary_comparison.csv`",
            "- `tables/supplementary_table_s16_finngen_r12_full_panel_wald_mr.csv`",
            "- `tables/supplementary_table_s17_finngen_r12_primary_comparison.csv`",
        ]
    )
    write_text(PROJECT_ROOT / "finngen_full_replication_log_2026-05-27.md", "\n".join(lines) + "\n")


def main() -> None:
    all_effects: list[dict[str, Any]] = []
    stats_by_outcome: dict[str, dict[str, int]] = {}
    for outcome_id, config in FINNGEN.items():
        effects, stats = extract_finngen_effects(outcome_id, config)
        stats_by_outcome[outcome_id] = stats
        all_effects.extend(effects)
        write_csv(
            OUT_DATA / f"{outcome_id}_outcome_effects_full_panel.csv",
            effects,
            [
                "protein",
                "protein_name",
                "panel",
                "outcome_id",
                "phenocode",
                "outcome_name",
                "SNP",
                "chr",
                "pos",
                "finngen_chr",
                "finngen_pos",
                "finngen_rsids",
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
                "outcome_af_alt_cases",
                "outcome_af_alt_controls",
                "nearest_genes",
                "beta_exposure",
                "se_exposure",
                "f_stat",
            ],
        )

    mr_rows = run_mr(all_effects)
    mr_fields = [
        "protein",
        "protein_name",
        "panel",
        "outcome_id",
        "phenocode",
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
    write_csv(OUT_RESULTS / "finngen_r12_full_panel_wald_mr.csv", mr_rows, mr_fields)
    write_csv(TABLE_DIR / "supplementary_table_s16_finngen_r12_full_panel_wald_mr.csv", mr_rows, mr_fields)

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
    write_csv(OUT_RESULTS / "finngen_r12_primary_comparison.csv", compare_rows, compare_fields)
    write_csv(TABLE_DIR / "supplementary_table_s17_finngen_r12_primary_comparison.csv", compare_rows, compare_fields)

    summarize(mr_rows, compare_rows, stats_by_outcome)
    update_status(stats_by_outcome)
    update_manifest()
    write_log(mr_rows, compare_rows, stats_by_outcome)
    print(f"FinnGen MR rows: {len(mr_rows)}")
    print(f"Primary comparison rows: {len(compare_rows)}")


if __name__ == "__main__":
    main()
