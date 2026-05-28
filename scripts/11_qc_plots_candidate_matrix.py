"""Generate instrument QC, MR figures, and FGF5/LPA evidence matrix."""

from __future__ import annotations

import csv
import html
import math
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[1]
EXPOSURE = PROJECT_ROOT / "data" / "processed" / "ukbppp_inflammation_cis_pqtl_instruments_lead_no_mhc.csv"
OUTCOME_DIR = PROJECT_ROOT / "data" / "outcomes"
MR_DIR = PROJECT_ROOT / "results" / "mr"
REPLICATION = PROJECT_ROOT / "results" / "replication" / "finngen_r12_candidate_wald_mr.csv"
COLOC = PROJECT_ROOT / "results" / "coloc_inputs" / "finngen_pheweb_candidate_pqtl_disease_coloc.csv"
AUDIT = MR_DIR / "candidate_allele_audit.csv"
SHARED = MR_DIR / "shared_candidate_preliminary.csv"

QC_DIR = PROJECT_ROOT / "results" / "qc"
FIG_DIR = PROJECT_ROOT / "results" / "figures"
CAND_DIR = PROJECT_ROOT / "results" / "candidates"
TABLE_DIR = PROJECT_ROOT / "tables"

CANDIDATES = ["FGF5", "LPA"]
OUTCOMES = {
    "AF": "Atrial fibrillation",
    "HF": "Heart failure",
}


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
        if value in {"", None}:
            return None
        number = float(value)
        if number != number:
            return None
        return number
    except Exception:
        return None


def fmt(value: Any, digits: int = 3) -> str:
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
    if number < 0.001:
        return f"{number:.2e}"
    return f"{number:.4f}".rstrip("0").rstrip(".")


def percentile(values: list[float], pct: float) -> float:
    if not values:
        return float("nan")
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    pos = (len(ordered) - 1) * pct / 100
    lower = math.floor(pos)
    upper = math.ceil(pos)
    if lower == upper:
        return ordered[int(pos)]
    weight = pos - lower
    return ordered[lower] * (1 - weight) + ordered[upper] * weight


def is_palindromic(row: dict[str, str]) -> bool:
    alleles = {row.get("effect_allele", "").upper(), row.get("other_allele", "").upper()}
    return alleles in [{"A", "T"}, {"C", "G"}]


def marker_for_gene(row: dict[str, str]) -> str:
    annotation = row.get("locus_annotation", "")
    protein = row.get("protein", "")
    if annotation in {"", "-", protein}:
        return "target_or_blank"
    return "nearest_locus_not_target_review"


def make_qc() -> None:
    exposure = read_csv(EXPOSURE)
    exposure_by_snp = {row["SNP"]: row for row in exposure}
    exposure_qc = []
    for row in exposure:
        f_stat = fnum(row.get("f_stat"))
        exposure_qc.append(
            {
                **row,
                "palindromic_alleles": "TRUE" if is_palindromic(row) else "FALSE",
                "weak_instrument_f_le_10": "TRUE" if f_stat is not None and f_stat <= 10 else "FALSE",
                "cis_gene_match_ok": "TRUE" if row.get("cis_gene_match") == row.get("protein") else "FALSE",
                "locus_annotation_flag": marker_for_gene(row),
                "candidate_flag": "TRUE" if row.get("protein") in CANDIDATES else "FALSE",
            }
        )

    qc_fields = list(exposure_qc[0].keys()) if exposure_qc else []
    write_csv(QC_DIR / "exposure_instrument_qc.csv", exposure_qc, qc_fields)
    write_csv(TABLE_DIR / "supplementary_table_exposure_instruments_qc.csv", exposure_qc, qc_fields)

    f_stats = [fnum(row.get("f_stat")) for row in exposure if fnum(row.get("f_stat")) is not None]
    f_stats = [x for x in f_stats if x is not None]
    panel_counts = Counter(row.get("panel") for row in exposure)
    locus_flags = Counter(row["locus_annotation_flag"] for row in exposure_qc)
    pal_count = sum(1 for row in exposure_qc if row["palindromic_alleles"] == "TRUE")
    weak_count = sum(1 for row in exposure_qc if row["weak_instrument_f_le_10"] == "TRUE")
    cis_ok_count = sum(1 for row in exposure_qc if row["cis_gene_match_ok"] == "TRUE")

    outcome_qc = []
    harmonise_lines = []
    for outcome_id, outcome_name in OUTCOMES.items():
        outcome_rows = read_csv(OUTCOME_DIR / f"{outcome_id}_outcome_effects_local.csv")
        mr_rows = read_csv(MR_DIR / f"{outcome_id}_wald_mr.csv")
        outcome_snps = {row["SNP"] for row in outcome_rows}
        mr_snps = {row["SNP"] for row in mr_rows}
        actions = Counter(row.get("harmonise_action", "") for row in mr_rows)
        fdr_sig = sum(1 for row in mr_rows if (fnum(row.get("fdr")) or 1) < 0.05)
        nominal = sum(1 for row in mr_rows if (fnum(row.get("pval_mr")) or 1) < 0.05)
        row = {
            "outcome_id": outcome_id,
            "outcome_name": outcome_name,
            "exposure_instruments": len(exposure_by_snp),
            "outcome_effects_extracted": len(outcome_snps),
            "harmonised_mr_rows": len(mr_rows),
            "missing_outcome_effects": len(exposure_by_snp) - len(outcome_snps),
            "extracted_but_not_harmonised": len(outcome_snps - mr_snps),
            "nominal_p_lt_0_05": nominal,
            "fdr_lt_0_05": fdr_sig,
            "aligned": actions.get("aligned", 0),
            "flipped": actions.get("flipped", 0),
            "aligned_complement": actions.get("aligned_complement", 0),
            "flipped_complement": actions.get("flipped_complement", 0),
        }
        outcome_qc.append(row)
        harmonise_lines.append(
            f"- {outcome_id}: extracted {row['outcome_effects_extracted']}/{row['exposure_instruments']} SNPs; "
            f"harmonised MR rows {row['harmonised_mr_rows']}; not harmonised after extraction {row['extracted_but_not_harmonised']}; "
            f"aligned={row['aligned']}, flipped={row['flipped']}, "
            f"aligned_complement={row['aligned_complement']}, flipped_complement={row['flipped_complement']}."
        )

    write_csv(QC_DIR / "outcome_harmonisation_qc.csv", outcome_qc, list(outcome_qc[0].keys()))
    write_csv(TABLE_DIR / "supplementary_table_outcome_harmonisation_qc.csv", outcome_qc, list(outcome_qc[0].keys()))

    candidate_rows = [row for row in exposure_qc if row.get("protein") in CANDIDATES]
    write_csv(QC_DIR / "candidate_instrument_qc.csv", candidate_rows, qc_fields)

    lines = [
        "# Instrument Quality Control Summary",
        "",
        "## Exposure instruments",
        "",
        f"- Lead cis-pQTL instruments after MHC exclusion: {len(exposure)}",
        f"- Panel counts: Inflammation={panel_counts.get('Inflammation', 0)}, Inflammation_II={panel_counts.get('Inflammation_II', 0)}",
        f"- Weak instruments by F <= 10: {weak_count}",
        f"- F statistic: min={fmt(min(f_stats))}, Q1={fmt(percentile(f_stats, 25))}, median={fmt(median(f_stats))}, Q3={fmt(percentile(f_stats, 75))}, max={fmt(max(f_stats))}",
        f"- cis gene match equals target protein gene: {cis_ok_count}/{len(exposure)}",
        f"- Palindromic allele instruments: {pal_count}/{len(exposure)}",
        f"- Locus annotation flags: target_or_blank={locus_flags.get('target_or_blank', 0)}, nearest_locus_not_target_review={locus_flags.get('nearest_locus_not_target_review', 0)}",
        "",
        "## Candidate instrument strength",
        "",
        "| Protein | SNP | Effect allele | Other allele | EAF | Beta | SE | F statistic | Locus flag |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for row in candidate_rows:
        lines.append(
            f"| {row['protein']} | {row['SNP']} | {row['effect_allele']} | {row['other_allele']} | "
            f"{fmt(row['eaf'])} | {fmt(row['beta'])} | {fmt(row['se'])} | {fmt(row['f_stat'])} | {row['locus_annotation_flag']} |"
        )
    lines += [
        "",
        "## Outcome extraction and harmonisation",
        "",
        *harmonise_lines,
        "",
        "## Pleiotropy note",
        "",
        "Because the main design uses one lead cis-pQTL per protein, MR-Egger, weighted median, MR-PRESSO and leave-one-out are not statistically applicable to the primary single-variant estimates. The practical pleiotropy controls at this stage are cis restriction, MHC exclusion, strong-instrument filtering, allele harmonisation audit, replication, and formal colocalization. Formal colocalization remains the key next check for distinguishing a shared causal variant from LD-driven association.",
        "",
        "Output tables:",
        f"- `{QC_DIR / 'exposure_instrument_qc.csv'}`",
        f"- `{QC_DIR / 'outcome_harmonisation_qc.csv'}`",
        f"- `{QC_DIR / 'candidate_instrument_qc.csv'}`",
    ]
    write_text(QC_DIR / "instrument_quality_summary.md", "\n".join(lines))


def svg_header(width: int, height: int) -> list[str]:
    return [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        "<style>text{font-family:Arial,Helvetica,sans-serif}.small{font-size:12px}.axis{stroke:#252a31;stroke-width:1}.grid{stroke:#e7e9ee;stroke-width:1}.note{fill:#666}.title{font-size:20px;font-weight:700}.label{font-size:11px}</style>",
        '<rect width="100%" height="100%" fill="#ffffff"/>',
    ]


def esc(text: Any) -> str:
    return html.escape(str(text), quote=True)


def make_volcano(outcome_id: str) -> None:
    rows = read_csv(MR_DIR / f"{outcome_id}_wald_mr.csv")
    data = []
    for row in rows:
        beta = fnum(row.get("beta_mr"))
        p = fnum(row.get("pval_mr"))
        fdr = fnum(row.get("fdr"))
        if beta is None or p is None or p <= 0:
            continue
        data.append({**row, "beta": beta, "neglogp": -math.log10(p), "fdr_num": fdr or 1.0})
    betas = [row["beta"] for row in data]
    central = max(abs(percentile(betas, 2.5)), abs(percentile(betas, 97.5)), 0.25)
    x_limit = min(max(central * 1.15, 0.35), 2.0)
    y_max = max(row["neglogp"] for row in data) * 1.08
    y_max = max(y_max, 5.0)

    width, height = 980, 680
    left, right, top, bottom = 90, 40, 70, 90
    plot_w = width - left - right
    plot_h = height - top - bottom

    def sx(beta: float) -> float:
        clipped = max(min(beta, x_limit), -x_limit)
        return left + (clipped + x_limit) / (2 * x_limit) * plot_w

    def sy(neglogp: float) -> float:
        return top + (1 - neglogp / y_max) * plot_h

    lines = svg_header(width, height)
    title = f"{outcome_id} MR volcano plot"
    lines.append(f'<text x="{left}" y="36" class="title">{esc(title)}</text>')
    lines.append(f'<text x="{left}" y="56" class="small note">x-axis: Wald ratio log(OR); y-axis: -log10(P). X-axis is winsorized at +/-{x_limit:.2f} for readability.</text>')

    for i in range(6):
        y_val = y_max * i / 5
        y = sy(y_val)
        lines.append(f'<line x1="{left}" y1="{y:.1f}" x2="{left + plot_w}" y2="{y:.1f}" class="grid"/>')
        lines.append(f'<text x="{left - 12}" y="{y + 4:.1f}" text-anchor="end" class="small">{y_val:.1f}</text>')
    for x_val in [-x_limit, -x_limit / 2, 0, x_limit / 2, x_limit]:
        x = sx(x_val)
        lines.append(f'<line x1="{x:.1f}" y1="{top}" x2="{x:.1f}" y2="{top + plot_h}" class="grid"/>')
        lines.append(f'<text x="{x:.1f}" y="{top + plot_h + 22}" text-anchor="middle" class="small">{x_val:.2f}</text>')
    lines.append(f'<line x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" y2="{top + plot_h}" class="axis"/>')
    lines.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_h}" class="axis"/>')
    sig_y = sy(-math.log10(0.05))
    lines.append(f'<line x1="{left}" y1="{sig_y:.1f}" x2="{left + plot_w}" y2="{sig_y:.1f}" stroke="#9aa0a6" stroke-dasharray="5 5"/>')
    lines.append(f'<line x1="{sx(0):.1f}" y1="{top}" x2="{sx(0):.1f}" y2="{top + plot_h}" stroke="#9aa0a6" stroke-dasharray="5 5"/>')

    for row in sorted(data, key=lambda r: r["fdr_num"] < 0.05):
        if row["fdr_num"] < 0.05 and row["beta"] > 0:
            color = "#b12a34"
        elif row["fdr_num"] < 0.05 and row["beta"] < 0:
            color = "#2f6fbb"
        else:
            color = "#9aa0a6"
        radius = 4.5 if row["fdr_num"] < 0.05 or row["protein"] in CANDIDATES else 3
        lines.append(
            f'<circle cx="{sx(row["beta"]):.1f}" cy="{sy(row["neglogp"]):.1f}" r="{radius}" fill="{color}" fill-opacity="0.82"/>'
        )

    label_rows = [row for row in data if row["fdr_num"] < 0.05 or row["protein"] in CANDIDATES]
    label_rows = sorted(label_rows, key=lambda r: r["neglogp"], reverse=True)
    for idx, row in enumerate(label_rows):
        x = sx(row["beta"])
        y = sy(row["neglogp"])
        dx = 8 if row["beta"] >= 0 else -8
        anchor = "start" if row["beta"] >= 0 else "end"
        dy = -6 + (idx % 3) * 10
        lines.append(f'<text x="{x + dx:.1f}" y="{y + dy:.1f}" text-anchor="{anchor}" class="label">{esc(row["protein"])}</text>')

    lines.append(f'<text x="{left + plot_w / 2}" y="{height - 32}" text-anchor="middle" class="small">Wald ratio log(OR)</text>')
    lines.append(f'<text x="22" y="{top + plot_h / 2}" text-anchor="middle" transform="rotate(-90 22 {top + plot_h / 2})" class="small">-log10(P)</text>')
    lines.append(f'<circle cx="{width - 265}" cy="34" r="5" fill="#b12a34"/><text x="{width - 252}" y="38" class="small">FDR&lt;0.05, risk-increasing</text>')
    lines.append(f'<circle cx="{width - 265}" cy="54" r="5" fill="#2f6fbb"/><text x="{width - 252}" y="58" class="small">FDR&lt;0.05, protective</text>')
    lines.append("</svg>")
    write_text(FIG_DIR / f"{outcome_id.lower()}_wald_mr_volcano.svg", "\n".join(lines))


def forest_ticks(low: float, high: float) -> list[float]:
    candidates = [0.5, 0.67, 0.75, 0.9, 1.0, 1.1, 1.25, 1.5, 2.0]
    ticks = [tick for tick in candidates if low <= tick <= high]
    if 1.0 not in ticks:
        ticks.append(1.0)
    return sorted(set(ticks))


def make_forest(rows: list[dict[str, Any]], path: Path, title: str) -> None:
    rows = [row for row in rows if fnum(row.get("or")) and fnum(row.get("or_lci95")) and fnum(row.get("or_uci95"))]
    if not rows:
        return
    low = min(fnum(row["or_lci95"]) for row in rows if fnum(row["or_lci95"]) is not None)
    high = max(fnum(row["or_uci95"]) for row in rows if fnum(row["or_uci95"]) is not None)
    low = min(low * 0.92, 0.99)
    high = max(high * 1.08, 1.01)
    log_low = math.log(low)
    log_high = math.log(high)
    width = 1120
    row_h = 38
    top = 92
    height = top + 58 + row_h * len(rows)
    left_text = 28
    plot_left = 390
    plot_w = 430
    right_text = 840

    def sx(or_value: float) -> float:
        return plot_left + (math.log(or_value) - log_low) / (log_high - log_low) * plot_w

    lines = svg_header(width, height)
    lines.append(f'<text x="{left_text}" y="36" class="title">{esc(title)}</text>')
    lines.append(f'<text x="{left_text}" y="62" class="small note">Squares show OR; horizontal lines show 95% CI. Vertical line marks OR=1.</text>')
    lines.append(f'<text x="{left_text}" y="{top - 18}" class="small" font-weight="700">Protein / outcome</text>')
    lines.append(f'<text x="{right_text}" y="{top - 18}" class="small" font-weight="700">OR (95% CI), P/FDR</text>')
    for tick in forest_ticks(low, high):
        x = sx(tick)
        lines.append(f'<line x1="{x:.1f}" y1="{top - 8}" x2="{x:.1f}" y2="{height - 44}" class="grid"/>')
        lines.append(f'<text x="{x:.1f}" y="{height - 18}" text-anchor="middle" class="small">{tick:g}</text>')
    lines.append(f'<line x1="{sx(1.0):.1f}" y1="{top - 8}" x2="{sx(1.0):.1f}" y2="{height - 44}" stroke="#252a31" stroke-dasharray="4 4"/>')

    for idx, row in enumerate(rows):
        y = top + idx * row_h
        if idx % 2 == 0:
            lines.append(f'<rect x="18" y="{y - 20}" width="{width - 36}" height="{row_h}" fill="#f8f9fb"/>')
        or_value = fnum(row["or"]) or 1
        lci = fnum(row["or_lci95"]) or or_value
        uci = fnum(row["or_uci95"]) or or_value
        color = row.get("color", "#b12a34" if or_value >= 1 else "#2f6fbb")
        lines.append(f'<line x1="{sx(lci):.1f}" y1="{y:.1f}" x2="{sx(uci):.1f}" y2="{y:.1f}" stroke="{color}" stroke-width="2"/>')
        lines.append(f'<rect x="{sx(or_value) - 5:.1f}" y="{y - 5:.1f}" width="10" height="10" fill="{color}"/>')
        lines.append(f'<text x="{left_text}" y="{y + 4}" class="small">{esc(row["label"])}</text>')
        stat = f"{fmt(or_value)} ({fmt(lci)}-{fmt(uci)}), P={fmt_p(row.get('pval_mr'))}"
        if row.get("fdr") not in {"", None}:
            stat += f", FDR={fmt_p(row.get('fdr'))}"
        lines.append(f'<text x="{right_text}" y="{y + 4}" class="small">{esc(stat)}</text>')
    lines.append(f'<text x="{plot_left + plot_w / 2}" y="{height - 2}" text-anchor="middle" class="small">Odds ratio</text>')
    lines.append("</svg>")
    write_text(path, "\n".join(lines))


def make_figures() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    for outcome_id in OUTCOMES:
        make_volcano(outcome_id)

    fdr_rows = []
    for outcome_id in OUTCOMES:
        for row in read_csv(MR_DIR / f"{outcome_id}_wald_mr.csv"):
            if (fnum(row.get("fdr")) or 1) < 0.05:
                fdr_rows.append(
                    {
                        **row,
                        "label": f"{row['protein']} ({outcome_id})",
                        "color": "#b12a34" if (fnum(row.get("or")) or 1) >= 1 else "#2f6fbb",
                    }
                )
    fdr_rows.sort(key=lambda row: (row["outcome_id"], fnum(row.get("fdr")) or 1))
    make_forest(fdr_rows, FIG_DIR / "fdr_significant_wald_mr_forest.svg", "FDR-significant primary MR results")

    candidate_rows = []
    primary = read_csv(MR_DIR / "wald_mr_all_outcomes.csv")
    for row in primary:
        if row.get("protein") in CANDIDATES:
            candidate_rows.append(
                {
                    **row,
                    "label": f"{row['protein']} primary {row['outcome_id']}",
                    "color": "#8a3ffc" if row.get("protein") == "FGF5" else "#00856f",
                }
            )
    for row in read_csv(REPLICATION):
        candidate_rows.append(
            {
                **row,
                "label": f"{row['protein']} FinnGen {row['phenocode']}",
                "fdr": "",
                "color": "#8a3ffc" if row.get("protein") == "FGF5" else "#00856f",
            }
        )
    order = {"FGF5": 0, "LPA": 1}
    candidate_rows.sort(key=lambda row: (order.get(row.get("protein"), 9), row.get("label", "")))
    make_forest(candidate_rows, FIG_DIR / "fgf5_lpa_primary_replication_forest.svg", "FGF5/LPA primary and FinnGen replication")


def result_string(row: dict[str, str] | None, include_fdr: bool = True) -> str:
    if not row:
        return "NA"
    text = f"OR {fmt(row.get('or'))} ({fmt(row.get('or_lci95'))}-{fmt(row.get('or_uci95'))}), P={fmt_p(row.get('pval_mr'))}"
    if include_fdr and row.get("fdr") not in {"", None}:
        text += f", FDR={fmt_p(row.get('fdr'))}"
    return text


def make_evidence_matrix() -> None:
    shared = {row["protein"]: row for row in read_csv(SHARED)}
    primary = defaultdict(dict)
    for row in read_csv(MR_DIR / "wald_mr_all_outcomes.csv"):
        if row.get("protein") in CANDIDATES:
            primary[row["protein"]][row["outcome_id"]] = row
    replication = defaultdict(dict)
    for row in read_csv(REPLICATION):
        replication[row["protein"]][row["phenocode"]] = row
    audit = defaultdict(dict)
    for row in read_csv(AUDIT):
        audit[row["protein"]][row["outcome_id"]] = row
    coloc = defaultdict(list)
    for row in read_csv(COLOC):
        coloc[row["protein"]].append(row)

    rows = []
    for protein in CANDIDATES:
        af_audit = audit[protein].get("AF")
        hf_audit = audit[protein].get("HF")
        allele = af_audit.get("exposure_effect_allele", "") if af_audit else ""
        audit_text = (
            f"{allele} allele increases {protein}; AF {af_audit.get('harmonise_action') if af_audit else 'NA'}, "
            f"HF {hf_audit.get('harmonise_action') if hf_audit else 'NA'}; both primary outcomes point to higher risk."
        )
        coloc_rows = coloc.get(protein, [])
        if coloc_rows:
            coloc_text = "; ".join(
                f"{row['finnGen_disease']} {row['pqtl_source_display']} CLPP={fmt(row['clpp'])}, CLPA={fmt(row['clpa'])}"
                for row in coloc_rows
            )
        else:
            coloc_text = "No FinnGen I9_AF/I9_HEARTFAIL pQTL-disease coloc record detected."

        if protein == "FGF5":
            interpretation = "AF-forward evidence is strong; HF signal is primary nominal and needs exact replication or formal coloc."
        else:
            interpretation = "HF-forward evidence is strong; AF replication is directionally consistent; formal coloc remains important at the LPA locus."

        row = {
            "protein": protein,
            "protein_name": primary[protein].get("AF", primary[protein].get("HF", {})).get("protein_name", ""),
            "lead_snp": shared.get(protein, {}).get("SNP", ""),
            "primary_af": result_string(primary[protein].get("AF")),
            "primary_hf": result_string(primary[protein].get("HF")),
            "shared_rule": shared.get(protein, {}).get("shared_rule", ""),
            "allele_direction_audit": audit_text,
            "finngen_af": result_string(replication[protein].get("I9_AF"), include_fdr=False),
            "finngen_hf": result_string(replication[protein].get("I9_HEARTFAIL"), include_fdr=False),
            "finnGen_precomputed_coloc": coloc_text,
            "formal_primary_coloc_status": "Pending dense regional UKB-PPP pQTL summary statistics.",
            "current_interpretation": interpretation,
        }
        rows.append(row)

    fields = [
        "protein",
        "protein_name",
        "lead_snp",
        "primary_af",
        "primary_hf",
        "shared_rule",
        "allele_direction_audit",
        "finngen_af",
        "finngen_hf",
        "finnGen_precomputed_coloc",
        "formal_primary_coloc_status",
        "current_interpretation",
    ]
    write_csv(CAND_DIR / "fgf5_lpa_evidence_matrix.csv", rows, fields)
    write_csv(TABLE_DIR / "table_candidate_fgf5_lpa_evidence_matrix.csv", rows, fields)

    lines = [
        "# FGF5/LPA Evidence Matrix",
        "",
        "| Protein | Lead SNP | Primary AF | Primary HF | FinnGen AF | FinnGen HF | Coloc status | Current interpretation |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for row in rows:
        lines.append(
            f"| {row['protein']} | {row['lead_snp']} | {row['primary_af']} | {row['primary_hf']} | "
            f"{row['finngen_af']} | {row['finngen_hf']} | {row['finnGen_precomputed_coloc']} Formal primary coloc pending. | "
            f"{row['current_interpretation']} |"
        )
    lines += [
        "",
        "Allele audit:",
        "- FGF5: C allele raises FGF5; primary AF and HF point to higher risk after harmonisation.",
        "- LPA: T allele raises apolipoprotein(a); primary AF and HF point to higher risk after harmonisation.",
        "",
        f"CSV table: `{CAND_DIR / 'fgf5_lpa_evidence_matrix.csv'}`",
    ]
    write_text(CAND_DIR / "fgf5_lpa_evidence_matrix.md", "\n".join(lines))


def make_manifest() -> None:
    lines = [
        "# QC, Figure, and Candidate Matrix Outputs",
        "",
        "Generated outputs:",
        "",
        "## Quality control",
        "",
        f"- `{QC_DIR / 'instrument_quality_summary.md'}`",
        f"- `{QC_DIR / 'exposure_instrument_qc.csv'}`",
        f"- `{QC_DIR / 'outcome_harmonisation_qc.csv'}`",
        f"- `{QC_DIR / 'candidate_instrument_qc.csv'}`",
        "",
        "## Figures",
        "",
        f"- `{FIG_DIR / 'af_wald_mr_volcano.svg'}`",
        f"- `{FIG_DIR / 'hf_wald_mr_volcano.svg'}`",
        f"- `{FIG_DIR / 'fdr_significant_wald_mr_forest.svg'}`",
        f"- `{FIG_DIR / 'fgf5_lpa_primary_replication_forest.svg'}`",
        "",
        "## Candidate matrix",
        "",
        f"- `{CAND_DIR / 'fgf5_lpa_evidence_matrix.md'}`",
        f"- `{CAND_DIR / 'fgf5_lpa_evidence_matrix.csv'}`",
    ]
    write_text(PROJECT_ROOT / "qc_figures_candidate_matrix_log_2026-05-26.md", "\n".join(lines))


def main() -> None:
    make_qc()
    make_figures()
    make_evidence_matrix()
    make_manifest()
    print(f"Wrote QC outputs to {QC_DIR}")
    print(f"Wrote figures to {FIG_DIR}")
    print(f"Wrote candidate matrix to {CAND_DIR}")


if __name__ == "__main__":
    main()
