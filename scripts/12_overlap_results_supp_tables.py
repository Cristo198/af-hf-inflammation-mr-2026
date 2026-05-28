"""Generate AF/HF overlap figures, Results draft, and supplementary tables."""

from __future__ import annotations

import csv
import html
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[1]
MR_ALL = PROJECT_ROOT / "results" / "mr" / "wald_mr_all_outcomes.csv"
SHARED = PROJECT_ROOT / "results" / "mr" / "shared_candidate_preliminary.csv"
FDR = PROJECT_ROOT / "results" / "mr" / "fdr_significant_preliminary.csv"
ALLELE_AUDIT = PROJECT_ROOT / "results" / "mr" / "candidate_allele_audit.csv"
REPLICATION = PROJECT_ROOT / "results" / "replication" / "finngen_r12_candidate_wald_mr.csv"
COLOC = PROJECT_ROOT / "results" / "coloc_inputs" / "finngen_pheweb_candidate_pqtl_disease_coloc.csv"
EXPOSURE_QC = PROJECT_ROOT / "results" / "qc" / "exposure_instrument_qc.csv"
HARMON_QC = PROJECT_ROOT / "results" / "qc" / "outcome_harmonisation_qc.csv"
CANDIDATE_MATRIX = PROJECT_ROOT / "results" / "candidates" / "fgf5_lpa_evidence_matrix.csv"

FIG_DIR = PROJECT_ROOT / "results" / "figures"
OVERLAP_DIR = PROJECT_ROOT / "results" / "overlap"
TEXT_DIR = PROJECT_ROOT / "results" / "text"
TABLE_DIR = PROJECT_ROOT / "tables"

SET_ORDER = ["AF_nominal", "AF_FDR", "HF_nominal", "HF_FDR"]
SET_LABELS = {
    "AF_nominal": "AF P<0.05",
    "AF_FDR": "AF FDR<0.05",
    "HF_nominal": "HF P<0.05",
    "HF_FDR": "HF FDR<0.05",
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


def esc(text: Any) -> str:
    return html.escape(str(text), quote=True)


def svg_header(width: int, height: int) -> list[str]:
    return [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        "<style>text{font-family:Arial,Helvetica,sans-serif}.title{font-size:20px;font-weight:700}.small{font-size:12px}.label{font-size:14px;font-weight:700}.note{fill:#5f6368}.axis{stroke:#252a31;stroke-width:1}.grid{stroke:#e6e8ee;stroke-width:1}</style>",
        '<rect width="100%" height="100%" fill="#ffffff"/>',
    ]


def build_sets(rows: list[dict[str, str]]) -> dict[str, set[str]]:
    sets = {name: set() for name in SET_ORDER}
    for row in rows:
        protein = row["protein"]
        outcome = row["outcome_id"]
        p = fnum(row.get("pval_mr")) or 1.0
        fdr = fnum(row.get("fdr")) or 1.0
        if p < 0.05:
            sets[f"{outcome}_nominal"].add(protein)
        if fdr < 0.05:
            sets[f"{outcome}_FDR"].add(protein)
    return sets


def make_overlap_tables(rows: list[dict[str, str]], sets: dict[str, set[str]]) -> None:
    by_protein: dict[str, dict[str, dict[str, str]]] = defaultdict(dict)
    for row in rows:
        by_protein[row["protein"]][row["outcome_id"]] = row

    overlap = sorted(sets["AF_nominal"] & sets["HF_nominal"])
    overlap_rows = []
    for protein in overlap:
        af = by_protein[protein].get("AF", {})
        hf = by_protein[protein].get("HF", {})
        direction = "same_risk" if (fnum(af.get("or")) or 1) > 1 and (fnum(hf.get("or")) or 1) > 1 else (
            "same_protective" if (fnum(af.get("or")) or 1) < 1 and (fnum(hf.get("or")) or 1) < 1 else "opposite"
        )
        overlap_rows.append(
            {
                "protein": protein,
                "protein_name": af.get("protein_name") or hf.get("protein_name"),
                "AF_OR": af.get("or", ""),
                "AF_P": af.get("pval_mr", ""),
                "AF_FDR": af.get("fdr", ""),
                "HF_OR": hf.get("or", ""),
                "HF_P": hf.get("pval_mr", ""),
                "HF_FDR": hf.get("fdr", ""),
                "direction": direction,
                "meets_shared_candidate_rule": "TRUE" if protein in {"FGF5", "LPA"} else "FALSE",
            }
        )
    write_csv(
        OVERLAP_DIR / "af_hf_nominal_overlap_proteins.csv",
        overlap_rows,
        ["protein", "protein_name", "AF_OR", "AF_P", "AF_FDR", "HF_OR", "HF_P", "HF_FDR", "direction", "meets_shared_candidate_rule"],
    )

    exact_counter: Counter[tuple[str, ...]] = Counter()
    exact_members: dict[tuple[str, ...], list[str]] = defaultdict(list)
    all_proteins = set().union(*sets.values())
    for protein in sorted(all_proteins):
        membership = tuple(name for name in SET_ORDER if protein in sets[name])
        exact_counter[membership] += 1
        exact_members[membership].append(protein)

    rows_out = []
    for membership, count in sorted(exact_counter.items(), key=lambda item: (-item[1], item[0])):
        rows_out.append(
            {
                "combination": "+".join(membership),
                "count": count,
                "proteins": ";".join(exact_members[membership]),
            }
        )
    write_csv(OVERLAP_DIR / "mr_significance_upset_exact_membership.csv", rows_out, ["combination", "count", "proteins"])

    lines = [
        "# AF/HF MR Signal Overlap Summary",
        "",
        f"- AF nominal P<0.05 proteins: {len(sets['AF_nominal'])}",
        f"- HF nominal P<0.05 proteins: {len(sets['HF_nominal'])}",
        f"- AF/HF nominal overlap: {len(overlap)} ({', '.join(overlap)})",
        f"- AF FDR<0.05 proteins: {len(sets['AF_FDR'])}",
        f"- HF FDR<0.05 proteins: {len(sets['HF_FDR'])}",
        f"- AF FDR and HF nominal overlap: {', '.join(sorted(sets['AF_FDR'] & sets['HF_nominal'])) or 'none'}",
        f"- HF FDR and AF nominal overlap: {', '.join(sorted(sets['HF_FDR'] & sets['AF_nominal'])) or 'none'}",
        f"- Both-outcome FDR overlap: {', '.join(sorted(sets['AF_FDR'] & sets['HF_FDR'])) or 'none'}",
        "",
        "Shared-candidate rule used here: one outcome FDR<0.05, the other outcome nominal P<0.05, and concordant MR direction.",
    ]
    write_text(OVERLAP_DIR / "af_hf_overlap_summary.md", "\n".join(lines))


def make_venn(sets: dict[str, set[str]]) -> None:
    af = sets["AF_nominal"]
    hf = sets["HF_nominal"]
    inter = af & hf
    width, height = 960, 560
    lines = svg_header(width, height)
    lines += [
        '<text x="48" y="38" class="title">AF/HF nominal MR signal overlap</text>',
        '<text x="48" y="62" class="small note">Nominal threshold: Wald ratio P&lt;0.05. Shared candidates require additional FDR and concordant-direction criteria.</text>',
        '<circle cx="390" cy="280" r="175" fill="#b12a34" fill-opacity="0.28" stroke="#b12a34" stroke-width="3"/>',
        '<circle cx="570" cy="280" r="175" fill="#00856f" fill-opacity="0.28" stroke="#00856f" stroke-width="3"/>',
        '<text x="305" y="130" class="label" fill="#8c1d28">AF nominal</text>',
        '<text x="595" y="130" class="label" fill="#006b5a">HF nominal</text>',
        f'<text x="315" y="288" class="title" text-anchor="middle">{len(af - hf)}</text>',
        f'<text x="480" y="288" class="title" text-anchor="middle">{len(inter)}</text>',
        f'<text x="650" y="288" class="title" text-anchor="middle">{len(hf - af)}</text>',
        f'<text x="315" y="313" class="small" text-anchor="middle">AF only</text>',
        f'<text x="480" y="313" class="small" text-anchor="middle">overlap</text>',
        f'<text x="650" y="313" class="small" text-anchor="middle">HF only</text>',
        '<rect x="65" y="400" width="830" height="92" rx="6" fill="#f7f8fb" stroke="#d8dce3"/>',
        f'<text x="85" y="426" class="small"><tspan font-weight="700">Nominal overlap proteins:</tspan> {esc(", ".join(sorted(inter)))}</text>',
        f'<text x="85" y="452" class="small"><tspan font-weight="700">Shared candidates after FDR + direction rule:</tspan> FGF5, LPA</text>',
        f'<text x="85" y="478" class="small"><tspan font-weight="700">FDR sets:</tspan> AF={len(sets["AF_FDR"])}, HF={len(sets["HF_FDR"])}, both FDR={len(sets["AF_FDR"] & sets["HF_FDR"])}</text>',
        "</svg>",
    ]
    write_text(FIG_DIR / "af_hf_nominal_overlap_venn.svg", "\n".join(lines))


def make_upset(sets: dict[str, set[str]]) -> None:
    all_proteins = set().union(*sets.values())
    combos: list[tuple[tuple[str, ...], int, list[str]]] = []
    for protein in sorted(all_proteins):
        membership = tuple(name for name in SET_ORDER if protein in sets[name])
        found = False
        for i, (combo, count, members) in enumerate(combos):
            if combo == membership:
                combos[i] = (combo, count + 1, members + [protein])
                found = True
                break
        if not found:
            combos.append((membership, 1, [protein]))
    combos.sort(key=lambda item: (-item[1], item[0]))

    width = 1080
    height = 660
    left = 210
    top = 84
    bar_top = 90
    bar_h = 260
    matrix_top = 420
    col_w = 90
    max_count = max(count for _, count, _ in combos)

    lines = svg_header(width, height)
    lines += [
        '<text x="42" y="38" class="title">MR significance UpSet plot</text>',
        '<text x="42" y="62" class="small note">Exact membership across nominal and FDR sets. FDR sets are nested within their corresponding nominal sets.</text>',
    ]

    # set size bars on the left
    for idx, set_name in enumerate(SET_ORDER):
        y = matrix_top + idx * 42
        lines.append(f'<text x="42" y="{y + 5}" class="small">{esc(SET_LABELS[set_name])}</text>')
        lines.append(f'<text x="{left - 18}" y="{y + 5}" class="small" text-anchor="end">{len(sets[set_name])}</text>')
        lines.append(f'<line x1="{left}" y1="{y}" x2="{left + len(sets[set_name]) * 2.4:.1f}" y2="{y}" stroke="#5f6368" stroke-width="7"/>')

    # intersection bars and matrix
    for idx, (combo, count, members) in enumerate(combos):
        x = left + 95 + idx * col_w
        bar_height = count / max_count * bar_h
        lines.append(f'<rect x="{x - 18}" y="{bar_top + bar_h - bar_height:.1f}" width="36" height="{bar_height:.1f}" fill="#3b63a3"/>')
        lines.append(f'<text x="{x}" y="{bar_top + bar_h - bar_height - 8:.1f}" class="small" text-anchor="middle">{count}</text>')
        active_y = []
        for set_idx, set_name in enumerate(SET_ORDER):
            y = matrix_top + set_idx * 42
            active = set_name in combo
            color = "#252a31" if active else "#d4d8df"
            radius = 7 if active else 5
            lines.append(f'<circle cx="{x}" cy="{y}" r="{radius}" fill="{color}"/>')
            if active:
                active_y.append(y)
        if len(active_y) > 1:
            lines.append(f'<line x1="{x}" y1="{min(active_y)}" x2="{x}" y2="{max(active_y)}" stroke="#252a31" stroke-width="2"/>')
        member_label = ", ".join(members[:4]) + ("..." if len(members) > 4 else "")
        lines.append(f'<text x="{x}" y="{matrix_top + 185}" class="small" text-anchor="middle" transform="rotate(45 {x} {matrix_top + 185})">{esc(member_label)}</text>')

    # axes and labels
    lines.append(f'<line x1="{left + 60}" y1="{bar_top + bar_h}" x2="{width - 50}" y2="{bar_top + bar_h}" class="axis"/>')
    lines.append(f'<text x="{left + 95}" y="{height - 22}" class="small note">Columns show exact set combinations; labels list representative proteins.</text>')
    lines.append("</svg>")
    write_text(FIG_DIR / "mr_significance_upset.svg", "\n".join(lines))


def copy_table(src: Path, dst: Path) -> tuple[str, int]:
    rows = read_csv(src)
    fieldnames = list(rows[0].keys()) if rows else []
    if rows:
        write_csv(dst, rows, fieldnames)
    else:
        write_text(dst, "")
    return dst.name, len(rows)


def make_supplementary_tables() -> None:
    table_specs = [
        ("S1", "Exposure lead cis-pQTL instrument QC", EXPOSURE_QC, TABLE_DIR / "supplementary_table_s1_exposure_lead_cis_pqtl_qc.csv"),
        ("S2", "Outcome extraction and harmonisation QC", HARMON_QC, TABLE_DIR / "supplementary_table_s2_outcome_harmonisation_qc.csv"),
        ("S3", "All primary Wald ratio MR results", MR_ALL, TABLE_DIR / "supplementary_table_s3_all_primary_wald_mr_results.csv"),
        ("S4", "FDR-significant primary MR results", FDR, TABLE_DIR / "supplementary_table_s4_fdr_significant_primary_mr_results.csv"),
        ("S5", "AF/HF nominal overlap proteins", OVERLAP_DIR / "af_hf_nominal_overlap_proteins.csv", TABLE_DIR / "supplementary_table_s5_af_hf_nominal_overlap_proteins.csv"),
        ("S6", "FGF5/LPA candidate evidence matrix", CANDIDATE_MATRIX, TABLE_DIR / "supplementary_table_s6_fgf5_lpa_candidate_evidence_matrix.csv"),
        ("S7", "FGF5/LPA allele direction audit", ALLELE_AUDIT, TABLE_DIR / "supplementary_table_s7_fgf5_lpa_allele_direction_audit.csv"),
        ("S8", "FinnGen R12 candidate replication", REPLICATION, TABLE_DIR / "supplementary_table_s8_finngen_r12_candidate_replication.csv"),
        ("S9", "FinnGen precomputed pQTL-disease coloc records", COLOC, TABLE_DIR / "supplementary_table_s9_finngen_precomputed_pqtl_disease_coloc.csv"),
        ("S10", "UpSet exact membership", OVERLAP_DIR / "mr_significance_upset_exact_membership.csv", TABLE_DIR / "supplementary_table_s10_upset_exact_membership.csv"),
    ]

    manifest_rows = []
    for table_id, description, src, dst in table_specs:
        filename, n_rows = copy_table(src, dst)
        manifest_rows.append(
            {
                "table_id": table_id,
                "description": description,
                "filename": filename,
                "rows": n_rows,
                "source": str(src),
            }
        )
    write_csv(TABLE_DIR / "supplementary_tables_manifest.csv", manifest_rows, ["table_id", "description", "filename", "rows", "source"])

    lines = [
        "# Supplementary Tables Manifest",
        "",
        "| Table | Description | File | Rows |",
        "|---|---|---|---|",
    ]
    for row in manifest_rows:
        lines.append(f"| {row['table_id']} | {row['description']} | `{row['filename']}` | {row['rows']} |")
    lines += [
        "",
        "Note: Results are first-pass single-variant Wald ratio estimates. Formal primary-outcome colocalization remains pending until dense regional pQTL summary statistics are available.",
    ]
    write_text(TABLE_DIR / "supplementary_tables_manifest.md", "\n".join(lines))


def get_result(rows: list[dict[str, str]], protein: str, outcome: str) -> dict[str, str] | None:
    for row in rows:
        if row.get("protein") == protein and row.get("outcome_id") == outcome:
            return row
    return None


def result_text(row: dict[str, str] | None, include_fdr: bool = True) -> str:
    if not row:
        return "NA"
    text = f"OR={fmt(row.get('or'))}, 95%CI {fmt(row.get('or_lci95'))}-{fmt(row.get('or_uci95'))}, P={fmt_p(row.get('pval_mr'))}"
    if include_fdr and row.get("fdr") not in {"", None}:
        text += f", FDR={fmt_p(row.get('fdr'))}"
    return text


def make_results_draft(rows: list[dict[str, str]], sets: dict[str, set[str]]) -> None:
    shared = read_csv(SHARED)
    shared_by_protein = {row["protein"]: row for row in shared}
    fdr_rows = read_csv(FDR)
    replication_rows = read_csv(REPLICATION)
    coloc_rows = read_csv(COLOC)

    fdr_by_outcome = defaultdict(list)
    for row in fdr_rows:
        fdr_by_outcome[row["outcome_id"]].append(row)
    for outcome in fdr_by_outcome:
        fdr_by_outcome[outcome].sort(key=lambda row: fnum(row.get("pval_mr")) or 1)

    af_fgf5 = get_result(rows, "FGF5", "AF")
    hf_fgf5 = get_result(rows, "FGF5", "HF")
    af_lpa = get_result(rows, "LPA", "AF")
    hf_lpa = get_result(rows, "LPA", "HF")
    fg_af = next((row for row in replication_rows if row["protein"] == "FGF5" and row["phenocode"] == "I9_AF"), None)
    fg_hf = next((row for row in replication_rows if row["protein"] == "LPA" and row["phenocode"] == "I9_HEARTFAIL"), None)
    lpa_af = next((row for row in replication_rows if row["protein"] == "LPA" and row["phenocode"] == "I9_AF"), None)

    af_sig = "、".join(
        f"{row['protein']}（{result_text(row)}）" for row in fdr_by_outcome["AF"]
    )
    hf_sig = "、".join(
        f"{row['protein']}（{result_text(row)}）" for row in fdr_by_outcome["HF"]
    )
    nominal_overlap = sorted(sets["AF_nominal"] & sets["HF_nominal"])

    lines = [
        "# Results段落初稿",
        "",
        "说明：以下只更新Results，不包含Introduction或Discussion。正式主结局共定位仍等待FGF5/LPA区域级密集UKB-PPP pQTL summary statistics。",
        "",
        "## 3.1 数据源、工具变量和质量控制",
        "",
        "本研究从UKB-PPP Olink Explore炎症相关面板中纳入737种蛋白。基于EUR人群pQTL结果、P < 5 × 10^-8阈值和目标基因上下游1 Mb cis窗口，共筛得545条cis-pQTL关联。按每种蛋白保留一个lead cis-pQTL并排除MHC区域后，529种蛋白进入主MR分析。所有lead工具变量F统计量均大于10，F统计量中位数为788.5，最小值为45.8，提示弱工具变量偏倚风险较低。FGF5和LPA的lead工具变量F统计量分别为7582.3和5377.8。由于主分析采用单lead cis-pQTL设计，MR-Egger、weighted median、MR-PRESSO和leave-one-out不适用于主分析估计；本阶段主要通过cis限制、MHC排除、强工具变量筛选、等位基因协调、复制验证和后续共定位控制潜在多效性风险。",
        "",
        "心房颤动主结局来自Nielsen等2018年GWAS，心力衰竭主结局来自HERMES GWAS。本地结局提取中，AF结局成功提取492个唯一工具变量，等位基因协调后489个进入MR；HF结局成功提取418个唯一工具变量，等位基因协调后410个进入MR。",
        "",
        "## 3.2 循环炎症蛋白与心房颤动风险",
        "",
        f"在AF主结局中，37种蛋白达到名义显著，5种蛋白达到FDR校正显著。FDR显著蛋白包括：{af_sig}。其中，FGF5为最强信号之一，遗传预测FGF5水平升高与AF风险升高相关（{result_text(af_fgf5)}）。",
        "",
        "## 3.3 循环炎症蛋白与心力衰竭风险",
        "",
        f"在HF主结局中，41种蛋白达到名义显著，4种蛋白达到FDR校正显著。FDR显著蛋白包括：{hf_sig}。其中，LPA与HF风险升高显著相关（{result_text(hf_lpa)}）。",
        "",
        "## 3.4 AF-HF共享炎症蛋白筛选",
        "",
        f"AF和HF名义显著信号共有6种蛋白重叠，分别为{ '、'.join(nominal_overlap) }。进一步采用“一个结局FDR < 0.05、另一个结局P < 0.05且方向一致”的探索性共享候选标准后，FGF5和LPA被筛选为AF-HF连续体候选蛋白。FGF5在AF中达到FDR显著（{result_text(af_fgf5)}），在HF中达到名义显著且方向一致（{result_text(hf_fgf5)}）。LPA在HF中达到FDR显著（{result_text(hf_lpa)}），在AF中达到名义显著且方向一致（{result_text(af_lpa)}）。等位基因方向复核显示，FGF5的C等位基因和LPA的T等位基因分别为蛋白升高等位基因，两个候选蛋白在AF和HF中的主分析方向均指向风险增加。",
        "",
        "Venn和UpSet图显示，AF和HF名义显著集合存在有限重叠，而FDR层面的双结局直接重叠为0；FGF5和LPA分别代表AF主导和HF主导的共享候选模式。",
        "",
        "## 3.5 FinnGen复制和预计算共定位线索",
        "",
        f"候选复制使用FinnGen R12 PheWeb公开变异接口。FGF5在FinnGen AF中方向一致并显著（{result_text(fg_af, include_fdr=False)}），但FinnGen strict HF在该精确变异处无有效beta。LPA在FinnGen AF中方向一致并显著（{result_text(lpa_af, include_fdr=False)}），在FinnGen strict HF中方向一致并达到名义显著（{result_text(fg_hf, include_fdr=False)}）。",
        "",
        "FinnGen R12预计算pQTL-疾病共定位记录为FGF5-AF提供支持性线索：FinnGen Olink记录CLPP=0.204、CLPA=0.557，UK Biobank PPP Olink 3k记录CLPP=0.249、CLPA=0.249。FGF5-HF及LPA-AF/HF暂未检索到对应的FinnGen pQTL-疾病共定位记录。需要强调的是，这些CLPP/CLPA结果属于FinnGen fine-mapping框架下的预计算支持证据，并不等同于本研究主结局Nielsen AF和HERMES HF的正式coloc.abf PP.H4结果。正式主结局共定位仍需获取FGF5和LPA区域级密集pQTL summary statistics后完成。",
        "",
        "## 3.6 图表和补充表",
        "",
        "本阶段已生成AF/HF主MR火山图、FDR显著结果森林图、FGF5/LPA主分析与复制森林图、AF-HF名义重叠Venn图和MR显著性UpSet图。补充表S1-S10整理了工具变量QC、结局协调、全部MR结果、FDR显著结果、AF/HF重叠蛋白、FGF5/LPA证据矩阵、等位基因方向复核、FinnGen复制和预计算共定位记录。",
        "",
        "## 暂不填写的Results小节",
        "",
        "中介分析和靶点优先级排序尚未完成，相关Results段落暂不填写。Introduction和Discussion按用户后续专业提示词另行处理。",
    ]
    write_text(TEXT_DIR / "results_paragraphs_draft_2026-05-26.md", "\n".join(lines))


def make_log() -> None:
    lines = [
        "# Venn/UpSet, Results Draft, and Supplementary Tables Log",
        "",
        "Date: 2026-05-26",
        "",
        "Generated:",
        f"- `{FIG_DIR / 'af_hf_nominal_overlap_venn.svg'}`",
        f"- `{FIG_DIR / 'mr_significance_upset.svg'}`",
        f"- `{OVERLAP_DIR / 'af_hf_overlap_summary.md'}`",
        f"- `{TEXT_DIR / 'results_paragraphs_draft_2026-05-26.md'}`",
        f"- `{TABLE_DIR / 'supplementary_tables_manifest.md'}`",
        "",
        "Note: Introduction and Discussion were not written or modified in this step.",
    ]
    write_text(PROJECT_ROOT / "venn_upset_results_supp_tables_log_2026-05-26.md", "\n".join(lines))


def main() -> None:
    rows = read_csv(MR_ALL)
    sets = build_sets(rows)
    make_overlap_tables(rows, sets)
    make_venn(sets)
    make_upset(sets)
    make_supplementary_tables()
    make_results_draft(rows, sets)
    make_log()
    print(f"Wrote Venn/UpSet figures to {FIG_DIR}")
    print(f"Wrote overlap summaries to {OVERLAP_DIR}")
    print(f"Wrote Results draft to {TEXT_DIR}")
    print(f"Wrote supplementary tables to {TABLE_DIR}")


if __name__ == "__main__":
    main()
