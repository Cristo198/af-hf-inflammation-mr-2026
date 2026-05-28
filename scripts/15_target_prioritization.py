"""Build target-prioritization scorecard for FGF5/LPA candidates."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[1]

SHARED = PROJECT_ROOT / "results" / "mr" / "shared_candidate_preliminary.csv"
QC = PROJECT_ROOT / "results" / "qc" / "candidate_instrument_qc.csv"
REPLICATION = PROJECT_ROOT / "results" / "replication" / "finngen_r12_candidate_wald_mr.csv"
COLOC = PROJECT_ROOT / "results" / "coloc" / "formal_coloc_abf_summary.csv"
MEDIATION = PROJECT_ROOT / "results" / "mediation" / "af_mediated_effects_fgf5_lpa.csv"
CANDIDATE_MATRIX = PROJECT_ROOT / "results" / "candidates" / "fgf5_lpa_evidence_matrix.csv"
CANDIDATE_MATRIX_MD = PROJECT_ROOT / "results" / "candidates" / "fgf5_lpa_evidence_matrix.md"
PRIORITY_DIR = PROJECT_ROOT / "results" / "prioritization"
TABLE_DIR = PROJECT_ROOT / "tables"


DRUGGABILITY_NOTES = {
    "FGF5": {
        "clinical_precedence": "No mature FGF5-directed cardiovascular drug-development signal found in the targeted public scan; prior AF-focused MR literature exists.",
        "external_evidence": "Prior AF MR publication: PMID 39059473. No phase 2/3 FGF5-directed cardiovascular outcome trial identified in the targeted scan.",
        "tractability_score": 0.5,
    },
    "LPA": {
        "clinical_precedence": "High: multiple Lp(a)-lowering RNA therapies targeting apolipoprotein(a)/LPA are in phase 3 cardiovascular outcome trials.",
        "external_evidence": "ClinicalTrials.gov: pelacarsen Lp(a)HORIZON NCT04023552; olpasiran OCEAN(a)-Outcomes NCT05581303; lepodisiran ACCLAIM-Lp(a) NCT06292013.",
        "tractability_score": 2.0,
    },
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
        if value in {"", None, "NA"}:
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


def fmt_or(row: dict[str, str], prefix: str) -> str:
    return (
        f"OR {fmt(row.get(prefix + '_OR'))} "
        f"({fmt(row.get(prefix + '_LCI95'))}-{fmt(row.get(prefix + '_UCI95'))}), "
        f"P={fmt_p(row.get(prefix + '_P'))}, FDR={fmt_p(row.get(prefix + '_FDR'))}"
    )


def score_mr(row: dict[str, str]) -> tuple[float, str]:
    af_fdr = fnum(row.get("AF_FDR")) or 1.0
    hf_fdr = fnum(row.get("HF_FDR")) or 1.0
    af_p = fnum(row.get("AF_P")) or 1.0
    hf_p = fnum(row.get("HF_P")) or 1.0
    if af_fdr < 0.05 and hf_fdr < 0.05:
        return 4.0, "Both AF and HF reach FDR significance."
    if af_fdr < 0.05 and hf_p < 0.05:
        return 3.0, "AF reaches FDR significance and HF is nominally significant."
    if hf_fdr < 0.05 and af_p < 0.05:
        return 3.0, "HF reaches FDR significance and AF is nominally significant."
    if af_p < 0.05 and hf_p < 0.05:
        return 2.0, "Both outcomes are nominally significant."
    return 0.0, "Shared MR evidence is not supported."


def score_replication(protein: str, rows: list[dict[str, str]]) -> tuple[float, str]:
    matches = [row for row in rows if row.get("protein") == protein and row.get("replication_direction") == "same_risk"]
    outcomes = {row.get("phenocode") for row in matches if fnum(row.get("pval_mr")) is not None and (fnum(row.get("pval_mr")) or 1.0) < 0.05}
    if {"I9_AF", "I9_HEARTFAIL"}.issubset(outcomes):
        return 2.0, "FinnGen supports the same risk-increasing direction for AF and strict HF."
    if "I9_AF" in outcomes:
        return 1.0, "FinnGen supports the same risk-increasing direction for AF; strict HF exact-variant replication is unavailable or not significant."
    if "I9_HEARTFAIL" in outcomes:
        return 1.0, "FinnGen supports the same risk-increasing direction for strict HF only."
    return 0.0, "No significant same-direction FinnGen candidate replication."


def score_coloc(protein: str, rows: list[dict[str, str]]) -> tuple[float, str, str]:
    matches = [row for row in rows if row.get("protein") == protein]
    pp = {row.get("outcome_id"): fnum(row.get("PP.H4")) for row in matches}
    h3 = {row.get("outcome_id"): fnum(row.get("PP.H3")) for row in matches}
    af = pp.get("AF") or 0.0
    hf = pp.get("HF") or 0.0
    if af >= 0.8 and hf >= 0.8:
        score = 4.0
        note = "Strong formal coloc for both AF and HF."
    elif af >= 0.8 or hf >= 0.8:
        score = 3.0
        outcome = "AF" if af >= 0.8 else "HF"
        other = "HF" if outcome == "AF" else "AF"
        note = f"Strong formal coloc for {outcome}; {other} does not show strong coloc."
    elif af >= 0.5 or hf >= 0.5:
        score = 1.0
        note = "Partial or intermediate formal coloc support only."
    else:
        score = 0.0
        note = "No formal coloc support under the single-causal-variant coloc.abf model."
    detail = f"AF PP.H4={fmt(af)}, HF PP.H4={fmt(hf)}, AF PP.H3={fmt(h3.get('AF'))}, HF PP.H3={fmt(h3.get('HF'))}"
    return score, note, detail


def score_mediation(protein: str, rows: list[dict[str, str]]) -> tuple[float, str]:
    row = next((item for item in rows if item.get("protein") == protein), None)
    if not row:
        return 0.0, "AF-mediated effect analysis not available."
    p = fnum(row.get("p_indirect")) or 1.0
    prop = fnum(row.get("proportion_mediated"))
    indirect_or = fmt(row.get("or_indirect"))
    lci = fmt(row.get("or_indirect_lci95"))
    uci = fmt(row.get("or_indirect_uci95"))
    prop_text = "NA" if prop is None else f"{prop * 100:.1f}%"
    note = (
        f"Exploratory AF-mediated indirect effect: OR {indirect_or} "
        f"({lci}-{uci}), P={fmt_p(p)}, mediated proportion {prop_text}."
    )
    same_direction = str(row.get("same_direction_indirect_total", "")).lower() == "true"
    if p < 0.05 and same_direction and prop is not None and prop >= 0.25:
        return 1.0, note
    if p < 0.05 and same_direction:
        return 0.5, note
    return 0.0, note


def penalty(protein: str) -> tuple[float, str]:
    if protein == "FGF5":
        return -0.5, "HF support remains nominal and HF formal coloc is not supportive; FGF5-AF has prior MR literature, reducing novelty for an AF-only claim."
    if protein == "LPA":
        return -1.5, "Formal coloc is not supportive for AF/HF and the LPA region shows strong distinct-signal/LD concern, especially for HF."
    return 0.0, ""


def class_from_score(score: float) -> str:
    if score >= 8.5:
        return "Tier 1 - 遗传证据高优先级"
    if score >= 6.5:
        return "Tier 2 - 中等优先级/转化可开发性突出"
    return "Tier 3 - 探索性优先级"


def build_scorecard() -> list[dict[str, Any]]:
    shared = read_csv(SHARED)
    qc = {row["protein"]: row for row in read_csv(QC)}
    replication = read_csv(REPLICATION)
    coloc = read_csv(COLOC)
    mediation = read_csv(MEDIATION)

    rows: list[dict[str, Any]] = []
    for row in shared:
        protein = row["protein"]
        mr_score, mr_note = score_mr(row)
        replication_score, replication_note = score_replication(protein, replication)
        coloc_score, coloc_note, coloc_detail = score_coloc(protein, coloc)
        mediation_score, mediation_note = score_mediation(protein, mediation)
        penalty_score, penalty_note = penalty(protein)
        qcrow = qc.get(protein, {})
        instrument_score = 1.0 if (fnum(qcrow.get("f_stat")) or 0.0) > 10 else 0.0
        direction_score = 1.0 if row.get("direction") == "same_risk" else 0.0
        tractability_score = DRUGGABILITY_NOTES[protein]["tractability_score"]
        total = (
            mr_score
            + direction_score
            + replication_score
            + coloc_score
            + mediation_score
            + instrument_score
            + tractability_score
            + penalty_score
        )
        genetic_score = mr_score + direction_score + replication_score + coloc_score + mediation_score + instrument_score + penalty_score
        rows.append(
            {
                "rank": 0,
                "protein": protein,
                "protein_name": row.get("protein_name", ""),
                "uniprot": qcrow.get("uniprot", ""),
                "panel": row.get("panel", ""),
                "lead_snp": row.get("SNP", ""),
                "effect_allele_raises_protein": qcrow.get("effect_allele", ""),
                "primary_af": fmt_or(row, "AF"),
                "primary_hf": fmt_or(row, "HF"),
                "shared_rule": row.get("shared_rule", ""),
                "instrument_f": fmt(qcrow.get("f_stat")),
                "mr_score_0_4": mr_score,
                "direction_score_0_1": direction_score,
                "replication_score_0_2": replication_score,
                "formal_coloc_score_0_4": coloc_score,
                "mediation_score_0_1": mediation_score,
                "instrument_score_0_1": instrument_score,
                "tractability_score_0_2": tractability_score,
                "penalty": penalty_score,
                "genetic_score_after_penalty": round(genetic_score, 2),
                "overall_score": round(total, 2),
                "priority_class": class_from_score(total),
                "mr_note": mr_note,
                "replication_note": replication_note,
                "formal_coloc_note": coloc_note,
                "formal_coloc_detail": coloc_detail,
                "mediation_note": mediation_note,
                "clinical_precedence": DRUGGABILITY_NOTES[protein]["clinical_precedence"],
                "external_evidence": DRUGGABILITY_NOTES[protein]["external_evidence"],
                "penalty_note": penalty_note,
                "recommended_interpretation": "",
            }
        )

    rows.sort(key=lambda item: (item["overall_score"], item["genetic_score_after_penalty"]), reverse=True)
    for idx, row in enumerate(rows, start=1):
        row["rank"] = idx
        if row["protein"] == "FGF5":
            row["recommended_interpretation"] = (
                "当前遗传证据优先级最高的候选：以AF为主导，主MR、FinnGen AF复制、等位基因方向和FGF5-AF正式共定位一致；"
                "HF证据方向一致但仅为名义显著，且未获正式共定位支持。"
            )
        elif row["protein"] == "LPA":
            row["recommended_interpretation"] = (
                "次级候选：HF主MR和FinnGen AF/HF复制方向一致，且Lp(a)靶向药物开发基础很强；"
                "但本研究正式共定位未支持AF/HF共享单一因果变异。"
            )
    return rows


def write_scorecard(rows: list[dict[str, Any]]) -> None:
    fieldnames = [
        "rank",
        "protein",
        "protein_name",
        "uniprot",
        "panel",
        "lead_snp",
        "effect_allele_raises_protein",
        "primary_af",
        "primary_hf",
        "shared_rule",
        "instrument_f",
        "mr_score_0_4",
        "direction_score_0_1",
        "replication_score_0_2",
                "formal_coloc_score_0_4",
                "mediation_score_0_1",
                "instrument_score_0_1",
        "tractability_score_0_2",
        "penalty",
        "genetic_score_after_penalty",
        "overall_score",
        "priority_class",
        "mr_note",
        "replication_note",
                "formal_coloc_note",
                "formal_coloc_detail",
                "mediation_note",
                "clinical_precedence",
        "external_evidence",
        "penalty_note",
        "recommended_interpretation",
    ]
    write_csv(PRIORITY_DIR / "target_priority_scorecard.csv", rows, fieldnames)
    write_csv(TABLE_DIR / "supplementary_table_s12_target_priority_scorecard.csv", rows, fieldnames)

    lines = [
        "# Target Prioritization Scorecard",
        "",
        "Scope: FGF5 and LPA, the two exploratory AF-HF shared candidates selected by the rule of one primary outcome FDR < 0.05, the other primary outcome P < 0.05, and consistent risk direction.",
        "",
        "Scoring framework: primary MR evidence (0-4), same-direction AF-HF signal (0-1), FinnGen replication (0-2), formal coloc.abf evidence (0-4), exploratory AF-mediated pathway support (0-1), instrument confidence (0-1), clinical tractability (0-2), and predefined evidence penalty.",
        "",
        "| Rank | Protein | Priority class | Overall score | Genetic score | Key interpretation |",
        "|---|---|---|---:|---:|---|",
    ]
    for row in rows:
        lines.append(
            f"| {row['rank']} | {row['protein']} | {row['priority_class']} | {row['overall_score']} | "
            f"{row['genetic_score_after_penalty']} | {row['recommended_interpretation']} |"
        )
    lines.extend(
        [
            "",
            "## Evidence Notes",
            "",
        ]
    )
    for row in rows:
        lines.extend(
            [
                f"### {row['protein']}",
                "",
                f"- Primary MR: {row['mr_note']} AF: {row['primary_af']}; HF: {row['primary_hf']}.",
                f"- Replication: {row['replication_note']}",
                f"- Formal coloc: {row['formal_coloc_note']} {row['formal_coloc_detail']}.",
                f"- Mediation: {row['mediation_note']}",
                f"- Clinical tractability: {row['clinical_precedence']}",
                f"- Caution: {row['penalty_note']}",
                "",
            ]
        )
    lines.extend(
        [
            "## External Sources Used For Tractability/Precedence",
            "",
            "- Open Targets target-prioritisation/tractability documentation: https://platform-docs.opentargets.org/web-interface/target-prioritisation and https://platform-docs.opentargets.org/target/tractability",
            "- FGF5 AF MR prior literature: https://pubmed.ncbi.nlm.nih.gov/39059473/",
            "- Pelacarsen Lp(a)HORIZON: https://clinicaltrials.gov/study/NCT04023552",
            "- Olpasiran OCEAN(a)-Outcomes: https://clinicaltrials.gov/study/NCT05581303",
            "- Lepodisiran ACCLAIM-Lp(a): https://clinicaltrials.gov/study/NCT06292013",
            "",
            f"CSV scorecard: `{PRIORITY_DIR / 'target_priority_scorecard.csv'}`",
        ]
    )
    write_text(PRIORITY_DIR / "target_priority_scorecard.md", "\n".join(lines) + "\n")


def update_candidate_matrix() -> None:
    rows = read_csv(CANDIDATE_MATRIX)
    if not rows:
        return
    mediation = {row["protein"]: row for row in read_csv(MEDIATION)}
    for row in rows:
        med = mediation.get(row["protein"], {})
        if med:
            prop = fnum(med.get("proportion_mediated"))
            prop_text = "NA" if prop is None else f"{prop * 100:.1f}%"
            row["af_mediation_status"] = (
                f"Indirect OR {fmt(med.get('or_indirect'))} "
                f"({fmt(med.get('or_indirect_lci95'))}-{fmt(med.get('or_indirect_uci95'))}), "
                f"P={fmt_p(med.get('p_indirect'))}, mediated proportion {prop_text}."
            )
        else:
            row["af_mediation_status"] = "NA"
        if row["protein"] == "FGF5":
            row["formal_primary_coloc_status"] = (
                "FGF5-AF strong formal coloc: PP.H4=0.9868, PP.H3=0.0132; "
                "FGF5-HF no strong coloc: PP.H4=0.0544, PP.H3=0.0682."
            )
            row["current_interpretation"] = (
                "AF-forward evidence is strongest: primary MR, FinnGen AF replication, allele audit, formal FGF5-AF coloc, and exploratory AF mediation align; "
                "HF remains nominal and not directly coloc-supported."
            )
        elif row["protein"] == "LPA":
            row["formal_primary_coloc_status"] = (
                "LPA-AF no strong coloc: PP.H4=0.0987, PP.H3=0.4054; "
                "LPA-HF no strong coloc with high distinct-signal support: PP.H4=0.0396, PP.H3=0.9604."
            )
            row["current_interpretation"] = (
                "HF-forward MR, FinnGen replication, and a smaller AF-mediated component are supportive, and LPA is highly tractable clinically; "
                "formal coloc currently argues against a shared single causal variant in the tested AF/HF loci."
            )
    fieldnames = list(rows[0].keys())
    write_csv(CANDIDATE_MATRIX, rows, fieldnames)
    write_csv(TABLE_DIR / "supplementary_table_s6_fgf5_lpa_candidate_evidence_matrix.csv", rows, fieldnames)

    lines = [
        "# FGF5/LPA Evidence Matrix",
        "",
        "| Protein | Lead SNP | Primary AF | Primary HF | FinnGen AF | FinnGen HF | Formal primary coloc status | AF mediation | Current interpretation |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for row in rows:
        lines.append(
            f"| {row['protein']} | {row['lead_snp']} | {row['primary_af']} | {row['primary_hf']} | "
            f"{row['finngen_af']} | {row['finngen_hf']} | {row['formal_primary_coloc_status']} | "
            f"{row['af_mediation_status']} | {row['current_interpretation']} |"
        )
    lines.extend(
        [
            "",
            "Allele direction:",
            "",
            "- FGF5: C allele raises FGF5; primary AF and HF point to higher risk after harmonisation.",
            "- LPA: T allele raises apolipoprotein(a); primary AF and HF point to higher risk after harmonisation.",
            "",
            f"CSV table: `{CANDIDATE_MATRIX}`",
        ]
    )
    write_text(CANDIDATE_MATRIX_MD, "\n".join(lines) + "\n")


def update_manifest() -> None:
    rows = read_csv(TABLE_DIR / "supplementary_tables_manifest.csv")
    rows = [row for row in rows if row.get("table_id") != "S12"]
    rows.append(
        {
            "table_id": "S12",
            "description": "FGF5/LPA target-prioritization scorecard",
            "filename": "supplementary_table_s12_target_priority_scorecard.csv",
            "rows": "2",
            "source": str(PRIORITY_DIR / "target_priority_scorecard.csv"),
        }
    )
    rows.sort(key=lambda row: int(str(row.get("table_id", "S0")).replace("S", "") or 0))
    fieldnames = ["table_id", "description", "filename", "rows", "source"]
    write_csv(TABLE_DIR / "supplementary_tables_manifest.csv", rows, fieldnames)

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
            "Note: Results are first-pass single-variant Wald ratio estimates plus formal coloc.abf, target prioritization, feasible sensitivity checks, candidate reverse MR, replication analyses, and exploratory AF-mediated effect estimates.",
        ]
    )
    write_text(TABLE_DIR / "supplementary_tables_manifest.md", "\n".join(lines) + "\n")


def write_log(rows: list[dict[str, Any]]) -> None:
    top = rows[0]
    second = rows[1]
    text = f"""# 靶点优先级排序记录

日期：2026-05-27

## 输入证据

- 主MR：`results/mr/shared_candidate_preliminary.csv`
- 工具变量质量：`results/qc/candidate_instrument_qc.csv`
- FinnGen R12候选复制：`results/replication/finngen_r12_candidate_wald_mr.csv`
- 正式主结局共定位：`results/coloc/formal_coloc_abf_summary.csv`
- AF介导分析：`results/mediation/af_mediated_effects_fgf5_lpa.csv`
- 外部可成药性/临床转化资料：Open Targets靶点优先级/tractability框架、ClinicalTrials.gov Lp(a) outcome trials、FGF5-AF既往MR文献。

## 评分规则

本轮评分是候选优先级排序，不是最终生物学结论。评分维度为：主MR证据0-4分，AF/HF方向一致0-1分，FinnGen复制0-2分，正式共定位0-4分，探索性AF介导支持0-1分，工具变量与等位基因可信度0-1分，临床可成药性0-2分，并根据关键局限预设扣分。

## 排序结论

1. {top['protein']}：{top['priority_class']}，总分{top['overall_score']}。{top['recommended_interpretation']}
2. {second['protein']}：{second['priority_class']}，总分{second['overall_score']}。{second['recommended_interpretation']}

## 稿件解释建议

- FGF5应写作“当前遗传证据优先级最高的候选”，核心理由是AF主MR、FinnGen AF复制、等位基因方向、正式FGF5-AF共定位和AF介导路径支持较一致。限制是HF只达到名义显著且HF正式共定位不足，所以不应写成“已证实AF-HF双结局共定位靶点”。
- LPA应写作“转化可开发性最高但遗传定位需谨慎的候选”。核心理由是HF主MR、AF/HF FinnGen复制、较小但名义显著的AF介导间接效应和Lp(a)药物开发链条强，但本研究正式共定位未支持AF/HF共享单一因果变异，尤其HF显示高PP.H3。
- Introduction和Discussion初稿已完成，后续应随目标期刊和全文英文统一进一步精修。

## 输出文件

- `results/prioritization/target_priority_scorecard.md`
- `results/prioritization/target_priority_scorecard.csv`
- `tables/supplementary_table_s12_target_priority_scorecard.csv`
- `results/candidates/fgf5_lpa_evidence_matrix.md`（已同步正式共定位状态）
"""
    write_text(PROJECT_ROOT / "target_prioritization_log_2026-05-27.md", text)


def main() -> None:
    rows = build_scorecard()
    write_scorecard(rows)
    update_candidate_matrix()
    update_manifest()
    write_log(rows)
    print(f"Wrote target prioritization scorecard for {len(rows)} candidates.")


if __name__ == "__main__":
    main()
