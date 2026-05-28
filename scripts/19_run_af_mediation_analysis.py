"""Run AF-mediated effects for inflammatory proteins on HF.

The mediation path is:

    genetically predicted protein level -> AF -> HF

The AF -> HF step is estimated from local AF and HF GWAS summary statistics
using genome-wide significant AF instruments and distance pruning, because no
LD reference panel is available in this workspace.
"""

from __future__ import annotations

import csv
import gzip
import io
import math
import zipfile
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW_OUTCOMES = PROJECT_ROOT / "data" / "raw" / "outcomes"
MR_ALL = PROJECT_ROOT / "results" / "mr" / "wald_mr_all_outcomes.csv"
SHARED_CANDIDATES = PROJECT_ROOT / "results" / "mr" / "shared_candidate_preliminary.csv"
NOMINAL_OVERLAP = PROJECT_ROOT / "results" / "overlap" / "af_hf_nominal_overlap_proteins.csv"
MED_DIR = PROJECT_ROOT / "results" / "mediation"
TEXT_DIR = PROJECT_ROOT / "results" / "text"
TABLE_DIR = PROJECT_ROOT / "tables"

P_THRESHOLD = 5e-8
DISTANCE_PRUNE_BP = 10_000_000

OUTCOMES = {
    "AF": {
        "name": "Atrial fibrillation",
        "file": RAW_OUTCOMES / "nielsen-thorolfsdottir-willer-NG2018-AFib-gwas-summary-statistics.tbl.gz",
        "snp": "rs_dbSNP147",
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
        "name": "Heart failure",
        "file": RAW_OUTCOMES / "HERMES_Jan2019_HeartFailure_summary_data.txt.zip",
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

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
PALINDROMIC = {frozenset({"A", "T"}), frozenset({"C", "G"})}


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
        if value in {"", None, "NA", "nan"}:
            return None
        number = float(value)
        if number != number or math.isinf(number):
            return None
        return number
    except Exception:
        return None


def norm_allele(value: Any) -> str:
    return str(value or "").strip().upper()


def is_snp_allele(value: Any) -> bool:
    return norm_allele(value) in {"A", "C", "G", "T"}


def is_palindromic(a1: Any, a2: Any) -> bool:
    return frozenset({norm_allele(a1), norm_allele(a2)}) in PALINDROMIC


def complement(value: Any) -> str:
    allele = norm_allele(value)
    if len(allele) != 1:
        return allele
    return COMPLEMENT.get(allele, allele)


def p_from_z(z: float) -> float:
    if z != z:
        return float("nan")
    return math.erfc(abs(z) / math.sqrt(2))


def ci(beta: float, se: float) -> tuple[float, float]:
    return beta - 1.96 * se, beta + 1.96 * se


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
    if number < 0.001:
        return f"{number:.2e}"
    return f"{number:.4f}".rstrip("0").rstrip(".")


def open_text(path: Path) -> Iterable[str]:
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            yield from handle
    elif path.suffix == ".zip":
        with zipfile.ZipFile(path) as archive:
            members = [name for name in archive.namelist() if not name.endswith("/")]
            if not members:
                raise FileNotFoundError(f"No file inside {path}")
            with archive.open(members[0]) as raw:
                text = io.TextIOWrapper(raw, encoding="utf-8", errors="replace")
                yield from text
    else:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            yield from handle


def detect_delimiter(line: str) -> str:
    return "\t" if "\t" in line else ","


def _chain_header(header: str, iterator: Iterable[str]) -> Iterable[str]:
    yield header
    yield from iterator


def dict_reader(path: Path) -> Iterable[dict[str, str]]:
    iterator = iter(open_text(path))
    header = next(iterator)
    delimiter = detect_delimiter(header)
    return csv.DictReader(_chain_header(header, iterator), delimiter=delimiter)


def extract_disease_instruments(outcome_id: str, config: dict[str, Any]) -> tuple[list[dict[str, Any]], int]:
    raw: list[dict[str, Any]] = []
    for row in dict_reader(config["file"]):
        p = fnum(row.get(config["pval"]))
        if p is None or p > P_THRESHOLD:
            continue
        chr_ = str(row.get(config["chr"], "")).replace("chr", "").replace("CHR", "")
        if chr_ not in {str(i) for i in range(1, 23)}:
            continue
        pos = fnum(row.get(config["pos"]))
        beta = fnum(row.get(config["beta"]))
        se = fnum(row.get(config["se"]))
        ea = norm_allele(row.get(config["effect_allele"], ""))
        oa = norm_allele(row.get(config["other_allele"], ""))
        if pos is None or beta is None or se is None or se <= 0:
            continue
        if not is_snp_allele(ea) or not is_snp_allele(oa):
            continue
        snp = row.get(config["snp"], "")
        if not snp or not snp.startswith("rs"):
            continue
        raw.append(
            {
                "outcome_id": outcome_id,
                "outcome_name": config["name"],
                "SNP": snp,
                "chr": chr_,
                "pos": int(pos),
                "effect_allele": ea,
                "other_allele": oa,
                "eaf": row.get(config["eaf"], ""),
                "beta": beta,
                "se": se,
                "pval": p,
                "abs_z": abs(beta / se),
                "f_stat": (beta / se) ** 2,
            }
        )

    selected: list[dict[str, Any]] = []
    by_chr_selected: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in sorted(raw, key=lambda item: (item["pval"], -item["abs_z"])):
        chr_selected = by_chr_selected[row["chr"]]
        if any(abs(row["pos"] - kept["pos"]) <= DISTANCE_PRUNE_BP for kept in chr_selected):
            continue
        selected.append(row)
        chr_selected.append(row)
    return selected, len(raw)


def allele_score(exposure: dict[str, Any], outcome: dict[str, Any]) -> int:
    exp_ea = norm_allele(exposure["effect_allele"])
    exp_oa = norm_allele(exposure["other_allele"])
    out_ea = norm_allele(outcome["effect_allele"])
    out_oa = norm_allele(outcome["other_allele"])
    if exp_ea == out_ea and exp_oa == out_oa:
        return 2
    if exp_ea == out_oa and exp_oa == out_ea:
        return 2
    if exp_ea in {out_ea, out_oa} or exp_oa in {out_ea, out_oa}:
        return 1
    return 0


def extract_outcome_rows_for_instruments(
    config: dict[str, Any], instruments: list[dict[str, Any]]
) -> dict[str, dict[str, Any]]:
    by_snp = {row["SNP"]: row for row in instruments}
    chrpos_to_snp = {f"{row['chr']}:{row['pos']}": row["SNP"] for row in instruments}
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)

    for row in dict_reader(config["file"]):
        snp = row.get(config["snp"], "")
        pos_key = f"{str(row.get(config['chr'], '')).replace('chr', '').replace('CHR', '')}:{row.get(config['pos'], '')}"
        matched_snp = snp if snp in by_snp else chrpos_to_snp.get(pos_key, "")
        if not matched_snp:
            continue
        beta = fnum(row.get(config["beta"]))
        se = fnum(row.get(config["se"]))
        if beta is None or se is None or se <= 0:
            continue
        out = {
            "SNP": matched_snp,
            "chr": str(row.get(config["chr"], "")).replace("chr", "").replace("CHR", ""),
            "pos": row.get(config["pos"], ""),
            "effect_allele": norm_allele(row.get(config["effect_allele"], "")),
            "other_allele": norm_allele(row.get(config["other_allele"], "")),
            "eaf": row.get(config["eaf"], ""),
            "beta": beta,
            "se": se,
            "pval": fnum(row.get(config["pval"])),
        }
        grouped[matched_snp].append(out)

    selected: dict[str, dict[str, Any]] = {}
    for snp, group in grouped.items():
        exposure = by_snp[snp]

        def key(outcome: dict[str, Any]) -> tuple[int, float]:
            score = allele_score(exposure, outcome)
            p = outcome.get("pval")
            p = p if isinstance(p, float) else 1.0
            return (score, -p)

        selected[snp] = max(group, key=key)
    return selected


def harmonise_outcome_to_exposure(
    exposure: dict[str, Any], outcome: dict[str, Any]
) -> tuple[float | None, str]:
    exp_ea = norm_allele(exposure["effect_allele"])
    exp_oa = norm_allele(exposure["other_allele"])
    out_ea = norm_allele(outcome["effect_allele"])
    out_oa = norm_allele(outcome["other_allele"])
    beta = fnum(outcome.get("beta"))
    if beta is None:
        return None, "missing_outcome_beta"
    if exp_ea == out_ea and exp_oa == out_oa:
        return beta, "aligned"
    if exp_ea == out_oa and exp_oa == out_ea:
        return -beta, "flipped"
    if complement(exp_ea) == out_ea and complement(exp_oa) == out_oa:
        return beta, "aligned_complement"
    if complement(exp_ea) == out_oa and complement(exp_oa) == out_ea:
        return -beta, "flipped_complement"
    return None, "allele_mismatch"


def build_af_to_hf_harmonised() -> tuple[list[dict[str, Any]], dict[str, Any]]:
    af_instruments, raw_count = extract_disease_instruments("AF", OUTCOMES["AF"])
    hf_rows = extract_outcome_rows_for_instruments(OUTCOMES["HF"], af_instruments)

    harmonised: list[dict[str, Any]] = []
    skipped_pal = 0
    allele_mismatch = 0
    for inst in af_instruments:
        hf = hf_rows.get(inst["SNP"])
        if not hf:
            continue
        pal = is_palindromic(inst["effect_allele"], inst["other_allele"])
        eaf = fnum(inst.get("eaf"))
        if pal and eaf is not None and 0.42 <= eaf <= 0.58:
            skipped_pal += 1
            continue
        hf_beta_aligned, action = harmonise_outcome_to_exposure(inst, hf)
        if hf_beta_aligned is None:
            allele_mismatch += 1
            continue
        harmonised.append(
            {
                "SNP": inst["SNP"],
                "chr": inst["chr"],
                "pos": inst["pos"],
                "af_effect_allele": inst["effect_allele"],
                "af_other_allele": inst["other_allele"],
                "af_eaf": inst["eaf"],
                "af_beta": inst["beta"],
                "af_se": inst["se"],
                "af_p": inst["pval"],
                "af_f_stat": inst["f_stat"],
                "hf_effect_allele": hf["effect_allele"],
                "hf_other_allele": hf["other_allele"],
                "hf_beta_aligned": hf_beta_aligned,
                "hf_se": hf["se"],
                "hf_p": hf["pval"],
                "harmonise_action": action,
            }
        )
    metadata = {
        "raw_af_p5e_8_variants": raw_count,
        "distance_pruned_af_instruments": len(af_instruments),
        "hf_matched_instruments": len(hf_rows),
        "harmonised_instruments": len(harmonised),
        "skipped_ambiguous_palindromic": skipped_pal,
        "allele_mismatches": allele_mismatch,
    }
    return harmonised, metadata


def ivw_af_to_hf(rows: list[dict[str, Any]]) -> dict[str, Any]:
    denom = 0.0
    numer = 0.0
    for row in rows:
        bx = float(row["af_beta"])
        by = float(row["hf_beta_aligned"])
        sy = float(row["hf_se"])
        w = 1 / (sy * sy)
        denom += bx * bx * w
        numer += bx * by * w
    if denom == 0:
        raise ValueError("No denominator for AF->HF IVW")
    beta = numer / denom
    se = math.sqrt(1 / denom)
    p = p_from_z(beta / se)
    q = sum(((float(row["hf_beta_aligned"]) - beta * float(row["af_beta"])) ** 2) / (float(row["hf_se"]) ** 2) for row in rows)
    lo, hi = ci(beta, se)
    return {
        "exposure": "AF genetic liability",
        "outcome": "HF risk",
        "method": "Fixed-effect IVW",
        "n_instruments": len(rows),
        "beta": beta,
        "se": se,
        "or": math.exp(beta),
        "or_lci95": math.exp(lo),
        "or_uci95": math.exp(hi),
        "p": p,
        "q_heterogeneity": q,
    }


def mr_rows_by_protein() -> dict[str, dict[str, dict[str, str]]]:
    grouped: dict[str, dict[str, dict[str, str]]] = defaultdict(dict)
    for row in read_csv(MR_ALL):
        grouped[row["protein"]][row["outcome_id"]] = row
    return grouped


def group_labels() -> tuple[set[str], set[str]]:
    shared = {row["protein"] for row in read_csv(SHARED_CANDIDATES)}
    nominal = {row["protein"] for row in read_csv(NOMINAL_OVERLAP)}
    return shared, nominal


def mediation_rows(af_to_hf: dict[str, Any]) -> list[dict[str, Any]]:
    grouped = mr_rows_by_protein()
    shared, nominal = group_labels()
    rows: list[dict[str, Any]] = []
    b = float(af_to_hf["beta"])
    se_b = float(af_to_hf["se"])
    for protein, by_outcome in sorted(grouped.items()):
        af = by_outcome.get("AF")
        hf = by_outcome.get("HF")
        if not af or not hf:
            continue
        a = fnum(af.get("beta_mr"))
        se_a = fnum(af.get("se_mr"))
        c_total = fnum(hf.get("beta_mr"))
        se_c = fnum(hf.get("se_mr"))
        if a is None or se_a is None or c_total is None or se_c is None:
            continue
        indirect = a * b
        se_indirect = math.sqrt((b * b * se_a * se_a) + (a * a * se_b * se_b))
        p_indirect = p_from_z(indirect / se_indirect) if se_indirect > 0 else float("nan")
        ind_lo, ind_hi = ci(indirect, se_indirect)
        direct = c_total - indirect
        se_direct = math.sqrt(se_c * se_c + se_indirect * se_indirect)
        direct_lo, direct_hi = ci(direct, se_direct)
        prop = indirect / c_total if abs(c_total) > 1e-12 else float("nan")
        same_direction_total = (indirect > 0 and c_total > 0) or (indirect < 0 and c_total < 0)
        if protein in shared:
            label = "FGF5/LPA shared candidate"
        elif protein in nominal:
            label = "AF-HF nominal overlap"
        else:
            label = "all primary proteins"
        if p_indirect < 0.05 and same_direction_total:
            interp = "AF-mediated pathway supported"
        elif p_indirect < 0.05:
            interp = "indirect effect significant but opposite to total"
        else:
            interp = "no nominal AF-mediated support"
        rows.append(
            {
                "protein": protein,
                "protein_name": af.get("protein_name", ""),
                "panel": af.get("panel", ""),
                "SNP": af.get("SNP", ""),
                "analysis_group": label,
                "beta_protein_to_af": a,
                "se_protein_to_af": se_a,
                "or_protein_to_af": fnum(af.get("or")),
                "p_protein_to_af": fnum(af.get("pval_mr")),
                "fdr_protein_to_af": fnum(af.get("fdr")),
                "beta_af_to_hf": b,
                "se_af_to_hf": se_b,
                "or_af_to_hf": af_to_hf["or"],
                "p_af_to_hf": af_to_hf["p"],
                "beta_total_protein_to_hf": c_total,
                "se_total_protein_to_hf": se_c,
                "or_total_protein_to_hf": fnum(hf.get("or")),
                "p_total_protein_to_hf": fnum(hf.get("pval_mr")),
                "fdr_total_protein_to_hf": fnum(hf.get("fdr")),
                "beta_indirect": indirect,
                "se_indirect": se_indirect,
                "or_indirect": math.exp(indirect),
                "or_indirect_lci95": math.exp(ind_lo),
                "or_indirect_uci95": math.exp(ind_hi),
                "p_indirect": p_indirect,
                "beta_direct": direct,
                "se_direct": se_direct,
                "or_direct": math.exp(direct),
                "or_direct_lci95": math.exp(direct_lo),
                "or_direct_uci95": math.exp(direct_hi),
                "proportion_mediated": prop,
                "same_direction_indirect_total": same_direction_total,
                "interpretation": interp,
            }
        )
    return rows


MEDIATION_FIELDS = [
    "protein",
    "protein_name",
    "panel",
    "SNP",
    "analysis_group",
    "beta_protein_to_af",
    "se_protein_to_af",
    "or_protein_to_af",
    "p_protein_to_af",
    "fdr_protein_to_af",
    "beta_af_to_hf",
    "se_af_to_hf",
    "or_af_to_hf",
    "p_af_to_hf",
    "beta_total_protein_to_hf",
    "se_total_protein_to_hf",
    "or_total_protein_to_hf",
    "p_total_protein_to_hf",
    "fdr_total_protein_to_hf",
    "beta_indirect",
    "se_indirect",
    "or_indirect",
    "or_indirect_lci95",
    "or_indirect_uci95",
    "p_indirect",
    "beta_direct",
    "se_direct",
    "or_direct",
    "or_direct_lci95",
    "or_direct_uci95",
    "proportion_mediated",
    "same_direction_indirect_total",
    "interpretation",
]


AF_HF_FIELDS = [
    "SNP",
    "chr",
    "pos",
    "af_effect_allele",
    "af_other_allele",
    "af_eaf",
    "af_beta",
    "af_se",
    "af_p",
    "af_f_stat",
    "hf_effect_allele",
    "hf_other_allele",
    "hf_beta_aligned",
    "hf_se",
    "hf_p",
    "harmonise_action",
]


def write_mediation_summary(
    af_to_hf: dict[str, Any],
    metadata: dict[str, Any],
    all_rows: list[dict[str, Any]],
    candidate_rows: list[dict[str, Any]],
    nominal_rows: list[dict[str, Any]],
) -> None:
    lines = [
        "# AF-Mediated Protein-HF Analysis Summary",
        "",
        "Path tested: genetically predicted inflammatory protein level -> AF -> HF.",
        "",
        "## AF -> HF MR",
        "",
        f"- AF instruments: {metadata['raw_af_p5e_8_variants']} genome-wide significant variants; {metadata['distance_pruned_af_instruments']} after {DISTANCE_PRUNE_BP:,} bp distance pruning.",
        f"- Harmonised AF instruments with HF GWAS: {metadata['harmonised_instruments']} (HF matched: {metadata['hf_matched_instruments']}; ambiguous palindromic skipped: {metadata['skipped_ambiguous_palindromic']}; allele mismatches: {metadata['allele_mismatches']}).",
        f"- AF -> HF IVW estimate: OR {fmt(af_to_hf['or'])} (95% CI {fmt(af_to_hf['or_lci95'])}-{fmt(af_to_hf['or_uci95'])}), beta {fmt(af_to_hf['beta'])}, SE {fmt(af_to_hf['se'])}, P={fmt_p(af_to_hf['p'])}.",
        f"- Cochran-style Q for the AF -> HF instrument set: {fmt(af_to_hf['q_heterogeneity'])}.",
        "",
        "## Candidate Mediation Results",
        "",
        "| Protein | Total protein -> HF OR | Indirect OR via AF | Indirect P | Proportion mediated | Interpretation |",
        "|---|---:|---:|---:|---:|---|",
    ]
    for row in candidate_rows:
        prop = row["proportion_mediated"]
        prop_text = "NA" if prop != prop else fmt(prop)
        lines.append(
            f"| {row['protein']} | {fmt(row['or_total_protein_to_hf'])} | {fmt(row['or_indirect'])} "
            f"({fmt(row['or_indirect_lci95'])}-{fmt(row['or_indirect_uci95'])}) | "
            f"{fmt_p(row['p_indirect'])} | {prop_text} | {row['interpretation']} |"
        )
    if not candidate_rows:
        lines.append("| NA | NA | NA | NA | NA | no candidate rows |")
    lines.extend(
        [
            "",
            "## Scope",
            "",
            f"- All proteins with both AF and HF primary MR estimates: {len(all_rows)}.",
            f"- AF-HF nominal-overlap proteins: {len(nominal_rows)}.",
            f"- FGF5/LPA shared candidates: {len(candidate_rows)}.",
            "",
            "## Caveats",
            "",
            "- The AF -> HF step uses distance-pruned genome-wide significant AF variants because a formal LD reference panel is not available in this workspace.",
            "- Product-of-coefficients mediation is reported on the log-odds scale and should be treated as exploratory for binary outcomes.",
            "- The indirect effect assumes the protein -> AF and AF -> HF estimates are approximately independent and does not replace multivariable MR with individual-level data.",
            "",
            f"AF-HF instrument table: `{MED_DIR / 'af_to_hf_harmonised_instruments.csv'}`",
            f"All mediation table: `{MED_DIR / 'af_mediated_effects_all_primary.csv'}`",
            f"Candidate mediation table: `{MED_DIR / 'af_mediated_effects_fgf5_lpa.csv'}`",
        ]
    )
    write_text(MED_DIR / "af_mediated_effects_summary.md", "\n".join(lines) + "\n")


def write_results_draft(candidate_rows: list[dict[str, Any]], af_to_hf: dict[str, Any], metadata: dict[str, Any]) -> None:
    by_protein = {row["protein"]: row for row in candidate_rows}
    fgf5 = by_protein.get("FGF5")
    lpa = by_protein.get("LPA")

    def sentence_for(row: dict[str, Any] | None) -> str:
        if row is None:
            return ""
        prop = row["proportion_mediated"]
        prop_text = "not estimable" if prop != prop else f"{prop * 100:.1f}%"
        return (
            f"For {row['protein']}, the AF-mediated indirect effect had OR {fmt(row['or_indirect'])} "
            f"(95% CI {fmt(row['or_indirect_lci95'])}-{fmt(row['or_indirect_uci95'])}; "
            f"P={fmt_p(row['p_indirect'])}), corresponding to an exploratory mediated proportion of {prop_text}."
        )

    lines = [
        "# Results Paragraph Draft: AF-Mediated Effects",
        "",
        "### 3.7 Mediation analysis",
        "",
        (
            f"To explore whether AF may mediate part of the inflammatory protein-HF association, we first estimated the "
            f"genetic effect of AF liability on HF using {metadata['harmonised_instruments']} harmonised AF instruments "
            f"selected at P < 5 x 10^-8 and distance-pruned at {DISTANCE_PRUNE_BP // 1_000_000} Mb. "
            f"The fixed-effect IVW estimate supported a positive AF -> HF association "
            f"(OR {fmt(af_to_hf['or'])}, 95% CI {fmt(af_to_hf['or_lci95'])}-{fmt(af_to_hf['or_uci95'])}; "
            f"P={fmt_p(af_to_hf['p'])})."
        ),
        "",
        (
            sentence_for(fgf5)
            + (" " if fgf5 and lpa else "")
            + sentence_for(lpa)
        ),
        "",
        (
            "These estimates suggest that AF-mediated pathways may contribute to the observed FGF5/LPA-HF associations, "
            "but the analysis is exploratory because it uses summary-level binary outcome estimates, distance-pruned "
            "AF instruments without an LD reference panel, and product-of-coefficients assumptions. The mediation "
            "results were therefore used as supportive evidence in target prioritization rather than as definitive "
            "proof of a mechanistic AF-mediated pathway."
        ),
    ]
    write_text(TEXT_DIR / "results_mediation_draft_2026-05-28.md", "\n".join(lines) + "\n")


def update_manifest(all_rows: list[dict[str, Any]], af_hf_rows: list[dict[str, Any]]) -> None:
    manifest = TABLE_DIR / "supplementary_tables_manifest.csv"
    rows = read_csv(manifest)
    rows = [row for row in rows if row.get("table_id") not in {"S20", "S21"}]
    rows.extend(
        [
            {
                "table_id": "S20",
                "description": "AF-to-HF MR harmonised instruments for mediation analysis",
                "filename": "supplementary_table_s20_af_to_hf_mr_instruments.csv",
                "rows": str(len(af_hf_rows)),
                "source": str(MED_DIR / "af_to_hf_harmonised_instruments.csv"),
            },
            {
                "table_id": "S21",
                "description": "AF-mediated protein-to-HF effects for primary MR proteins",
                "filename": "supplementary_table_s21_af_mediated_effects.csv",
                "rows": str(len(all_rows)),
                "source": str(MED_DIR / "af_mediated_effects_all_primary.csv"),
            },
        ]
    )
    rows.sort(key=lambda row: int(str(row.get("table_id", "S0")).replace("S", "") or 0))
    fieldnames = ["table_id", "description", "filename", "rows", "source"]
    write_csv(manifest, rows, fieldnames)
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


def write_log(af_to_hf: dict[str, Any], metadata: dict[str, Any], candidate_rows: list[dict[str, Any]]) -> None:
    lines = [
        "# AF介导分析记录",
        "",
        "日期：2026-05-28",
        "",
        "## 已完成",
        "",
        f"- 使用AF GWAS中P < 5e-8的变异作为AF遗传工具变量，按{DISTANCE_PRUNE_BP:,} bp距离剪枝。",
        f"- 原始AF显著变异{metadata['raw_af_p5e_8_variants']}个，剪枝后{metadata['distance_pruned_af_instruments']}个，与HF GWAS协调后{metadata['harmonised_instruments']}个。",
        f"- AF -> HF固定效应IVW：OR {fmt(af_to_hf['or'])}，95%CI {fmt(af_to_hf['or_lci95'])}-{fmt(af_to_hf['or_uci95'])}，P={fmt_p(af_to_hf['p'])}。",
        "- 对所有同时具有AF和HF主MR结果的蛋白计算了蛋白 -> AF -> HF间接效应。",
        "",
        "## FGF5/LPA候选结果",
        "",
        "| Protein | Indirect OR | P | Proportion mediated | Interpretation |",
        "|---|---:|---:|---:|---|",
    ]
    for row in candidate_rows:
        prop = row["proportion_mediated"]
        prop_text = "NA" if prop != prop else fmt(prop)
        lines.append(
            f"| {row['protein']} | {fmt(row['or_indirect'])} "
            f"({fmt(row['or_indirect_lci95'])}-{fmt(row['or_indirect_uci95'])}) | "
            f"{fmt_p(row['p_indirect'])} | {prop_text} | {row['interpretation']} |"
        )
    lines.extend(
        [
            "",
            "## 解释限制",
            "",
            "- 这是summary-level二分类结局下的探索性乘积法中介分析。",
            "- AF -> HF工具变量未使用正式LD参考面板clumping，仅使用距离剪枝。",
            "- 中介结果可作为靶点优先级支持证据，但不能单独证明机制链条。",
        ]
    )
    write_text(PROJECT_ROOT / "mediation_analysis_log_2026-05-28.md", "\n".join(lines) + "\n")


def main() -> None:
    af_hf_rows, metadata = build_af_to_hf_harmonised()
    af_to_hf = ivw_af_to_hf(af_hf_rows)
    all_rows = mediation_rows(af_to_hf)
    shared, nominal = group_labels()
    candidate_rows = [row for row in all_rows if row["protein"] in shared]
    nominal_rows = [row for row in all_rows if row["protein"] in nominal]

    write_csv(MED_DIR / "af_to_hf_harmonised_instruments.csv", af_hf_rows, AF_HF_FIELDS)
    write_csv(TABLE_DIR / "supplementary_table_s20_af_to_hf_mr_instruments.csv", af_hf_rows, AF_HF_FIELDS)
    write_csv(MED_DIR / "af_mediated_effects_all_primary.csv", all_rows, MEDIATION_FIELDS)
    write_csv(MED_DIR / "af_mediated_effects_nominal_overlap.csv", nominal_rows, MEDIATION_FIELDS)
    write_csv(MED_DIR / "af_mediated_effects_fgf5_lpa.csv", candidate_rows, MEDIATION_FIELDS)
    write_csv(TABLE_DIR / "supplementary_table_s21_af_mediated_effects.csv", all_rows, MEDIATION_FIELDS)
    write_mediation_summary(af_to_hf, metadata, all_rows, candidate_rows, nominal_rows)
    write_results_draft(candidate_rows, af_to_hf, metadata)
    update_manifest(all_rows, af_hf_rows)
    write_log(af_to_hf, metadata, candidate_rows)

    print(f"AF->HF harmonised instruments: {len(af_hf_rows)}")
    print(f"AF->HF OR: {af_to_hf['or']:.6g}, P={af_to_hf['p']:.4g}")
    print(f"Mediation rows: {len(all_rows)}")
    print(f"Candidate rows: {len(candidate_rows)}")


if __name__ == "__main__":
    main()
