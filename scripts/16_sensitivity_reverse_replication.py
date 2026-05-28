"""Run feasible sensitivity checks and candidate reverse MR analyses."""

from __future__ import annotations

import csv
import gzip
import io
import math
import tarfile
import zipfile
from collections import defaultdict
from pathlib import Path
from statistics import NormalDist
from typing import Any, Iterable


PROJECT_ROOT = Path(__file__).resolve().parents[1]

EXPOSURE = PROJECT_ROOT / "data" / "processed" / "ukbppp_inflammation_cis_pqtl_instruments_lead_no_mhc.csv"
MR_ALL = PROJECT_ROOT / "results" / "mr" / "wald_mr_all_outcomes.csv"
RAW_OUTCOMES = PROJECT_ROOT / "data" / "raw" / "outcomes"
RAW_PGWAS = PROJECT_ROOT / "data" / "raw" / "ukbppp_pgwas"
SENS_DIR = PROJECT_ROOT / "results" / "sensitivity"
REVERSE_DIR = PROJECT_ROOT / "results" / "reverse_mr"
REPL_DIR = PROJECT_ROOT / "results" / "replication"
TABLE_DIR = PROJECT_ROOT / "tables"

P_THRESHOLD = 5e-8
DISTANCE_PRUNE_BP = 10_000_000
CANDIDATES = {"FGF5", "LPA"}

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

PGWAS = {
    "FGF5": {
        "protein_name": "Fibroblast growth factor 5",
        "tar_pattern": "*fgf5*.tar",
    },
    "LPA": {
        "protein_name": "Apolipoprotein(a)",
        "tar_pattern": "*lpa*.tar",
    },
}

TARGET_CIS_EXCLUDE_GRCH37 = {
    "FGF5": {"chr": "4", "center": 81182554, "window": 1_000_000},
    "LPA": {"chr": "6", "center": 161089307, "window": 1_000_000},
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


def is_snp_allele(value: str) -> bool:
    return norm_allele(value) in {"A", "C", "G", "T"}


def is_palindromic(a1: str, a2: str) -> bool:
    return frozenset({norm_allele(a1), norm_allele(a2)}) in PALINDROMIC


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


def dict_reader(path: Path) -> Iterable[dict[str, str]]:
    iterator = iter(open_text(path))
    header = next(iterator)
    delimiter = detect_delimiter(header)
    return csv.DictReader(_chain_header(header, iterator), delimiter=delimiter)


def _chain_header(header: str, iterator: Iterable[str]) -> Iterable[str]:
    yield header
    yield from iterator


def build_sensitivity() -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    exposure = {(row["protein"], row["SNP"]): row for row in read_csv(EXPOSURE)}
    rows = read_csv(MR_ALL)
    steiger_rows: list[dict[str, Any]] = []
    candidate_rows: list[dict[str, Any]] = []
    for row in rows:
        exp = exposure.get((row["protein"], row["SNP"]))
        if not exp:
            continue
        eaf = fnum(exp.get("eaf"))
        bx = fnum(row.get("beta_exposure"))
        by = fnum(row.get("beta_outcome"))
        if eaf is None or bx is None or by is None:
            continue
        exposure_r2 = 2 * eaf * (1 - eaf) * bx * bx
        outcome_r2_proxy = 2 * eaf * (1 - eaf) * by * by
        direction = "protein_to_disease_supported" if exposure_r2 > outcome_r2_proxy else "direction_uncertain"
        out = {
            "protein": row["protein"],
            "protein_name": row["protein_name"],
            "outcome_id": row["outcome_id"],
            "SNP": row["SNP"],
            "effect_allele": exp["effect_allele"],
            "other_allele": exp["other_allele"],
            "palindromic_alleles": is_palindromic(exp["effect_allele"], exp["other_allele"]),
            "instrument_f": exp["f_stat"],
            "exposure_r2_proxy": f"{exposure_r2:.8g}",
            "outcome_r2_proxy": f"{outcome_r2_proxy:.8g}",
            "steiger_proxy_direction": direction,
            "harmonise_action": row["harmonise_action"],
            "mr_or": row["or"],
            "mr_p": row["pval_mr"],
            "mr_fdr": row["fdr"],
        }
        steiger_rows.append(out)
        if row["protein"] in CANDIDATES:
            candidate_rows.append(out)

    fieldnames = [
        "protein",
        "protein_name",
        "outcome_id",
        "SNP",
        "effect_allele",
        "other_allele",
        "palindromic_alleles",
        "instrument_f",
        "exposure_r2_proxy",
        "outcome_r2_proxy",
        "steiger_proxy_direction",
        "harmonise_action",
        "mr_or",
        "mr_p",
        "mr_fdr",
    ]
    write_csv(SENS_DIR / "proxy_steiger_directionality_all_primary_mr.csv", steiger_rows, fieldnames)
    write_csv(SENS_DIR / "candidate_sensitivity_summary.csv", candidate_rows, fieldnames)
    write_csv(TABLE_DIR / "supplementary_table_s13_proxy_steiger_directionality.csv", steiger_rows, fieldnames)

    total = len(steiger_rows)
    supported = sum(1 for row in steiger_rows if row["steiger_proxy_direction"] == "protein_to_disease_supported")
    palindromic = sum(1 for row in steiger_rows if row["palindromic_alleles"])
    weak = sum(1 for row in steiger_rows if (fnum(row["instrument_f"]) or 0.0) <= 10)
    candidate_supported = sum(1 for row in candidate_rows if row["steiger_proxy_direction"] == "protein_to_disease_supported")
    lines = [
        "# Sensitivity Analysis Summary",
        "",
        "This project uses one lead cis-pQTL per protein for the primary analysis. Therefore MR-Egger, weighted median, MR-PRESSO and leave-one-out are not statistically applicable to the main single-variant estimates.",
        "",
        "Feasible sensitivity checks completed in this step:",
        "",
        f"- Strong-instrument check: {weak} of {total} harmonised primary MR rows had F <= 10.",
        f"- Allele ambiguity check: {palindromic} of {total} harmonised primary MR rows used palindromic alleles; FGF5 and LPA candidate instruments are not palindromic.",
        f"- Proxy-Steiger directionality: {supported} of {total} harmonised primary MR rows had exposure R2 proxy greater than outcome R2 proxy.",
        f"- Candidate proxy-Steiger directionality: {candidate_supported} of {len(candidate_rows)} FGF5/LPA candidate rows supported the protein-to-disease direction.",
        "- LD/confounding sensitivity: formal coloc.abf was completed for FGF5/LPA against AF and HF primary outcomes; FGF5-AF was strongly supported, whereas FGF5-HF, LPA-AF and LPA-HF were not strongly supported.",
        "",
        "Important caveat: proxy-Steiger values here are approximate because binary outcome variance is represented on the observed log-odds scale and not transformed to liability-scale R2.",
        "",
        f"All-row table: `{SENS_DIR / 'proxy_steiger_directionality_all_primary_mr.csv'}`",
        f"Candidate table: `{SENS_DIR / 'candidate_sensitivity_summary.csv'}`",
    ]
    write_text(SENS_DIR / "sensitivity_summary.md", "\n".join(lines) + "\n")
    return steiger_rows, candidate_rows


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
        if pos is None or beta is None or se is None or not is_snp_allele(ea) or not is_snp_allele(oa):
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
                "abs_z": abs(beta / se) if se else 0.0,
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


def find_tar(pattern: str) -> Path | None:
    matches = sorted(RAW_PGWAS.glob(pattern))
    return matches[0] if matches else None


def tar_member_by_chr(tar: tarfile.TarFile) -> dict[str, tarfile.TarInfo]:
    out = {}
    for member in tar.getmembers():
        if not member.isfile() or not member.name.endswith(".gz"):
            continue
        name = member.name.lower()
        marker = "discovery_chr"
        if marker not in name:
            continue
        chrom = name.split(marker, 1)[1].split("_", 1)[0].upper()
        out[chrom] = member
    return out


def harmonise_reverse(disease: dict[str, Any], pqtl: dict[str, str]) -> tuple[float | None, str]:
    disease_ea = norm_allele(disease["effect_allele"])
    disease_oa = norm_allele(disease["other_allele"])
    pqtl_ea = norm_allele(pqtl["ALLELE1"])
    pqtl_oa = norm_allele(pqtl["ALLELE0"])
    beta = fnum(pqtl.get("BETA"))
    if beta is None:
        return None, "missing_pqtl_beta"
    if disease_ea == pqtl_ea and disease_oa == pqtl_oa:
        return beta, "aligned"
    if disease_ea == pqtl_oa and disease_oa == pqtl_ea:
        return -beta, "flipped"
    if complement(disease_ea) == pqtl_ea and complement(disease_oa) == pqtl_oa:
        return beta, "aligned_complement"
    if complement(disease_ea) == pqtl_oa and complement(disease_oa) == pqtl_ea:
        return -beta, "flipped_complement"
    return None, "allele_mismatch"


def pqtl_id_pos37(pqtl: dict[str, str]) -> int | None:
    parts = (pqtl.get("ID") or "").split(":")
    if len(parts) >= 2:
        pos = fnum(parts[1])
        if pos is not None:
            return int(pos)
    pos = fnum(pqtl.get("GENPOS"))
    return int(pos) if pos is not None else None


def reverse_matches_for_protein(protein: str, instruments: list[dict[str, Any]]) -> list[dict[str, Any]]:
    tar_path = find_tar(PGWAS[protein]["tar_pattern"])
    if not tar_path:
        return []
    by_chr: dict[str, dict[int, list[dict[str, Any]]]] = defaultdict(lambda: defaultdict(list))
    for inst in instruments:
        by_chr[inst["chr"]][inst["pos"]].append(inst)
    needed_chroms = set(by_chr)
    matches: list[dict[str, Any]] = []
    with tarfile.open(tar_path) as tar:
        members = tar_member_by_chr(tar)
        for chrom in sorted(needed_chroms, key=lambda x: int(x)):
            member = members.get(chrom)
            if member is None:
                continue
            raw = tar.extractfile(member)
            if raw is None:
                continue
            with gzip.open(raw, "rt", encoding="utf-8", errors="replace") as handle:
                header = handle.readline().strip().split()
                for line in handle:
                    parts = line.strip().split()
                    if len(parts) != len(header):
                        continue
                    pqtl = dict(zip(header, parts))
                    pos_int = pqtl_id_pos37(pqtl)
                    if pos_int is None:
                        continue
                    if pos_int not in by_chr[chrom]:
                        continue
                    for disease in by_chr[chrom][pos_int]:
                        pal = is_palindromic(disease["effect_allele"], disease["other_allele"])
                        eaf = fnum(disease.get("eaf"))
                        if pal and eaf is not None and 0.42 <= eaf <= 0.58:
                            continue
                        aligned_beta, action = harmonise_reverse(disease, pqtl)
                        if aligned_beta is None:
                            continue
                        se = fnum(pqtl.get("SE"))
                        if se is None or se <= 0:
                            continue
                        matches.append(
                            {
                                "protein": protein,
                                "protein_name": PGWAS[protein]["protein_name"],
                                "outcome_id": disease["outcome_id"],
                                "outcome_name": disease["outcome_name"],
                                "SNP": disease["SNP"],
                                "chr": disease["chr"],
                                "pos": disease["pos"],
                                "disease_effect_allele": disease["effect_allele"],
                                "disease_other_allele": disease["other_allele"],
                                "disease_beta": disease["beta"],
                                "disease_se": disease["se"],
                                "disease_p": disease["pval"],
                                "protein_beta_aligned": aligned_beta,
                                "protein_se": se,
                                "harmonise_action": action,
                                "pqtl_allele1": pqtl.get("ALLELE1", ""),
                                "pqtl_allele0": pqtl.get("ALLELE0", ""),
                                "pqtl_a1freq": pqtl.get("A1FREQ", ""),
                                "pqtl_log10p": pqtl.get("LOG10P", ""),
                            }
                        )
    return matches


def ivw_reverse(rows: list[dict[str, Any]]) -> dict[str, Any] | None:
    if not rows:
        return None
    denom = 0.0
    numer = 0.0
    for row in rows:
        bx = float(row["disease_beta"])
        by = float(row["protein_beta_aligned"])
        sy = float(row["protein_se"])
        w = 1 / (sy * sy)
        denom += bx * bx * w
        numer += bx * by * w
    if denom == 0:
        return None
    beta = numer / denom
    se = math.sqrt(1 / denom)
    z = beta / se if se else float("nan")
    normal = NormalDist()
    p = 2 * (1 - normal.cdf(abs(z))) if z == z else float("nan")
    q = sum(((float(row["protein_beta_aligned"]) - beta * float(row["disease_beta"])) ** 2) / (float(row["protein_se"]) ** 2) for row in rows)
    method = "Wald ratio" if len(rows) == 1 else "Fixed-effect IVW"
    return {
        "method": method,
        "n_instruments": len(rows),
        "beta_reverse": beta,
        "se_reverse": se,
        "p_reverse": p,
        "q_heterogeneity": q,
    }


def outside_target_cis(row: dict[str, Any]) -> bool:
    target = TARGET_CIS_EXCLUDE_GRCH37.get(row["protein"])
    if not target:
        return True
    if str(row["chr"]) != target["chr"]:
        return True
    return abs(int(row["pos"]) - int(target["center"])) > int(target["window"])


def run_reverse_mr() -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    all_instruments: list[dict[str, Any]] = []
    instrument_summary: list[dict[str, Any]] = []
    for outcome_id, config in OUTCOMES.items():
        instruments, raw_count = extract_disease_instruments(outcome_id, config)
        all_instruments.extend(instruments)
        instrument_summary.append(
            {
                "outcome_id": outcome_id,
                "outcome_name": config["name"],
                "raw_p5e_8_variants": raw_count,
                "distance_pruned_instruments": len(instruments),
                "distance_prune_window_bp": DISTANCE_PRUNE_BP,
            }
        )

    match_rows: list[dict[str, Any]] = []
    for protein in sorted(CANDIDATES):
        match_rows.extend(reverse_matches_for_protein(protein, all_instruments))

    result_rows: list[dict[str, Any]] = []
    grouped: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in match_rows:
        grouped[(row["protein"], row["outcome_id"])].append(row)
    for (protein, outcome_id), rows in sorted(grouped.items()):
        analysis_sets = [
            ("all_distance_pruned", rows),
            ("target_cis_excluded", [row for row in rows if outside_target_cis(row)]),
        ]
        cis_excluded_count = len(rows) - len(analysis_sets[1][1])
        for analysis_set, analysis_rows in analysis_sets:
            res = ivw_reverse(analysis_rows)
            if not res:
                result_rows.append(
                    {
                        "analysis_set": analysis_set,
                        "exposure": f"{outcome_id} genetic liability",
                        "outcome_protein": protein,
                        "protein_name": PGWAS[protein]["protein_name"],
                        "method": "not_applicable",
                        "n_instruments": 0,
                        "n_target_cis_excluded": cis_excluded_count if analysis_set == "target_cis_excluded" else 0,
                        "beta_reverse": "",
                        "se_reverse": "",
                        "p_reverse": "",
                        "q_heterogeneity": "",
                        "interpretation": "no_harmonised_instruments",
                    }
                )
                continue
            result_rows.append(
                {
                    "analysis_set": analysis_set,
                    "exposure": f"{outcome_id} genetic liability",
                    "outcome_protein": protein,
                    "protein_name": PGWAS[protein]["protein_name"],
                    "method": res["method"],
                    "n_instruments": res["n_instruments"],
                    "n_target_cis_excluded": cis_excluded_count if analysis_set == "target_cis_excluded" else 0,
                    "beta_reverse": f"{res['beta_reverse']:.8g}",
                    "se_reverse": f"{res['se_reverse']:.8g}",
                    "p_reverse": f"{res['p_reverse']:.8g}",
                    "q_heterogeneity": f"{res['q_heterogeneity']:.8g}",
                    "interpretation": "nominal_reverse_signal" if res["p_reverse"] < 0.05 else "no_reverse_signal",
                }
            )

    instrument_fields = ["outcome_id", "outcome_name", "raw_p5e_8_variants", "distance_pruned_instruments", "distance_prune_window_bp"]
    match_fields = [
        "protein",
        "protein_name",
        "outcome_id",
        "outcome_name",
        "SNP",
        "chr",
        "pos",
        "disease_effect_allele",
        "disease_other_allele",
        "disease_beta",
        "disease_se",
        "disease_p",
        "protein_beta_aligned",
        "protein_se",
        "harmonise_action",
        "pqtl_allele1",
        "pqtl_allele0",
        "pqtl_a1freq",
        "pqtl_log10p",
    ]
    result_fields = [
        "analysis_set",
        "exposure",
        "outcome_protein",
        "protein_name",
        "method",
        "n_instruments",
        "n_target_cis_excluded",
        "beta_reverse",
        "se_reverse",
        "p_reverse",
        "q_heterogeneity",
        "interpretation",
    ]
    write_csv(REVERSE_DIR / "disease_instrument_summary.csv", instrument_summary, instrument_fields)
    write_csv(REVERSE_DIR / "candidate_reverse_mr_harmonised_instruments.csv", match_rows, match_fields)
    write_csv(REVERSE_DIR / "candidate_reverse_mr_results.csv", result_rows, result_fields)
    write_csv(TABLE_DIR / "supplementary_table_s14_candidate_reverse_mr_results.csv", result_rows, result_fields)

    lines = [
        "# Candidate Reverse MR Summary",
        "",
        f"Disease instruments were selected at P < {P_THRESHOLD:g} from the local AF and HF GWAS files and distance-pruned at {DISTANCE_PRUNE_BP:,} bp because an LD reference panel is not available in this workspace.",
        "",
        "Two result sets are reported: all distance-pruned disease instruments, and a more conservative target-cis-excluded set that removes disease instruments within +/-1 Mb of the target protein cis locus.",
        "",
        "| Analysis set | Exposure | Outcome protein | Method | Instruments | Cis excluded | Beta | SE | P | Interpretation |",
        "|---|---|---|---|---:|---:|---:|---:|---:|---|",
    ]
    for row in result_rows:
        lines.append(
            f"| {row['analysis_set']} | {row['exposure']} | {row['outcome_protein']} | {row['method']} | "
            f"{row['n_instruments']} | {row['n_target_cis_excluded']} | {fmt(row['beta_reverse'])} | "
            f"{fmt(row['se_reverse'])} | {fmt_p(row['p_reverse'])} | {row['interpretation']} |"
        )
    if not result_rows:
        lines.append("| NA | NA | NA | 0 | NA | NA | NA | no harmonised reverse-MR instruments |")
    lines.extend(
        [
            "",
            "Caveat: this is a candidate-level reverse MR using distance-pruned disease instruments and the available UKB-PPP pGWAS files for FGF5 and LPA. It should be treated as a screening check rather than a definitive bidirectional MR.",
            "",
            f"Instrument summary: `{REVERSE_DIR / 'disease_instrument_summary.csv'}`",
            f"Harmonised instruments: `{REVERSE_DIR / 'candidate_reverse_mr_harmonised_instruments.csv'}`",
            f"Reverse MR results: `{REVERSE_DIR / 'candidate_reverse_mr_results.csv'}`",
        ]
    )
    write_text(REVERSE_DIR / "candidate_reverse_mr_summary.md", "\n".join(lines) + "\n")
    return instrument_summary, match_rows, result_rows


def write_replication_status() -> None:
    rows = [
        {
            "replication_layer": "FinnGen R12 candidate exact-variant replication",
            "status": "completed",
            "scope": "FGF5 and LPA candidate variants for I9_AF and I9_HEARTFAIL",
            "output": "results/replication/finngen_r12_candidate_wald_mr.csv",
            "note": "FGF5 replicated for AF; LPA replicated for AF and strict HF in the same risk-increasing direction.",
        },
        {
            "replication_layer": "FinnGen R12 full-panel replication",
            "status": "not_completed",
            "scope": "All 529 lead cis-pQTL instruments against FinnGen AF/HF",
            "output": "",
            "note": "Requires network/API access for all variants or downloaded FinnGen summary statistics.",
        },
        {
            "replication_layer": "UKB outcome replication",
            "status": "not_completed",
            "scope": "AF/HF replication in an independent UKB or UKB-derived outcome GWAS",
            "output": "",
            "note": "No independent UKB AF/HF summary-statistics file is available locally; UKB overlap with UKB-PPP exposure also needs careful handling.",
        },
    ]
    fields = ["replication_layer", "status", "scope", "output", "note"]
    write_csv(REPL_DIR / "finngen_ukb_replication_status.csv", rows, fields)
    write_csv(TABLE_DIR / "supplementary_table_s15_replication_status.csv", rows, fields)
    lines = [
        "# FinnGen/UKB Replication Status",
        "",
        "| Layer | Status | Scope | Note |",
        "|---|---|---|---|",
    ]
    for row in rows:
        lines.append(f"| {row['replication_layer']} | {row['status']} | {row['scope']} | {row['note']} |")
    lines.extend(
        [
            "",
            "Interpretation: candidate FinnGen replication has been completed. Full-panel FinnGen replication and UKB outcome replication should remain pending until the required external outcome data or API access is available.",
        ]
    )
    write_text(REPL_DIR / "finngen_ukb_replication_status.md", "\n".join(lines) + "\n")


def update_manifest() -> None:
    manifest = TABLE_DIR / "supplementary_tables_manifest.csv"
    rows = read_csv(manifest)
    rows = [row for row in rows if row.get("table_id") not in {"S13", "S14", "S15"}]
    rows.extend(
        [
            {
                "table_id": "S13",
                "description": "Proxy-Steiger directionality and sensitivity checks for primary MR",
                "filename": "supplementary_table_s13_proxy_steiger_directionality.csv",
                "rows": str(len(read_csv(SENS_DIR / "proxy_steiger_directionality_all_primary_mr.csv"))),
                "source": str(SENS_DIR / "proxy_steiger_directionality_all_primary_mr.csv"),
            },
            {
                "table_id": "S14",
                "description": "Candidate reverse MR results for AF/HF liability to FGF5/LPA protein levels",
                "filename": "supplementary_table_s14_candidate_reverse_mr_results.csv",
                "rows": str(len(read_csv(REVERSE_DIR / "candidate_reverse_mr_results.csv"))),
                "source": str(REVERSE_DIR / "candidate_reverse_mr_results.csv"),
            },
            {
                "table_id": "S15",
                "description": "FinnGen/UKB replication status",
                "filename": "supplementary_table_s15_replication_status.csv",
                "rows": "3",
                "source": str(REPL_DIR / "finngen_ukb_replication_status.csv"),
            },
        ]
    )
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
            "Note: Results are first-pass single-variant Wald ratio estimates plus formal coloc.abf, target prioritization, feasible sensitivity checks, candidate reverse MR, and replication status tracking.",
        ]
    )
    write_text(TABLE_DIR / "supplementary_tables_manifest.md", "\n".join(lines) + "\n")


def write_log(
    steiger_rows: list[dict[str, Any]],
    candidate_rows: list[dict[str, Any]],
    instrument_summary: list[dict[str, Any]],
    reverse_results: list[dict[str, Any]],
) -> None:
    supported = sum(1 for row in steiger_rows if row["steiger_proxy_direction"] == "protein_to_disease_supported")
    candidate_supported = sum(1 for row in candidate_rows if row["steiger_proxy_direction"] == "protein_to_disease_supported")
    reverse_nominal_all = sum(
        1
        for row in reverse_results
        if row.get("analysis_set") == "all_distance_pruned" and row["interpretation"] == "nominal_reverse_signal"
    )
    reverse_total_all = sum(1 for row in reverse_results if row.get("analysis_set") == "all_distance_pruned")
    reverse_nominal_cis_excluded = sum(
        1
        for row in reverse_results
        if row.get("analysis_set") == "target_cis_excluded" and row["interpretation"] == "nominal_reverse_signal"
    )
    reverse_total_cis_excluded = sum(1 for row in reverse_results if row.get("analysis_set") == "target_cis_excluded")
    lines = [
        "# 第3周敏感性分析、反向MR和复制状态记录",
        "",
        "日期：2026-05-27",
        "",
        "## 顺序判断",
        "",
        "第3周分析应先于正式Methods定稿完成；Methods v0.1可以在本轮分析后撰写，并把仍缺数据的复制分析写成待完成或预设流程。",
        "",
        "## 已完成",
        "",
        f"- Proxy-Steiger方向性：{supported}/{len(steiger_rows)}条主MR记录支持蛋白到疾病方向；FGF5/LPA候选为{candidate_supported}/{len(candidate_rows)}条。",
        "- 单SNP MR限制说明：主分析每个蛋白使用一个lead cis-pQTL，因此MR-Egger、weighted median、MR-PRESSO和leave-one-out不适用于主估计。",
        "- 正式共定位已作为LD混杂敏感性证据纳入：FGF5-AF支持，其他候选-结局组合不支持强共定位。",
        "- 候选反向MR已完成：AF/HF遗传工具变量来自本地主结局GWAS，蛋白结局来自FGF5/LPA UKB-PPP全基因组pGWAS。",
        f"- 候选反向MR中名义显著反向信号数量：全距离剪枝集{reverse_nominal_all}/{reverse_total_all}；target-cis剔除集{reverse_nominal_cis_excluded}/{reverse_total_cis_excluded}。",
        "- FinnGen候选复制此前已完成；全量FinnGen和UKB复制仍需要额外数据或网络/API访问。",
        "",
        "## 疾病工具变量概况",
        "",
        "| Outcome | Raw P<5e-8 variants | Distance-pruned instruments |",
        "|---|---:|---:|",
    ]
    for row in instrument_summary:
        lines.append(f"| {row['outcome_id']} | {row['raw_p5e_8_variants']} | {row['distance_pruned_instruments']} |")
    lines.extend(
        [
            "",
            "## 输出文件",
            "",
            "- `results/sensitivity/sensitivity_summary.md`",
            "- `results/sensitivity/proxy_steiger_directionality_all_primary_mr.csv`",
            "- `results/reverse_mr/candidate_reverse_mr_summary.md`",
            "- `results/reverse_mr/candidate_reverse_mr_results.csv`",
            "- `results/replication/finngen_ukb_replication_status.md`",
            "- `tables/supplementary_table_s13_proxy_steiger_directionality.csv`",
            "- `tables/supplementary_table_s14_candidate_reverse_mr_results.csv`",
            "- `tables/supplementary_table_s15_replication_status.csv`",
        ],
    )
    write_text(PROJECT_ROOT / "week3_sensitivity_reverse_replication_log_2026-05-27.md", "\n".join(lines) + "\n")


def main() -> None:
    steiger_rows, candidate_rows = build_sensitivity()
    instrument_summary, _matches, reverse_results = run_reverse_mr()
    write_replication_status()
    update_manifest()
    write_log(steiger_rows, candidate_rows, instrument_summary, reverse_results)
    print(f"Sensitivity rows: {len(steiger_rows)}")
    print(f"Candidate sensitivity rows: {len(candidate_rows)}")
    print(f"Reverse MR results: {len(reverse_results)}")


if __name__ == "__main__":
    main()
