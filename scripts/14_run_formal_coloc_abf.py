"""Run formal coloc.abf-style colocalization for FGF5/LPA and AF/HF."""

from __future__ import annotations

import csv
import gzip
import io
import math
import tarfile
import zipfile
from pathlib import Path
from typing import Any, Iterable


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW_PGWAS = PROJECT_ROOT / "data" / "raw" / "ukbppp_pgwas"
RAW_META = PROJECT_ROOT / "data" / "raw" / "ukbppp_metadata"
RAW_OUTCOMES = PROJECT_ROOT / "data" / "raw" / "outcomes"
OUT_INPUT = PROJECT_ROOT / "data" / "coloc_inputs"
OUT_RESULTS = PROJECT_ROOT / "results" / "coloc"

P1 = 1e-4
P2 = 1e-4
P12 = 1e-5
W_PQTL = 0.15**2
W_CC = 0.2**2

CANDIDATES = {
    "FGF5": {
        "protein_name": "Fibroblast growth factor 5",
        "chr": "4",
        "lead_rsid": "rs12509595",
        "lead_pos38": 80261400,
        "tar_pattern": "*fgf5*.tar",
        "chrom_member_key": "discovery_chr4_",
        "map_file": "olink_rsid_map_mac5_info03_b0_7_chr4_patched_v2.tsv.gz",
    },
    "LPA": {
        "protein_name": "Apolipoprotein(a)",
        "chr": "6",
        "lead_rsid": "rs56393506",
        "lead_pos38": 160668275,
        "tar_pattern": "*lpa*.tar",
        "chrom_member_key": "discovery_chr6_",
        "map_file": "olink_rsid_map_mac5_info03_b0_7_chr6_patched_v2.tsv.gz",
    },
}

OUTCOMES = {
    "AF": {
        "name": "Atrial fibrillation",
        "file": "nielsen-thorolfsdottir-willer-NG2018-AFib-gwas-summary-statistics.tbl.gz",
        "effect_allele": "A2",
        "other_allele": "A1",
        "eaf": "Freq_A2",
        "beta": "Effect_A2",
        "se": "StdErr",
        "pval": "Pvalue",
        "snp": "rs_dbSNP147",
        "chr": "CHR",
        "pos": "POS_GRCh37",
        "n": "",
        "n_cases": 60620,
        "n_controls": 970216,
    },
    "HF": {
        "name": "Heart failure",
        "file": "HERMES_Jan2019_HeartFailure_summary_data.txt.zip",
        "effect_allele": "A1",
        "other_allele": "A2",
        "eaf": "freq",
        "beta": "b",
        "se": "se",
        "pval": "p",
        "snp": "SNP",
        "chr": "CHR",
        "pos": "BP",
        "n": "N",
        "n_cases": 47309,
        "n_controls": 930014,
    },
}

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def open_text(path: Path) -> Iterable[str]:
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            yield from handle
    elif path.suffix == ".zip":
        with zipfile.ZipFile(path) as archive:
            members = [name for name in archive.namelist() if not name.endswith("/")]
            with archive.open(members[0]) as raw:
                text = io.TextIOWrapper(raw, encoding="utf-8", errors="replace")
                yield from text
    else:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            yield from handle


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


def fmt(value: Any, digits: int = 4) -> str:
    number = fnum(value)
    if number is None:
        return "NA"
    if number == 0:
        return "0"
    if abs(number) < 0.001 or abs(number) >= 10000:
        return f"{number:.{digits}e}"
    return f"{number:.{digits}f}".rstrip("0").rstrip(".")


def norm_allele(value: str) -> str:
    return (value or "").upper()


def complement(value: str) -> str:
    value = norm_allele(value)
    if len(value) == 1:
        return COMPLEMENT.get(value, value)
    return value


def find_tar(pattern: str) -> Path:
    hits = sorted(RAW_PGWAS.glob(pattern))
    if not hits:
        raise FileNotFoundError(f"No pGWAS tar matching {pattern} in {RAW_PGWAS}")
    return hits[0]


def read_map_region(map_file: str, start38: int, end38: int) -> dict[str, dict[str, str]]:
    mapping: dict[str, dict[str, str]] = {}
    path = RAW_META / map_file
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            pos38 = int(row["POS38"])
            if start38 <= pos38 <= end38:
                mapping[row["ID"]] = row
    return mapping


def read_pqtl_region(protein: str, meta: dict[str, Any], window: int = 1_000_000) -> list[dict[str, Any]]:
    start38 = int(meta["lead_pos38"]) - window
    end38 = int(meta["lead_pos38"]) + window
    id_map = read_map_region(meta["map_file"], start38, end38)
    tar_path = find_tar(meta["tar_pattern"])
    rows: list[dict[str, Any]] = []

    with tarfile.open(tar_path, "r") as tar:
        member = next(m for m in tar.getmembers() if meta["chrom_member_key"] in m.name)
        raw = tar.extractfile(member)
        if raw is None:
            raise FileNotFoundError(member.name)
        with gzip.GzipFile(fileobj=raw) as gz:
            header = gz.readline().decode("utf-8").strip().split()
            for line in gz:
                parts = line.decode("utf-8", errors="replace").strip().split()
                if not parts:
                    continue
                row = dict(zip(header, parts))
                pos38 = int(row["GENPOS"])
                if not (start38 <= pos38 <= end38):
                    continue
                mapped = id_map.get(row["ID"])
                if not mapped:
                    continue
                rows.append(
                    {
                        "protein": protein,
                        "protein_name": meta["protein_name"],
                        "variant_id": row["ID"],
                        "rsid": mapped["rsid"],
                        "chr": row["CHROM"],
                        "pos38": row["GENPOS"],
                        "pos19": mapped["POS19"],
                        "effect_allele_pqtl": row["ALLELE1"],
                        "other_allele_pqtl": row["ALLELE0"],
                        "eaf_pqtl": row["A1FREQ"],
                        "beta_pqtl": row["BETA"],
                        "se_pqtl": row["SE"],
                        "n_pqtl": row["N"],
                        "info_pqtl": row["INFO"],
                        "log10p_pqtl": row["LOG10P"],
                    }
                )
    return rows


def outcome_rows_for_rsids(outcome_id: str, rsids: set[str]) -> list[dict[str, str]]:
    meta = OUTCOMES[outcome_id]
    path = RAW_OUTCOMES / meta["file"]
    iterator = iter(open_text(path))
    header = ""
    for line in iterator:
        if line.strip() and not line.startswith("#"):
            header = line
            break
    delimiter = "\t" if "\t" in header else ","
    fieldnames = next(csv.reader([header], delimiter=delimiter))
    reader = csv.DictReader(iterator, fieldnames=fieldnames, delimiter=delimiter)
    rows: list[dict[str, str]] = []
    for row in reader:
        raw_snp = row.get(meta["snp"], "")
        candidates = [value.strip() for value in raw_snp.replace(",", ";").split(";") if value.strip()]
        hits = [snp for snp in candidates if snp in rsids]
        if not hits and raw_snp in rsids:
            hits = [raw_snp]
        for snp in hits:
            rows.append(
                {
                    "outcome_id": outcome_id,
                    "outcome_name": meta["name"],
                    "rsid": snp,
                    "chr": row.get(meta["chr"], ""),
                    "pos19": row.get(meta["pos"], ""),
                    "effect_allele_outcome": row.get(meta["effect_allele"], ""),
                    "other_allele_outcome": row.get(meta["other_allele"], ""),
                    "eaf_outcome": row.get(meta["eaf"], ""),
                    "beta_outcome_raw": row.get(meta["beta"], ""),
                    "se_outcome": row.get(meta["se"], ""),
                    "p_outcome": row.get(meta["pval"], ""),
                    "n_outcome": row.get(meta["n"], "") if meta["n"] else str(meta["n_cases"] + meta["n_controls"]),
                    "n_cases": str(meta["n_cases"]),
                    "n_controls": str(meta["n_controls"]),
                }
            )
    return rows


def harmonise(pqtl: dict[str, Any], out: dict[str, str]) -> tuple[float | None, float | None, str]:
    ea_p = norm_allele(str(pqtl["effect_allele_pqtl"]))
    oa_p = norm_allele(str(pqtl["other_allele_pqtl"]))
    ea_o = norm_allele(out["effect_allele_outcome"])
    oa_o = norm_allele(out["other_allele_outcome"])
    beta = fnum(out["beta_outcome_raw"])
    eaf = fnum(out["eaf_outcome"])
    if beta is None:
        return None, None, "missing_beta"
    if ea_p == ea_o and oa_p == oa_o:
        return beta, eaf, "aligned"
    if ea_p == oa_o and oa_p == ea_o:
        return -beta, (1 - eaf) if eaf is not None else None, "flipped"
    if complement(ea_p) == ea_o and complement(oa_p) == oa_o:
        return beta, eaf, "aligned_complement"
    if complement(ea_p) == oa_o and complement(oa_p) == ea_o:
        return -beta, (1 - eaf) if eaf is not None else None, "flipped_complement"
    return None, None, "allele_mismatch"


def choose_outcome_row(rows: list[dict[str, str]], pqtl: dict[str, Any]) -> dict[str, str] | None:
    scored = []
    for row in rows:
        beta, _eaf, action = harmonise(pqtl, row)
        compatible = 1 if beta is not None else 0
        pos_match = 1 if row.get("pos19") == str(pqtl.get("pos19")) else 0
        p = fnum(row.get("p_outcome")) or 1.0
        scored.append((compatible, pos_match, -p, row))
    if not scored:
        return None
    return max(scored, key=lambda item: item[:3])[3]


def make_merged_input(protein: str, outcome_id: str, pqtl_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rsids = {str(row["rsid"]) for row in pqtl_rows if row.get("rsid")}
    outcome_by_rsid: dict[str, list[dict[str, str]]] = {}
    for row in outcome_rows_for_rsids(outcome_id, rsids):
        outcome_by_rsid.setdefault(row["rsid"], []).append(row)

    merged: list[dict[str, Any]] = []
    for pqtl in pqtl_rows:
        rsid = str(pqtl["rsid"])
        outcome = choose_outcome_row(outcome_by_rsid.get(rsid, []), pqtl)
        if outcome is None:
            continue
        beta_aligned, eaf_aligned, action = harmonise(pqtl, outcome)
        if beta_aligned is None:
            continue
        row = {
            **pqtl,
            **outcome,
            "beta_outcome": f"{beta_aligned:.10g}",
            "eaf_outcome_aligned": "" if eaf_aligned is None else f"{eaf_aligned:.10g}",
            "harmonise_action": action,
        }
        merged.append(row)

    # Drop duplicate rsIDs after keeping the strongest pQTL signal.
    best: dict[str, dict[str, Any]] = {}
    for row in merged:
        key = row["rsid"]
        score = fnum(row.get("log10p_pqtl")) or 0.0
        if key not in best or score > (fnum(best[key].get("log10p_pqtl")) or 0.0):
            best[key] = row
    return list(best.values())


def log_abf(beta: float, se: float, prior_var: float) -> float:
    var = se * se
    if var <= 0:
        return float("-inf")
    z = beta / se
    r = prior_var / (var + prior_var)
    return 0.5 * (math.log(1 - r) + r * z * z)


def logsumexp(values: list[float]) -> float:
    finite = [value for value in values if value != float("-inf")]
    if not finite:
        return float("-inf")
    max_value = max(finite)
    return max_value + math.log(sum(math.exp(value - max_value) for value in finite))


def logdiffexp(a: float, b: float) -> float:
    if b >= a:
        return float("-inf")
    return a + math.log1p(-math.exp(b - a))


def coloc_abf(rows: list[dict[str, Any]]) -> dict[str, Any]:
    usable = []
    for row in rows:
        bp = fnum(row.get("beta_pqtl"))
        sp = fnum(row.get("se_pqtl"))
        bo = fnum(row.get("beta_outcome"))
        so = fnum(row.get("se_outcome"))
        if None in {bp, sp, bo, so}:
            continue
        labf1 = log_abf(bp, sp, W_PQTL)
        labf2 = log_abf(bo, so, W_CC)
        row["lABF_pqtl"] = f"{labf1:.10g}"
        row["lABF_outcome"] = f"{labf2:.10g}"
        row["lABF_sum"] = f"{labf1 + labf2:.10g}"
        usable.append((row, labf1, labf2))

    if not usable:
        return {"n_snps": 0, "error": "no_usable_snps"}

    l1 = [x[1] for x in usable]
    l2 = [x[2] for x in usable]
    l12 = [x[1] + x[2] for x in usable]
    log_s1 = logsumexp(l1)
    log_s2 = logsumexp(l2)
    log_s12 = logsumexp(l12)
    log_h0 = 0.0
    log_h1 = math.log(P1) + log_s1
    log_h2 = math.log(P2) + log_s2
    log_h3 = math.log(P1) + math.log(P2) + logdiffexp(log_s1 + log_s2, log_s12)
    log_h4 = math.log(P12) + log_s12
    logs = [log_h0, log_h1, log_h2, log_h3, log_h4]
    denom = logsumexp(logs)
    pp = [math.exp(value - denom) for value in logs]

    best_row, _best1, _best2 = max(usable, key=lambda item: item[1] + item[2])
    best_snp_pp_h4 = math.exp((fnum(best_row["lABF_sum"]) or 0.0) - log_s12)
    return {
        "n_snps": len(usable),
        "PP.H0": pp[0],
        "PP.H1": pp[1],
        "PP.H2": pp[2],
        "PP.H3": pp[3],
        "PP.H4": pp[4],
        "lead_shared_rsid": best_row["rsid"],
        "lead_shared_variant_id": best_row["variant_id"],
        "lead_shared_snp_pp_h4": best_snp_pp_h4,
        "lead_pqtl_log10p": best_row["log10p_pqtl"],
        "lead_outcome_p": best_row["p_outcome"],
    }


def main() -> None:
    summary_rows: list[dict[str, Any]] = []
    for protein, meta in CANDIDATES.items():
        pqtl_rows = read_pqtl_region(protein, meta)
        write_csv(
            OUT_INPUT / f"{protein}_ukbppp_pqtl_region_1mb.csv",
            pqtl_rows,
            list(pqtl_rows[0].keys()) if pqtl_rows else [],
        )
        for outcome_id in OUTCOMES:
            merged = make_merged_input(protein, outcome_id, pqtl_rows)
            coloc = coloc_abf(merged)
            input_path = OUT_INPUT / f"{protein}_{outcome_id}_coloc_input_harmonised.csv"
            fieldnames = list(merged[0].keys()) if merged else []
            if merged and "lABF_pqtl" not in fieldnames:
                # coloc_abf mutates usable rows, but keep field order stable.
                fieldnames += ["lABF_pqtl", "lABF_outcome", "lABF_sum"]
            write_csv(input_path, merged, fieldnames)
            summary_rows.append(
                {
                    "protein": protein,
                    "protein_name": meta["protein_name"],
                    "outcome_id": outcome_id,
                    "outcome_name": OUTCOMES[outcome_id]["name"],
                    "region": f"chr{meta['chr']}:{meta['lead_pos38'] - 1_000_000}-{meta['lead_pos38'] + 1_000_000} (GRCh38)",
                    "lead_rsid": meta["lead_rsid"],
                    "pqtl_region_rows": len(pqtl_rows),
                    "harmonised_overlap_snps": coloc.get("n_snps", 0),
                    "PP.H0": fmt(coloc.get("PP.H0")),
                    "PP.H1": fmt(coloc.get("PP.H1")),
                    "PP.H2": fmt(coloc.get("PP.H2")),
                    "PP.H3": fmt(coloc.get("PP.H3")),
                    "PP.H4": fmt(coloc.get("PP.H4")),
                    "lead_shared_rsid": coloc.get("lead_shared_rsid", ""),
                    "lead_shared_variant_id": coloc.get("lead_shared_variant_id", ""),
                    "lead_shared_snp_pp_h4": fmt(coloc.get("lead_shared_snp_pp_h4")),
                    "lead_pqtl_log10p": coloc.get("lead_pqtl_log10p", ""),
                    "lead_outcome_p": coloc.get("lead_outcome_p", ""),
                    "input_file": str(input_path),
                }
            )

    fields = [
        "protein",
        "protein_name",
        "outcome_id",
        "outcome_name",
        "region",
        "lead_rsid",
        "pqtl_region_rows",
        "harmonised_overlap_snps",
        "PP.H0",
        "PP.H1",
        "PP.H2",
        "PP.H3",
        "PP.H4",
        "lead_shared_rsid",
        "lead_shared_variant_id",
        "lead_shared_snp_pp_h4",
        "lead_pqtl_log10p",
        "lead_outcome_p",
        "input_file",
    ]
    write_csv(OUT_RESULTS / "formal_coloc_abf_summary.csv", summary_rows, fields)

    lines = [
        "# Formal Coloc ABF Summary",
        "",
        "Method: coloc.abf-style single causal variant colocalization using dense UKB-PPP European discovery pQTL regional summary statistics and local AF/HF GWAS summary statistics. Priors: p1=1e-4, p2=1e-4, p12=1e-5. Prior variance: pQTL=0.15^2, case-control outcome=0.2^2.",
        "",
        "| Protein | Outcome | SNPs | PP.H3 | PP.H4 | Lead shared SNP | SNP.PP.H4 | Interpretation |",
        "|---|---|---:|---:|---:|---|---:|---|",
    ]
    for row in summary_rows:
        pp4 = fnum(row["PP.H4"]) or 0.0
        pp3 = fnum(row["PP.H3"]) or 0.0
        if pp4 >= 0.8:
            interpretation = "strong colocalization"
        elif pp4 >= 0.5:
            interpretation = "moderate colocalization"
        elif pp3 > pp4:
            interpretation = "distinct-signal/LD concern"
        else:
            interpretation = "limited colocalization support"
        lines.append(
            f"| {row['protein']} | {row['outcome_id']} | {row['harmonised_overlap_snps']} | "
            f"{row['PP.H3']} | {row['PP.H4']} | {row['lead_shared_rsid']} | "
            f"{row['lead_shared_snp_pp_h4']} | {interpretation} |"
        )
    lines += [
        "",
        f"CSV summary: `{OUT_RESULTS / 'formal_coloc_abf_summary.csv'}`",
    ]
    write_text(OUT_RESULTS / "formal_coloc_abf_summary.md", "\n".join(lines))
    print(f"Wrote {OUT_RESULTS / 'formal_coloc_abf_summary.csv'}")
    print(f"Wrote {OUT_RESULTS / 'formal_coloc_abf_summary.md'}")


if __name__ == "__main__":
    main()
