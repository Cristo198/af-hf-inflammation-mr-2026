"""Prepare candidate colocalization inputs and FinnGen precomputed coloc evidence.

The local AF/HF files can support lead-SNP checks, but formal coloc with the
primary Nielsen/HERMES outcomes still requires dense regional pQTL summary
statistics. As an immediately available complement, this script also queries
FinnGen R12 PheWeb precomputed pQTL-disease colocalization records.
"""

from __future__ import annotations

import csv
import json
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[1]
EXPOSURE = PROJECT_ROOT / "data" / "processed" / "ukbppp_inflammation_cis_pqtl_instruments_lead_no_mhc.csv"
PQTL = PROJECT_ROOT / "data" / "processed" / "ukbppp_inflammation_panel_cis_pqtl_p5e-8_1mb.csv"
OUTCOMES = {
    "AF": PROJECT_ROOT / "data" / "outcomes" / "AF_outcome_effects_local.csv",
    "HF": PROJECT_ROOT / "data" / "outcomes" / "HF_outcome_effects_local.csv",
}
OUT_DIR = PROJECT_ROOT / "results" / "coloc_inputs"
CANDIDATES = {"FGF5", "LPA"}
FINNGEN_PHENOS = {"I9_AF", "I9_HEARTFAIL"}
FINNGEN_GENE_PQTL_COLOC = "https://r12.finngen.fi/api/gene_pqtl_colocalization/{gene}"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def get_json(url: str) -> Any:
    request = urllib.request.Request(url, headers={"User-Agent": "cardio-no-lab-project/0.1"})
    with urllib.request.urlopen(request, timeout=60) as response:
        return json.loads(response.read().decode("utf-8"))


def flatten_finngen_colocs(gene: str) -> list[dict[str, Any]]:
    url = FINNGEN_GENE_PQTL_COLOC.format(gene=gene)
    try:
        payload = get_json(url)
    except (urllib.error.URLError, TimeoutError):
        return []

    rows: list[dict[str, Any]] = []
    for pqtl in payload:
        for group in pqtl.get("disease_colocalizations") or []:
            for coloc in group:
                if coloc.get("phenotype1") not in FINNGEN_PHENOS:
                    continue
                if coloc.get("phenotype2_description") != gene and coloc.get("phenotype2") != gene:
                    continue
                rows.append(
                    {
                        "protein": gene,
                        "finnGen_disease": coloc.get("phenotype1"),
                        "disease_description": coloc.get("phenotype1_description"),
                        "pqtl_trait": coloc.get("phenotype2"),
                        "pqtl_description": coloc.get("phenotype2_description"),
                        "pqtl_source": coloc.get("source2"),
                        "pqtl_source_display": coloc.get("source2_displayname"),
                        "disease_variant": (
                            f"{coloc.get('locus_id1_chromosome')}:{coloc.get('locus_id1_position')}:"
                            f"{coloc.get('locus_id1_ref')}:{coloc.get('locus_id1_alt')}"
                        ),
                        "pqtl_variant": (
                            f"{coloc.get('locus_id2_chromosome')}:{coloc.get('locus_id2_position')}:"
                            f"{coloc.get('locus_id2_ref')}:{coloc.get('locus_id2_alt')}"
                        ),
                        "beta_disease": coloc.get("beta1"),
                        "beta_pqtl": coloc.get("beta2"),
                        "p_disease": coloc.get("pval1"),
                        "p_pqtl": coloc.get("pval2"),
                        "clpp": coloc.get("clpp"),
                        "clpa": coloc.get("clpa"),
                        "len_cs_disease": coloc.get("len_cs1"),
                        "len_cs_pqtl": coloc.get("len_cs2"),
                        "len_intersection": coloc.get("len_inter"),
                        "api_source": url,
                    }
                )
    return rows


def main() -> None:
    exposure = {row["protein"]: row for row in read_csv(EXPOSURE) if row["protein"] in CANDIDATES}
    pqtl = [row for row in read_csv(PQTL) if row["target_gene"] in CANDIDATES]
    coloc_rows: list[dict[str, Any]] = []

    lines = [
        "# Candidate Colocalization Input Status",
        "",
        "This folder contains the currently available candidate pQTL rows, matched primary AF/HF lead outcome rows, and FinnGen R12 precomputed pQTL-disease colocalization records where available.",
        "",
        "Important: this is not yet a full formal coloc run for the primary Nielsen AF and HERMES HF outcomes. Formal local coloc still requires dense regional UKB-PPP pQTL summary statistics, ideally all SNPs in +/-1 Mb around each candidate gene. The current UKB-PPP browser-derived pQTL table only contains P < 5e-8 cis-pQTL rows from the previous screening step.",
        "",
        "FinnGen PheWeb reports CLPP/CLPA from fine-mapped pQTL-disease colocalization. These values are useful supporting evidence, but they are not the same output scale as coloc PP.H4 from a locally run coloc.abf analysis.",
        "",
    ]
    for protein, exp in exposure.items():
        protein_pqtl = [row for row in pqtl if row["target_gene"] == protein]
        pqtl_fields = list(protein_pqtl[0].keys()) if protein_pqtl else ["target_gene"]
        write_csv(OUT_DIR / f"{protein}_ukbppp_candidate_pqtl_available.csv", protein_pqtl, pqtl_fields)
        lines.append(f"## {protein}")
        lines.append("")
        lines.append(f"- Lead pQTL: {exp['SNP']} ({exp['effect_allele']}/{exp['other_allele']}), beta={exp['beta']}, P={exp['pval']}")
        lines.append(f"- Available pQTL rows in current screening table: {len(protein_pqtl)}")
        for outcome_id, path in OUTCOMES.items():
            rows = [row for row in read_csv(path) if row["SNP"] == exp["SNP"]]
            if rows:
                write_csv(OUT_DIR / f"{protein}_{outcome_id}_lead_outcome_row.csv", rows, list(rows[0].keys()))
                row = rows[0]
                lines.append(
                    f"- {outcome_id} lead outcome row: beta={row['beta']}, SE={row['se']}, P={row['pval']}, alleles={row['effect_allele']}/{row['other_allele']}"
                )
            else:
                lines.append(f"- {outcome_id}: lead SNP not found in local outcome extract.")

        protein_colocs = flatten_finngen_colocs(protein)
        coloc_rows.extend(protein_colocs)
        if protein_colocs:
            lines.append("- FinnGen precomputed pQTL-disease colocalization rows:")
            for row in protein_colocs:
                lines.append(
                    f"  - {row['finnGen_disease']} with {row['pqtl_source_display']}: "
                    f"CLPP={row['clpp']}, CLPA={row['clpa']}, disease_variant={row['disease_variant']}, pQTL_variant={row['pqtl_variant']}"
                )
        else:
            lines.append("- FinnGen precomputed pQTL-disease colocalization rows for I9_AF/I9_HEARTFAIL: 0")
        lines.append("")

    coloc_fields = [
        "protein",
        "finnGen_disease",
        "disease_description",
        "pqtl_trait",
        "pqtl_description",
        "pqtl_source",
        "pqtl_source_display",
        "disease_variant",
        "pqtl_variant",
        "beta_disease",
        "beta_pqtl",
        "p_disease",
        "p_pqtl",
        "clpp",
        "clpa",
        "len_cs_disease",
        "len_cs_pqtl",
        "len_intersection",
        "api_source",
    ]
    write_csv(OUT_DIR / "finngen_pheweb_candidate_pqtl_disease_coloc.csv", coloc_rows, coloc_fields)

    lines += [
        "Next required input for formal primary-outcome coloc:",
        "- UKB-PPP dense regional pQTL summary statistics for FGF5 and LPA.",
        "- AF/HF dense regional GWAS extracts can be prepared from the downloaded local summary statistics once pQTL regional data are available.",
        "",
        f"FinnGen precomputed coloc table: `{OUT_DIR / 'finngen_pheweb_candidate_pqtl_disease_coloc.csv'}`",
    ]
    (OUT_DIR / "candidate_coloc_input_status.md").write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote coloc input status to {OUT_DIR / 'candidate_coloc_input_status.md'}")


if __name__ == "__main__":
    main()
