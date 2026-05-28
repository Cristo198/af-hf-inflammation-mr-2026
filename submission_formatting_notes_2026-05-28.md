# Submission Formatting Notes

Date: 2026-05-28

Target journal: not yet selected. Current formatting is prepared as a conservative biomedical/JCR Q1 working package and should be adapted after the journal is chosen.

## Reference Format

Current working style: numbered in-text references, with the reference seed library stored in `references_seed.md`.

Before submission:

1. Import all references into Zotero, EndNote, NoteExpress, or another reference manager.
2. Export references in the selected journal style.
3. Verify author order, article title, journal abbreviation, year, volume, issue, pages, DOI, PMID, PMCID, and data-access URLs.
4. Replace the working reference note in the manuscript with the formatted reference list.

Priority references that should remain in the final library include STROBE-MR, coloc.abf, UKB-PPP, Nielsen AF GWAS, HERMES HF GWAS, FinnGen, OpenGWAS/MR-Base, AF/HF guidelines, AF inflammation mechanisms, HF inflammation mechanisms, and adjacent AF/HF proteomic MR studies.

## Main Tables

Prepared main-text table:

| Main table | File | Status |
|---|---|---|
| Table 1. Data sources, sample sizes, ancestry, phenotype definitions, local files, and access links | `tables/table1_data_sources.md`; `tables/table1_data_sources.csv` | Prepared |

Planned main-text tables already listed in the manuscript:

| Main table | Content | Current source |
|---|---|---|
| Table 2 | Lead cis-pQTL instruments | `tables/supplementary_table_s1_exposure_lead_cis_pqtl_qc.csv` |
| Table 3 | Primary MR results | `tables/supplementary_table_s3_all_primary_wald_mr_results.csv`; `tables/supplementary_table_s4_fdr_significant_primary_mr_results.csv` |
| Table 4 | Replication and sensitivity analyses | `tables/supplementary_table_s13_proxy_steiger_directionality.csv`; `tables/supplementary_table_s16_finngen_r12_full_panel_wald_mr.csv`; `tables/supplementary_table_s18_ukb_opengwas_full_panel_wald_mr.csv` |
| Table 5 | Formal colocalization and AF-mediated effect estimates | `tables/supplementary_table_s11_formal_coloc_abf_summary.csv`; `tables/supplementary_table_s21_af_mediated_effects.csv` |
| Table 6 | Target-prioritization scorecard | `tables/supplementary_table_s12_target_priority_scorecard.csv` |

## Supplementary Tables

Supplementary table manifest: `tables/supplementary_tables_manifest.md` and `tables/supplementary_tables_manifest.csv`.

Newly added: Supplementary Table S22, frozen software environment and package status.

## Figure Package

Current generated figure files:

| Figure file | Size bytes | Current role |
|---|---:|---|
| `results/figures/af_wald_mr_volcano.svg` | 39360 | Primary AF MR volcano plot |
| `results/figures/hf_wald_mr_volcano.svg` | 33436 | Primary HF MR volcano plot |
| `results/figures/fdr_significant_wald_mr_forest.svg` | 4939 | FDR-significant primary MR forest plot |
| `results/figures/fgf5_lpa_primary_replication_forest.svg` | 3623 | FGF5/LPA primary and replication forest plot |
| `results/figures/af_hf_nominal_overlap_venn.svg` | 1873 | AF/HF nominal overlap Venn plot |
| `results/figures/mr_significance_upset.svg` | 5176 | AF/HF MR significance UpSet plot |

Before submission:

1. Convert SVG figures to the target journal's required format, usually PDF/EPS/TIFF/PNG.
2. Use at least 300 dpi for color raster figures and 600 dpi for line art if TIFF/PNG is required.
3. Ensure all labels, legends, and axis text remain readable at final journal column width.
4. Add final regional colocalization plots if the target journal expects locus-level visualization.
5. Prepare a separate figure legend file if required by the journal.

## Submission-Level Items Still Pending

1. Target journal selection.
2. Final corresponding-author email.
3. Author confirmation of funding and conflict-of-interest statements.
4. Final reference-manager export in the selected journal style.
5. Final figure conversion and journal-specific table placement.
