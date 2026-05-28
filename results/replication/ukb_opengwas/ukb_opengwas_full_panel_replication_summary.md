# UKB/OpenGWAS Full-Panel Replication Summary

Scope: all 529 lead cis-pQTL instruments from the UKB-PPP inflammation panel main analysis. Because the exposure data are from UKB-PPP, these analyses are treated as supplementary and may be affected by sample overlap.

## Extraction Coverage

| UKB/OpenGWAS outcome | Harmonised instruments | Position rows seen | Allele mismatches | Gzip EOF warning |
|---|---:|---:|---:|---:|
| ukb-b-964 | 7 | 7 | 0 | 0 |
| ukb-d-HEARTFAIL | 7 | 8 | 1 | 0 |

## Replication Counts

| Comparison | Count |
|---|---:|
| AF primary proteins matched in UKB/OpenGWAS | 7 |
| AF same-direction nominal UKB/OpenGWAS replication | 0 |
| AF same-direction FDR UKB/OpenGWAS replication | 0 |
| HF primary proteins matched in UKB/OpenGWAS | 7 |
| HF same-direction nominal UKB/OpenGWAS replication | 0 |
| HF same-direction FDR UKB/OpenGWAS replication | 0 |

## Primary FDR Signals

| Protein | Primary outcome | Primary OR/P/FDR | UKB OR/P/FDR | Replication status |
|---|---|---|---|---|
| FGF5 | AF | 1.0629/7.57e-09/3.70e-06 | NA/NA/NA | not_matched_in_ukb |
| NFATC1 | AF | 1.2774/2.78e-04/0.0272 | NA/NA/NA | not_matched_in_ukb |
| PKLR | AF | 0.8618/8.84e-05/0.0108 | NA/NA/NA | not_matched_in_ukb |
| SPON1 | AF | 1.1137/1.96e-06/3.19e-04 | NA/NA/NA | not_matched_in_ukb |
| TNFSF12 | AF | 0.8962/5.33e-08/1.30e-05 | NA/NA/NA | not_matched_in_ukb |
| ABO | HF | 1.0344/3.18e-06/4.35e-04 | NA/NA/NA | not_matched_in_ukb |
| APOA2 | HF | 0.6425/2.46e-04/0.0252 | NA/NA/NA | not_matched_in_ukb |
| CELSR2 | HF | 0.8827/3.27e-09/1.34e-06 | NA/NA/NA | not_matched_in_ukb |
| LPA | HF | 1.0797/1.10e-07/2.26e-05 | NA/NA/NA | not_matched_in_ukb |

## FGF5/LPA Shared Candidates

| Protein | Primary outcome | UKB outcome | Primary OR/P | UKB OR/P | Replication status |
|---|---|---|---|---|---|
| FGF5 | AF | UKB_AF | 1.0629/7.57e-09 | NA/NA | not_matched_in_ukb |
| FGF5 | HF | UKB_HF | 1.0293/0.0205 | NA/NA | not_matched_in_ukb |
| LPA | AF | UKB_AF | 1.0347/0.0062 | NA/NA | not_matched_in_ukb |
| LPA | HF | UKB_HF | 1.0797/1.10e-07 | NA/NA | not_matched_in_ukb |

Interpretation: UKB/OpenGWAS replication is completed as a supplementary check, but coverage is very low for the 529 lead cis-pQTL panel and the FGF5/LPA shared candidates were not matched. Because UKB-PPP exposure data and UKB outcome data may overlap, these results should not be interpreted as fully independent replication.

Full MR table: `D:\Backup\Documents\New project 2\cardio_no_lab_paper_project\results\replication\ukb_opengwas\ukb_opengwas_full_panel_wald_mr.csv`
Primary comparison table: `D:\Backup\Documents\New project 2\cardio_no_lab_paper_project\results\replication\ukb_opengwas\ukb_opengwas_primary_comparison.csv`
