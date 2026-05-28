# FinnGen R12 Full-Panel Replication Summary

Scope: all 529 lead cis-pQTL instruments from the UKB-PPP inflammation panel main analysis.

## Extraction Coverage

| FinnGen outcome | Harmonised instruments | Position rows seen | Allele mismatches | Gzip EOF warning |
|---|---:|---:|---:|---:|
| I9_AF | 496 | 507 | 5 | 0 |
| I9_HEARTFAIL | 496 | 507 | 5 | 0 |

## Replication Counts

| Comparison | Count |
|---|---:|
| AF primary proteins matched in FinnGen | 467 |
| AF same-direction nominal FinnGen replication | 50 |
| AF same-direction FDR FinnGen replication | 14 |
| HF primary proteins matched in FinnGen | 399 |
| HF same-direction nominal FinnGen replication | 24 |
| HF same-direction FDR FinnGen replication | 1 |

## Primary FDR Signals

| Protein | Primary outcome | Primary OR/P/FDR | FinnGen OR/P/FDR | Replication status |
|---|---|---|---|---|
| FGF5 | AF | 1.0629/7.57e-09/3.70e-06 | 1.1297/0/0 | same_direction_fdr |
| NFATC1 | AF | 1.2774/2.78e-04/0.0272 | 1.2603/0.0037/0.084 | same_direction_nominal |
| PKLR | AF | 0.8618/8.84e-05/0.0108 | 0.9522/0.2496/0.6911 | same_direction_not_significant |
| SPON1 | AF | 1.1137/1.96e-06/3.19e-04 | 1.1438/7.31e-06/9.07e-04 | same_direction_fdr |
| TNFSF12 | AF | 0.8962/5.33e-08/1.30e-05 | 0.9762/0.3701/0.7683 | same_direction_not_significant |
| ABO | HF | 1.0344/3.18e-06/4.35e-04 | 1.012/0.0887/0.524 | same_direction_not_significant |
| APOA2 | HF | 0.6425/2.46e-04/0.0252 | 1.1286/0.3323/0.7886 | opposite_direction_not_significant |
| CELSR2 | HF | 0.8827/3.27e-09/1.34e-06 | 0.9772/0.2841/0.7226 | same_direction_not_significant |
| LPA | HF | 1.0797/1.10e-07/2.26e-05 | 1.0302/0.0318/0.4569 | same_direction_nominal |

## FGF5/LPA Shared Candidates

| Protein | Primary outcome | FinnGen outcome | Primary OR/P | FinnGen OR/P | Replication status |
|---|---|---|---|---|---|
| FGF5 | AF | FG_AF | 1.0629/7.57e-09 | 1.1297/0 | same_direction_fdr |
| FGF5 | HF | FG_HF | 1.0293/0.0205 | 1.0236/0.0712 | same_direction_not_significant |
| LPA | AF | FG_AF | 1.0347/0.0062 | 1.0655/1.49e-05 | same_direction_fdr |
| LPA | HF | FG_HF | 1.0797/1.10e-07 | 1.0302/0.0318 | same_direction_nominal |

Interpretation: FinnGen R12 AF and HF full-panel replication is now complete for the 529 lead cis-pQTL instruments. UKB/OpenGWAS replication was completed separately as a low-coverage supplementary analysis and did not cover FGF5/LPA.

Full MR table: `D:\Backup\Documents\New project 2\cardio_no_lab_paper_project\results\replication\full_panel\finngen_r12_full_panel_wald_mr.csv`
Primary comparison table: `D:\Backup\Documents\New project 2\cardio_no_lab_paper_project\results\replication\full_panel\finngen_r12_primary_comparison.csv`
