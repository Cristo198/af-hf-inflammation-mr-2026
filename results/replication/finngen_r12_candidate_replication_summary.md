# FinnGen R12 Candidate Replication Summary

Source: FinnGen R12 public PheWeb variant API. FinnGen effect allele is the alternative allele; candidate variant IDs were formed as `chr-pos-other-effect`, so the FinnGen alternative allele is aligned to the UKB-PPP exposure effect allele for FGF5 and LPA.

| Protein | FinnGen endpoint | Variant | OR (95% CI) | P | Direction |
|---|---|---|---|---|---|
| FGF5 | I9_AF | 4:80261400:T:C | 1.1297251 (1.1000495-1.1602012) | 2.6804024e-19 | same_risk |
| FGF5 | I9_HEARTFAIL | 4:80261400:T:C | NA | NA | missing_effect |
| LPA | I9_AF | 6:160668275:C:T | 1.0654534 (1.0353081-1.0964765) | 1.4940323e-05 | same_risk |
| LPA | I9_HEARTFAIL | 6:160668275:C:T | 1.0301798 (1.0025919-1.0585267) | 0.031799479 | same_risk |

Interpretation:
- FGF5 replicated for FinnGen AF with the same risk-increasing direction; FinnGen strict HF does not provide an exact-variant beta through the public variant API.
- LPA replicated in the same risk-increasing direction for FinnGen AF and FinnGen strict HF; HF is nominally significant.

Effect table: `D:\Backup\Documents\New project 2\cardio_no_lab_paper_project\data\replication\finngen_r12_candidate_pheweb_effects.csv`
MR table: `D:\Backup\Documents\New project 2\cardio_no_lab_paper_project\results\replication\finngen_r12_candidate_wald_mr.csv`