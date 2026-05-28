# Instrument Quality Control Summary

## Exposure instruments

- Lead cis-pQTL instruments after MHC exclusion: 529
- Panel counts: Inflammation=270, Inflammation_II=259
- Weak instruments by F <= 10: 0
- F statistic: min=45.797, Q1=227.916, median=788.518, Q3=2958.653, max=2.239e+04
- cis gene match equals target protein gene: 526/529
- Palindromic allele instruments: 78/529
- Locus annotation flags: target_or_blank=381, nearest_locus_not_target_review=148

## Candidate instrument strength

| Protein | SNP | Effect allele | Other allele | EAF | Beta | SE | F statistic | Locus flag |
|---|---|---|---|---|---|---|---|---|
| FGF5 | rs12509595 | C | T | 0.293 | 0.682 | 0.008 | 7582.335 | target_or_blank |
| LPA | rs56393506 | T | C | 0.172 | 0.761 | 0.01 | 5377.806 | target_or_blank |

## Outcome extraction and harmonisation

- AF: extracted 492/529 SNPs; harmonised MR rows 489; not harmonised after extraction 3; aligned=201, flipped=285, aligned_complement=1, flipped_complement=2.
- HF: extracted 418/529 SNPs; harmonised MR rows 410; not harmonised after extraction 8; aligned=240, flipped=169, aligned_complement=1, flipped_complement=0.

## Pleiotropy note

Because the main design uses one lead cis-pQTL per protein, MR-Egger, weighted median, MR-PRESSO and leave-one-out are not statistically applicable to the primary single-variant estimates. The practical pleiotropy controls at this stage are cis restriction, MHC exclusion, strong-instrument filtering, allele harmonisation audit, replication, and formal colocalization. Formal colocalization remains the key next check for distinguishing a shared causal variant from LD-driven association.

Output tables:
- `D:\Backup\Documents\New project 2\cardio_no_lab_paper_project\results\qc\exposure_instrument_qc.csv`
- `D:\Backup\Documents\New project 2\cardio_no_lab_paper_project\results\qc\outcome_harmonisation_qc.csv`
- `D:\Backup\Documents\New project 2\cardio_no_lab_paper_project\results\qc\candidate_instrument_qc.csv`