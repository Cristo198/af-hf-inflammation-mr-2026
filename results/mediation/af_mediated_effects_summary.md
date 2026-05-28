# AF-Mediated Protein-HF Analysis Summary

Path tested: genetically predicted inflammatory protein level -> AF -> HF.

## AF -> HF MR

- AF instruments: 11650 genome-wide significant variants; 75 after 10,000,000 bp distance pruning.
- Harmonised AF instruments with HF GWAS: 73 (HF matched: 74; ambiguous palindromic skipped: 1; allele mismatches: 0).
- AF -> HF IVW estimate: OR 1.2653 (95% CI 1.233-1.2984), beta 0.2353, SE 0.0132, P=3.57e-71.
- Cochran-style Q for the AF -> HF instrument set: 157.2948.

## Candidate Mediation Results

| Protein | Total protein -> HF OR | Indirect OR via AF | Indirect P | Proportion mediated | Interpretation |
|---|---:|---:|---:|---:|---|
| FGF5 | 1.0293 | 1.0145 (1.0093-1.0197) | 3.87e-08 | 0.4969 | AF-mediated pathway supported |
| LPA | 1.0797 | 1.0081 (1.0022-1.014) | 0.0068 | 0.1048 | AF-mediated pathway supported |

## Scope

- All proteins with both AF and HF primary MR estimates: 408.
- AF-HF nominal-overlap proteins: 6.
- FGF5/LPA shared candidates: 2.

## Caveats

- The AF -> HF step uses distance-pruned genome-wide significant AF variants because a formal LD reference panel is not available in this workspace.
- Product-of-coefficients mediation is reported on the log-odds scale and should be treated as exploratory for binary outcomes.
- The indirect effect assumes the protein -> AF and AF -> HF estimates are approximately independent and does not replace multivariable MR with individual-level data.

AF-HF instrument table: `D:\Backup\Documents\New project 2\cardio_no_lab_paper_project\results\mediation\af_to_hf_harmonised_instruments.csv`
All mediation table: `D:\Backup\Documents\New project 2\cardio_no_lab_paper_project\results\mediation\af_mediated_effects_all_primary.csv`
Candidate mediation table: `D:\Backup\Documents\New project 2\cardio_no_lab_paper_project\results\mediation\af_mediated_effects_fgf5_lpa.csv`
