# Sensitivity Analysis Summary

This project uses one lead cis-pQTL per protein for the primary analysis. Therefore MR-Egger, weighted median, MR-PRESSO and leave-one-out are not statistically applicable to the main single-variant estimates.

Feasible sensitivity checks completed in this step:

- Strong-instrument check: 0 of 899 harmonised primary MR rows had F <= 10.
- Allele ambiguity check: 141 of 899 harmonised primary MR rows used palindromic alleles; FGF5 and LPA candidate instruments are not palindromic.
- Proxy-Steiger directionality: 896 of 899 harmonised primary MR rows had exposure R2 proxy greater than outcome R2 proxy.
- Candidate proxy-Steiger directionality: 4 of 4 FGF5/LPA candidate rows supported the protein-to-disease direction.
- LD/confounding sensitivity: formal coloc.abf was completed for FGF5/LPA against AF and HF primary outcomes; FGF5-AF was strongly supported, whereas FGF5-HF, LPA-AF and LPA-HF were not strongly supported.

Important caveat: proxy-Steiger values here are approximate because binary outcome variance is represented on the observed log-odds scale and not transformed to liability-scale R2.

All-row table: `D:\Backup\Documents\New project 2\cardio_no_lab_paper_project\results\sensitivity\proxy_steiger_directionality_all_primary_mr.csv`
Candidate table: `D:\Backup\Documents\New project 2\cardio_no_lab_paper_project\results\sensitivity\candidate_sensitivity_summary.csv`
