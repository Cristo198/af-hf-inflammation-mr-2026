# Candidate Reverse MR Summary

Disease instruments were selected at P < 5e-08 from the local AF and HF GWAS files and distance-pruned at 10,000,000 bp because an LD reference panel is not available in this workspace.

Two result sets are reported: all distance-pruned disease instruments, and a more conservative target-cis-excluded set that removes disease instruments within +/-1 Mb of the target protein cis locus.

| Analysis set | Exposure | Outcome protein | Method | Instruments | Cis excluded | Beta | SE | P | Interpretation |
|---|---|---|---|---:|---:|---:|---:|---:|---|
| all_distance_pruned | AF genetic liability | FGF5 | Fixed-effect IVW | 74 | 0 | 0.0807 | 0.0109 | 1.40e-13 | nominal_reverse_signal |
| target_cis_excluded | AF genetic liability | FGF5 | Fixed-effect IVW | 73 | 1 | 0.0241 | 0.0109 | 0.0278 | nominal_reverse_signal |
| all_distance_pruned | HF genetic liability | FGF5 | Fixed-effect IVW | 9 | 0 | 0.051 | 0.0423 | 0.2286 | no_reverse_signal |
| target_cis_excluded | HF genetic liability | FGF5 | Fixed-effect IVW | 9 | 0 | 0.051 | 0.0423 | 0.2286 | no_reverse_signal |
| all_distance_pruned | AF genetic liability | LPA | Fixed-effect IVW | 74 | 0 | -0.0051 | 0.0102 | 0.6154 | no_reverse_signal |
| target_cis_excluded | AF genetic liability | LPA | Fixed-effect IVW | 74 | 0 | -0.0051 | 0.0102 | 0.6154 | no_reverse_signal |
| all_distance_pruned | HF genetic liability | LPA | Fixed-effect IVW | 9 | 0 | 0.9249 | 0.0415 | 0.00e+00 | nominal_reverse_signal |
| target_cis_excluded | HF genetic liability | LPA | Fixed-effect IVW | 8 | 1 | -0.0066 | 0.0437 | 0.8795 | no_reverse_signal |

Caveat: this is a candidate-level reverse MR using distance-pruned disease instruments and the available UKB-PPP pGWAS files for FGF5 and LPA. It should be treated as a screening check rather than a definitive bidirectional MR.

Instrument summary: `D:\Backup\Documents\New project 2\cardio_no_lab_paper_project\results\reverse_mr\disease_instrument_summary.csv`
Harmonised instruments: `D:\Backup\Documents\New project 2\cardio_no_lab_paper_project\results\reverse_mr\candidate_reverse_mr_harmonised_instruments.csv`
Reverse MR results: `D:\Backup\Documents\New project 2\cardio_no_lab_paper_project\results\reverse_mr\candidate_reverse_mr_results.csv`
