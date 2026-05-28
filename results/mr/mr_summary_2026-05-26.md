# Preliminary Wald-ratio MR Summary

Date: 2026-05-26

Input:
- Exposure: 529 UKB-PPP inflammation-panel lead cis-pQTL instruments after MHC exclusion.
- Outcomes: Nielsen 2018 AF and HERMES 2019 HF local summary statistics.
- Method: single-variant Wald ratio; FDR by Benjamini-Hochberg within outcome.

Coverage:
- AF outcome effects extracted: 489 harmonised MR rows.
- HF outcome effects extracted: 410 harmonised MR rows.

Significance counts:
- AF: nominal P<0.05 = 37, FDR<0.05 = 5
- HF: nominal P<0.05 = 41, FDR<0.05 = 4

FDR-significant AF results:
- FGF5 (Fibroblast growth factor 5), OR 1.0628884 [1.0411235, 1.0851083], P=7.5693691e-09, FDR=3.70142e-06
- TNFSF12 (Tumor necrosis factor ligand superfamily member 12), OR 0.89624508 [0.8615618, 0.93232458], P=5.3280571e-08, FDR=1.30271e-05
- SPON1 (Spondin-1), OR 1.1137385 [1.0653928, 1.1642781], P=1.9593173e-06, FDR=0.000319369
- PKLR (Pyruvate kinase PKLR), OR 0.86176886 [0.79999741, 0.92830996], P=8.8442529e-05, FDR=0.0108121
- NFATC1 (Nuclear factor of activated T-cells, cytoplasmic 1), OR 1.2774074 [1.119434, 1.4576739], P=0.00027783496, FDR=0.0271723

FDR-significant HF results:
- CELSR2 (Cadherin EGF LAG seven-pass G-type receptor 2), OR 0.88273185 [0.84700552, 0.91996511], P=3.2682066e-09, FDR=1.33996e-06
- LPA (Apolipoprotein(a)), OR 1.0797214 [1.0495757, 1.1107329], P=1.1017337e-07, FDR=2.25855e-05
- ABO (Histo-blood group ABO system transferase), OR 1.0343982 [1.0197848, 1.0492211], P=3.1802168e-06, FDR=0.00043463
- APOA2 (Apolipoprotein A-II), OR 0.64252639 [0.50722419, 0.81392048], P=0.00024573278, FDR=0.0251876

Preliminary shared candidates:
- FGF5 (Fibroblast growth factor 5): AF OR 1.0628884 [1.0411235, 1.0851083], P=7.5693691e-09, FDR=3.70142e-06; HF OR 1.0293035 [1.0044669, 1.0547541], P=0.02046851, FDR=0.441689; AF_FDR+HF_nominal, same_risk.
- LPA (Apolipoprotein(a)): AF OR 1.0347384 [1.0097401, 1.0603555], P=0.0062032059, FDR=0.216669; HF OR 1.0797214 [1.0495757, 1.1107329], P=1.1017337e-07, FDR=2.25855e-05; HF_FDR+AF_nominal, same_risk.

Important interpretation note:
These are first-pass Wald-ratio results. Shared candidates require replication, allele-coding audit, LD/proxy review, and colocalization before being interpreted as causal targets.

FDR table: `D:\Backup\Documents\New project 2\cardio_no_lab_paper_project\results\mr\fdr_significant_preliminary.csv`
Shared candidate table: `D:\Backup\Documents\New project 2\cardio_no_lab_paper_project\results\mr\shared_candidate_preliminary.csv`