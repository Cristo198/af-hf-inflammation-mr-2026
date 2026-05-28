# Supplementary Table S22. Frozen Software Environment

Environment captured on 2026-05-28 in `D:\Backup\Documents\New project 2\cardio_no_lab_paper_project`.

| Component | Version or status | Role in project | Notes |
|---|---|---|---|
| Operating system | Windows 10.0.19045 | Local analysis platform | Captured in the project workspace. |
| Python | 3.12.4 | Primary analysis runtime | Executable: `C:\Users\Administrator\AppData\Local\Programs\Python\Python312\python.exe`. |
| Python standard library | `argparse`; `collections`; `csv`; `gzip`; `html`; `importlib`; `io`; `json`; `math`; `os`; `pathlib`; `re`; `ssl`; `statistics`; `subprocess`; `sys`; `tarfile`; `time`; `typing`; `urllib`; `zipfile`; `xml.etree.ElementTree` | Core data parsing, harmonization, MR, colocalization, replication, sensitivity, reverse MR, mediation, prioritization, and table generation | No external numerical-analysis package was required for the executed core pipeline. |
| `synapseclient` | 4.12.0 | Controlled UKB-PPP pGWAS download helper | Used for Synapse-controlled data download only; not required for the offline analytical scripts after data were downloaded. |
| `pandas` | not installed | Not used | Core scripts used `csv` and other standard-library parsing instead. |
| `numpy` | not installed | Not used | Core scripts used `math`/`statistics`/standard-library implementations. |
| `scipy` | not installed | Not used | Normal approximation and coloc calculations were implemented in project scripts. |
| `matplotlib` | not installed | Not used | Project figures were generated as script-written SVG/HTML/CSV artifacts rather than by `matplotlib`. |
| `seaborn` | not installed | Not used | Not required for current figures. |
| `statsmodels` | not installed | Not used | Not required for current MR or mediation calculations. |
| `openpyxl` | not installed | Not used | No Excel-writing dependency was used in the current pipeline. |
| R | not used | Not used | No R scripts or R package dependencies were required for the completed analyses. |

Recommended final submission note: the public code repository should include this table, all project scripts, and the exact input-file manifest. Raw restricted UKB-PPP pGWAS files should be accessed through UKB-PPP/Synapse rather than redistributed.
