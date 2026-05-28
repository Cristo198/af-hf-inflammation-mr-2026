# Genetically Supported Inflammatory Mediators of the AF-HF Continuum

This repository contains the reproducible analysis code, derived non-restricted tables, figure files, and manuscript-supporting materials for:

**Genetically Supported Inflammatory Mediators of the Atrial Fibrillation-Heart Failure Continuum: cis-pQTL Mendelian Randomization, Dual-Outcome Colocalization, and Mediation-Based Target Prioritization**

Repository: https://github.com/Cristo198/af-hf-inflammation-mr-2026

## Authors

Yangfeng Qin, Langping Township Health Center, Bobai County.

Runchuan Feng, Nanning Fourth People's Hospital.

Yangfeng Qin and Runchuan Feng contributed equally. Yangfeng Qin is the corresponding author for the manuscript.

## Overview

The project evaluates circulating inflammatory proteins as genetically supported mediators of the atrial fibrillation (AF)-heart failure (HF) continuum. The workflow integrates:

1. UKB-PPP Olink inflammatory-panel cis-pQTL instrument selection.
2. Two-sample cis-pQTL Mendelian randomization for AF and HF.
3. Full-panel FinnGen R12 replication.
4. Supplementary UKB/OpenGWAS replication.
5. Formal coloc.abf-style regional colocalization for FGF5 and LPA.
6. Feasible sensitivity checks, candidate reverse MR, and AF-mediated effect analysis.
7. Evidence-based target prioritization.

## Repository Contents

- `scripts/`: custom Python scripts for exposure construction, outcome extraction, MR, colocalization, replication, sensitivity analyses, reverse MR, mediation analysis, and prioritization.
- `tables/`: derived main and supplementary tables.
- `results/`: derived summaries and figure files.
- `DATA_ACCESS.md`: external data-source links and access notes.
- `requirements.txt`: minimal software dependency note.
- `CITATION.cff` and `.zenodo.json`: citation and Zenodo metadata templates.

## Data Availability

This repository does not redistribute raw or access-controlled GWAS/pQTL summary-statistics files. Users should obtain the required raw inputs from the original data providers listed in `DATA_ACCESS.md`.

The public release includes derived non-restricted result tables and figures generated for the manuscript. If a data provider's terms impose stricter restrictions on derived SNP-level tables, users should follow the original provider's terms.

## Software Environment

The completed local analysis used Python 3.12.4 on Windows 10.0.19045. Core analysis scripts relied on Python standard-library modules. `synapseclient` 4.12.0 was used only to assist controlled UKB-PPP pGWAS downloads.

See `tables/supplementary_table_s22_software_environment.md` for the frozen environment table.

## Reproducibility Notes

The scripts are numbered in the approximate order used during the analysis. Some steps require raw external inputs that are not included in this repository because of data-provider terms.

Before rerunning the pipeline, download or request the raw datasets described in `DATA_ACCESS.md`, place them in the expected local folders, and review file names in the scripts.

## Citation

If using this repository, please cite the manuscript and the archived Zenodo DOI once available. The code repository is available at https://github.com/Cristo198/af-hf-inflammation-mr-2026.

## License

Code is released under the MIT License unless otherwise stated. Derived tables and manuscript-supporting documentation should be cited with the associated manuscript and Zenodo record.
