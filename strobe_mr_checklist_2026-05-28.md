# STROBE-MR Checklist Review

Date: 2026-05-28

Guideline source: STROBE-MR statement and checklist by Skrivankova et al. [1]. Official checklist source: https://www.strobe-mr.org/download/strobe-mr-checklist/

Manuscript reviewed: `manuscript_draft_en_v0_1_2026-05-28.md`

Status key: Addressed = covered in the current English draft; Partial = covered but needs author/data confirmation or formatting; Pending = requires author input or final submission preparation.

| Item | STROBE-MR domain | Current status | Current manuscript location | Action before submission |
|---:|---|---|---|---|
| 1 | Title and abstract identify MR design | Addressed | Title; Abstract Methods | Keep "Mendelian randomization" and "cis-pQTL" in title/abstract. |
| 2 | Scientific background and rationale | Addressed | Introduction paragraphs 1-4 | None beyond final reference formatting. |
| 3 | Objectives and causal hypothesis | Addressed | Introduction final paragraph; Abstract Objectives | Keep hypothesis language cautious and assumption-dependent. |
| 4a | Study design and setting/data sources | Addressed | Methods: Study Design; Outcome GWAS Data; Exposure Data; Table 1 | Table 1 prepared with source, ancestry, sample size, phenotype definition, local files, and access links. |
| 4b | Participants and sample sizes | Addressed | Methods: Outcome GWAS Data; Results: Exposure Instruments and Outcome Coverage; Table 1 | Keep source-specific notes for FinnGen and UKB/OpenGWAS during final copyediting. |
| 4c | Genetic variant measurement, QC, and selection | Addressed | Methods: Exposure Data and Genetic Instruments; Allele Harmonization | None; optional add flow diagram. |
| 4d | Exposure/outcome definitions and diagnostic criteria | Addressed | Methods: Outcome GWAS Data; Table 1 | None beyond target-journal formatting. |
| 4e | Ethics and informed consent | Addressed | Methods: Ethics; Ethics Statement | Confirm whether local institutional exemption statement is required. |
| 5 | Core instrumental-variable assumptions | Addressed | Methods: Study Design | Keep relevance, independence, and exclusion restriction explicitly named. |
| 6a | Handling of quantitative variables, scales, and models | Addressed | Methods: Primary MR; AF-Mediated Effect Analysis | State protein units as genetically predicted protein abundance per UKB-PPP pQTL scale if final unit wording is available. |
| 6b | Handling and weighting of genetic variants | Addressed | Methods: Exposure Data; Primary MR; Reverse MR; Mediation | None. |
| 6c | MR estimator, covariates, and two-sample adjustment | Partial | Methods: Primary MR; Outcome GWAS Data | Add covariate-adjustment details from original GWAS/pQTL publications if target journal requires. |
| 6d | Missing data handling | Addressed | Methods: Allele Harmonization; Replication Analyses | None. |
| 6e | Multiple testing correction | Addressed | Methods: Primary MR | None. |
| 7 | Methods to assess MR assumptions | Addressed | Methods: Sensitivity Analyses; Formal Colocalization; Reverse MR | None. |
| 8 | Sensitivity and additional analyses | Addressed | Methods: Replication; Sensitivity; Reverse MR; AF-Mediated Effect; Target Prioritization | None. |
| 9a | Statistical software and packages | Addressed | Methods: Statistical Software and Reproducibility; Supplementary Table S22 | Keep the frozen environment table with the submission package. |
| 9b | Protocol registration | Addressed | Protocol Registration | Current wording: not prospectively registered. |
| 10a | Numbers at each stage and exclusion reasons | Addressed | Results: Exposure Instruments and Outcome Coverage; Tables S1-S3 | Add a figure/flow diagram for final submission. |
| 10b | Phenotypic summary statistics | Partial | Results; Table plan | Most public data are summary-level; add available phenotype summary details from source papers where possible. |
| 10c | Heterogeneity in meta-analytic sources | Partial | Methods/Results; Discussion limitations | Add source-level heterogeneity statements for Nielsen AF and HERMES HF if available. |
| 10d-i | Similarity of exposure and outcome genetic associations in two-sample MR | Partial | Methods; Discussion limitations | State primarily European ancestry and discuss pQTL/outcome ancestry comparability. |
| 10d-ii | Sample overlap between exposure and outcome samples | Partial | Methods: Replication; Discussion limitations | Explicitly state possible UKB overlap for UKB/OpenGWAS and uncertainty/expected overlap for primary outcome datasets. |
| 11a | SNP-exposure and SNP-outcome associations | Addressed | Supplementary tables S1-S3; Results | Ensure final supplement includes beta, SE, effect alleles, and P values. |
| 11b | MR estimates with uncertainty on interpretable scale | Addressed | Results; Tables S3-S6; S16-S21 | None. |
| 11c | Absolute risk translation | Not applicable / optional | NA | Not necessary for this genetic target-prioritization paper unless requested by journal. |
| 11d | Visualization of main results | Addressed | Figures already generated; Tables and Figures section | Insert figure files into final manuscript package. |
| 12a | Assessment of IV assumption validity | Addressed | Results: Sensitivity Analyses and Reverse MR; Discussion | None. |
| 12b | Additional assumption statistics | Addressed | Results: F statistics, proxy-Steiger, colocalization, reverse MR | Optional: report Q statistics for AF->HF mediation step in supplement. |
| 13a | Robustness/sensitivity analyses | Addressed | Results: Sensitivity Analyses and Reverse MR | None. |
| 13b | Additional analyses | Addressed | Replication, colocalization, mediation, target prioritization | None. |
| 13c | Direction of causal relationship | Addressed | Proxy-Steiger and candidate reverse MR | Keep caveats about binary outcome proxy-Steiger and distance-pruned reverse MR. |
| 13d | Comparison with non-MR estimates | Partial / optional | Discussion | Consider adding a short paragraph comparing with biomarker/proteomic literature if target journal expects it. |
| 13e | Additional plots | Partial | Figures already generated; colocalization regional plots planned | Generate final regional colocalization plots if required. |
| 14 | Key results linked to objectives | Addressed | Discussion first paragraph | None. |
| 15 | Limitations and bias direction/magnitude | Addressed | Discussion limitations paragraph | Consider adding a sentence on likely direction of sample-overlap bias if final overlap estimates are available. |
| 16a | Meaning and cautious interpretation | Addressed | Discussion | None. |
| 16b | Biological mechanisms and gene-environment equivalence | Addressed | Introduction and Discussion | Keep careful wording that genetic proxies do not equal short-term pharmacologic intervention. |
| 16c | Clinical/public health relevance | Addressed | Discussion; Target Prioritization; LPA tractability | None. |
| 17 | Generalizability across populations, timing, exposure levels | Addressed | Discussion limitations | None; could expand for non-European ancestry if journal requests. |
| 18 | Funding and funder role | Partial | Funding section; submission statements file | Default no-specific-funding statement drafted; author must replace if any funding applies. |
| 19 | Data and code availability | Addressed | Data Availability; Code Availability; Table 1 | GitHub repository URL and Zenodo DOI added. |
| 20 | Conflicts of interest | Partial | Conflicts of Interest; submission statements file | Default no-competing-interests statement drafted; author must confirm before submission. |

## Overall Assessment

The unified English draft substantially addresses the STROBE-MR structure. Table 1, software-environment freeze, data/code availability draft, protocol-registration wording, default funding/conflict statements, GitHub code repository URL, and Zenodo DOI have now been added. The main remaining gaps are submission-level rather than analytic: target-journal formatting, author confirmation of funding/conflict declarations, final correspondence email, and final reference-manager export.

## Priority Fixes Before Submission

1. Confirm the default no-specific-funding and no-competing-interests statements.
2. Add the final corresponding-author email.
3. Choose the target journal and export references in that journal style.
4. Insert final figures and supplementary tables into the journal-ready package.
5. Confirm whether the target journal requires a local institutional-review exemption statement.
