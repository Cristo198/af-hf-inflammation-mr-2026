# Genetically Supported Inflammatory Mediators of the Atrial Fibrillation-Heart Failure Continuum: cis-pQTL Mendelian Randomization, Dual-Outcome Colocalization, and Mediation-Based Target Prioritization

Yangfeng Qin^1* and Runchuan Feng^2*

^1 Langping Township Health Center, Bobai County, Guangxi, China.

^2 Nanning Fourth People's Hospital, Nanning, Guangxi, China.

*Yangfeng Qin and Runchuan Feng contributed equally to this work.

Correspondence: Yangfeng Qin, Langping Township Health Center, Bobai County, Guangxi, China. Email: [Yangfeng Qin email to be inserted].

Manuscript draft: v0.1 unified English version, 2026-05-28

## Abstract

### Background

Atrial fibrillation (AF) and heart failure (HF) frequently coexist and may reinforce each other through shared mechanisms, including inflammation, atrial remodeling, myocardial fibrosis, and neurohormonal activation. Observational biomarker studies have linked inflammatory mediators to AF and HF, but whether these proteins represent causal mediators, disease consequences, or correlated markers remains uncertain.

### Objectives

We aimed to identify genetically supported inflammatory mediators of the AF-HF continuum by integrating cis protein quantitative trait loci (cis-pQTLs), two-sample Mendelian randomization (MR), dual-outcome colocalization, AF-mediated effect analysis, replication, and target prioritization.

### Methods

We used UK Biobank Pharma Proteomics Project (UKB-PPP) Olink Explore inflammatory panels as the primary exposure source. Lead cis-pQTLs were selected for circulating inflammatory proteins at P < 5 x 10^-8 within +/-1 Mb of the encoding gene, with major histocompatibility complex variants excluded. Primary outcomes were AF from Nielsen et al. and HF from the HERMES consortium. Wald ratio MR was used for single-variant protein instruments. Shared candidates were defined by one outcome reaching false discovery rate (FDR) significance, the other reaching nominal significance, and concordant effect direction. Candidate evidence was evaluated using FinnGen R12 replication, UKB/OpenGWAS supplementary replication, formal coloc.abf colocalization with dense UKB-PPP pGWAS data, feasible sensitivity analyses, candidate reverse MR, AF-mediated effect analysis, and a structured target-prioritization scorecard.

### Results

Among 737 Olink inflammatory-panel proteins, 540 had at least one genome-wide significant cis-pQTL, and 529 lead non-MHC cis-pQTL instruments were retained for primary MR. Five proteins reached FDR significance for AF and four reached FDR significance for HF. FGF5 and LPA met the exploratory shared-candidate rule. FGF5 showed AF-dominant evidence, including primary AF MR significance, FinnGen AF replication, consistent allele direction, and strong FGF5-AF colocalization (PP.H4=0.987). FGF5-HF showed nominal MR support but not formal colocalization. LPA showed HF-dominant MR evidence and same-direction FinnGen replication for AF and strict HF, but formal colocalization did not support a single shared causal variant for either AF or HF. In exploratory mediation analysis, AF liability was associated with HF risk (OR 1.265, 95% CI 1.233-1.298; P=3.57 x 10^-71). The AF-mediated indirect effect was stronger for FGF5 (OR 1.0145, 95% CI 1.0093-1.0197; P=3.87 x 10^-8; mediated proportion 49.7%) than for LPA (OR 1.0081, 95% CI 1.0022-1.0140; P=0.0068; mediated proportion 10.5%). Target prioritization ranked FGF5 as a Tier 1 AF-dominant genetic candidate and LPA as a Tier 2 clinically tractable but genetically complex candidate.

### Conclusions

This study provides genetic evidence supporting selected inflammatory proteins in the AF-HF continuum. FGF5 showed the most coherent evidence as an AF-dominant candidate with exploratory AF-mediated HF relevance. LPA showed strong translational appeal but lacked formal colocalization support in the tested AF/HF loci. These findings support a tiered framework that combines cis-pQTL MR, replication, colocalization, mediation analysis, and tractability evidence before advancing inflammatory proteins as cardiovascular targets.

Keywords: atrial fibrillation; heart failure; inflammatory protein; Mendelian randomization; cis-pQTL; colocalization; mediation analysis; target prioritization

## Introduction

Atrial fibrillation (AF) and heart failure (HF) are increasingly encountered as an overlapping clinical syndrome rather than two isolated diagnoses. Contemporary cardiovascular statistics and major clinical guidelines continue to emphasize AF and HF as high-burden conditions that drive hospitalizations, long-term disability, and health-care use across ageing populations [16-18]. Their clinical interaction is particularly consequential: AF can worsen ventricular filling, promote rapid or irregular ventricular rates, and contribute to tachycardia-mediated cardiomyopathy, whereas HF promotes atrial stretch, fibrosis, neurohormonal activation, and electrical remodeling that favor AF initiation and maintenance. Once AF and HF coexist, patients often enter a self-reinforcing cycle of symptoms, decompensation, thromboembolic risk, and therapeutic complexity. This has shifted attention from single-disease models toward an AF-HF continuum shaped by shared upstream biology.

Inflammation is a plausible biological link within this continuum. In AF, inflammatory signaling in atrial cardiomyocytes and immune-cell interactions can influence calcium handling, ion-channel function, oxidative stress, endothelial activation, thrombogenicity, and atrial fibrosis [19]. In HF, innate and adaptive immune pathways contribute to myocardial injury, extracellular matrix remodeling, vascular dysfunction, and adverse ventricular remodeling across the ejection-fraction spectrum [20]. A recent American College of Cardiology scientific statement further highlights inflammation as a cross-cutting determinant of cardiovascular risk, biomarkers, imaging phenotypes, and therapeutic development [21]. Circulating inflammatory proteins therefore provide an attractive interface between systemic immune activation and cardiac structural or electrical disease. However, observational associations between inflammatory biomarkers and AF or HF are vulnerable to confounding by age, obesity, kidney dysfunction, diabetes, coronary artery disease, medication use, and disease severity. They also cannot reliably distinguish causal mediators from downstream consequences of established cardiovascular disease.

Mendelian randomization (MR) offers a complementary strategy for causal inference by using genetic variants as instruments for modifiable exposures, thereby reducing susceptibility to conventional confounding and reverse causation when core assumptions are met [1,2]. The rapid expansion of plasma proteomic genome-wide association studies has made protein quantitative trait loci (pQTLs) especially useful for target-oriented MR [5-8]. Cis-pQTLs, located near the gene encoding the measured protein, are particularly relevant because they are biologically closer to the encoded target and less likely than trans-pQTLs to act through broad upstream networks. Yet MR alone is insufficient for target prioritization: a protein-disease association may arise because the pQTL is in linkage disequilibrium with a nearby disease variant rather than because the same causal variant influences both protein abundance and disease risk. Bayesian colocalization can address this by estimating whether protein and disease association signals share a regional causal variant [3].

Several recent studies have examined inflammatory proteins in relation to AF or HF, and broader proteome-wide analyses have begun to nominate drug targets for these conditions [11-15]. These efforts provide important groundwork, but important gaps remain. First, many inflammation-focused MR studies considered AF or HF separately, rather than modeling them as connected outcomes within a shared pathophysiological continuum. Second, candidate signals have not always been evaluated with dense regional pQTL data and formal dual-outcome colocalization. Third, replication across independent or semi-independent outcome datasets remains variable. Fourth, statistical evidence must still be translated into a clinically interpretable target-prioritization framework that accounts for genetic support, allelic direction, colocalization, reproducibility, druggability, and existing therapeutic context.

We therefore designed a cis-pQTL MR study to identify genetically supported inflammatory mediators of the AF-HF continuum. We integrated UKB-PPP Olink inflammatory panels with large-scale AF and HF genome-wide association studies, followed by harmonized replication, formal regional colocalization for shared candidates, feasible sensitivity analyses, candidate reverse MR, AF-mediated effect analysis, and evidence-based target prioritization. We hypothesized that a subset of circulating inflammatory proteins would show directionally consistent genetic support across AF and HF, and that integrating MR with colocalization, replication, and mediation evidence would distinguish higher-confidence causal candidates from signals more likely to reflect linkage disequilibrium, pleiotropy, or disease-related reverse causation.

## Methods

### Study Design

This was a two-sample summary-level genetic study using cis-pQTL MR to evaluate the potential causal effects of circulating inflammatory proteins on AF and HF. The workflow included exposure selection, cis-pQTL instrument construction, outcome effect extraction, allele harmonization, primary MR, shared-candidate selection, replication, formal regional colocalization, feasible sensitivity analyses, candidate reverse MR, exploratory AF-mediated effect analysis, and target-prioritization scoring. All protein MR estimates were interpreted as the association between genetically predicted higher circulating protein abundance and disease risk.

The primary MR analysis relied on the three core instrumental-variable assumptions: genetic instruments are associated with the exposure (relevance), are not associated with confounders of the exposure-outcome relationship (independence), and affect the outcome only through the exposure of interest (exclusion restriction). We attempted to strengthen these assumptions by using genome-wide significant cis-pQTLs, excluding MHC-region instruments, assessing instrument strength, harmonizing alleles, conducting proxy-Steiger directionality checks, using colocalization as a local LD sensitivity analysis, and performing candidate reverse MR.

### Exposure Data and Genetic Instruments

The primary exposure source was the UKB-PPP Olink Explore 3072 proteomics pQTL resource. We locked the Olink Inflammation and Inflammation_II panels, comprising 737 proteins. Cis regions were defined as +/-1 Mb around the target protein-encoding gene. For the main analysis, European ancestry cis-pQTLs with P < 5 x 10^-8 were eligible, and one lead cis-pQTL was retained per protein. Instruments in the major histocompatibility complex region were excluded. Instrument strength was assessed using F statistics, with F <= 10 considered potentially weak.

In total, 540 proteins had at least one genome-wide significant cis-pQTL. After lead cis-pQTL selection and MHC exclusion, 529 proteins were included in the primary MR analysis.

### Outcome GWAS Data

The primary AF outcome used summary statistics from Nielsen et al., including 60,620 AF cases and 970,216 controls of European ancestry. The primary HF outcome used HERMES consortium summary statistics, including 47,309 HF cases and 930,014 controls of European ancestry. Both primary outcome datasets were downloaded locally and used for direct extraction of instrument-outcome associations. Data-source details, sample sizes, ancestry, phenotype definitions, local files, and access links are summarized in Table 1.

Replication analyses used FinnGen R12 summary statistics for I9_AF, a registry-based atrial fibrillation and flutter endpoint, and I9_HEARTFAIL, a strict registry-based heart failure endpoint. A supplementary UKB/OpenGWAS check used `ukb-b-964` and `ukb-d-HEARTFAIL` VCF.GZ files. The UKB/OpenGWAS analyses were considered supplementary because of low variant coverage and potential sample overlap with UKB-PPP.

### Allele Harmonization

Exposure and outcome data were harmonized by SNP, effect allele, non-effect allele, beta, standard error, and allele frequency when available. If the outcome effect allele was opposite to the protein-increasing allele, the outcome beta was flipped so that all MR estimates corresponded to genetically predicted higher protein abundance. Variants that could not be aligned were excluded. FGF5 and LPA candidate instruments underwent additional manual allele-direction auditing.

### Primary MR Analysis

Because the primary design retained one lead cis-pQTL per protein, the Wald ratio was used for the main protein-to-disease estimates. Binary outcome results were reported as odds ratios (ORs) with 95% confidence intervals (CIs). Multiple testing was controlled using the Benjamini-Hochberg FDR method. Primary statistical significance was defined as FDR < 0.05.

Exploratory AF-HF shared candidates were selected using a prespecified rule: one primary outcome had to reach FDR < 0.05, the other had to reach P < 0.05, and both estimates had to point in the same direction. FGF5 and LPA met this rule and were carried forward as shared-candidate proteins.

### Replication Analyses

Candidate exact-variant replication for FGF5 and LPA was performed in FinnGen R12 PheWeb data for I9_AF (atrial fibrillation and flutter) and I9_HEARTFAIL (heart failure, strict). FinnGen alleles were harmonized to the UKB-PPP protein-increasing allele before Wald ratio estimates were calculated.

Full-panel FinnGen replication was then performed locally for all 529 lead cis-pQTL instruments using downloaded FinnGen R12 AF and HF summary statistics. The FinnGen AF and HF files passed gzip reading checks, and each outcome provided 496 harmonized instruments. UKB/OpenGWAS supplementary replication was performed using locally downloaded VCF.GZ files, but coverage was low: only seven harmonized instruments were available for each UKB/OpenGWAS outcome, and neither FGF5 nor LPA was covered.

### Formal Regional Colocalization

For FGF5 and LPA, dense UKB-PPP European discovery pGWAS summary statistics were used for formal coloc.abf-style colocalization against the primary AF and HF GWAS results. For each candidate protein, a +/-1 Mb region around the lead pQTL was analyzed after rsID mapping and allele harmonization. The coloc.abf priors were p1=1 x 10^-4, p2=1 x 10^-4, and p12=1 x 10^-5. Prior variances were set to 0.15^2 for pQTL effects and 0.2^2 for binary outcome effects. The main colocalization metric was PP.H4, the posterior probability that the protein pQTL and disease GWAS signals share a causal variant. PP.H4 > 0.80 was interpreted as strong colocalization support.

### Sensitivity Analyses

Because the primary MR design used a single lead cis-pQTL per protein, MR-Egger, weighted median, MR-PRESSO, and leave-one-out analyses were not statistically applicable to the main estimates. Feasible sensitivity checks included instrument-strength assessment, palindromic allele checks, allele harmonization review, proxy-Steiger directionality checks, formal colocalization, and candidate reverse MR. The proxy-Steiger analysis compared exposure R2 proxy with outcome R2 proxy on the observed log-odds scale for binary outcomes and was interpreted as supportive rather than definitive directionality evidence.

### Candidate Reverse MR

Candidate reverse MR assessed whether AF or HF genetic liability might influence FGF5 or LPA protein levels. AF and HF instruments were selected from the local primary outcome GWAS files at P < 5 x 10^-8 and distance-pruned at 10 Mb because no LD reference panel was available. Protein-outcome effects were extracted from the genome-wide UKB-PPP pGWAS files for FGF5 and LPA. Fixed-effect IVW was used when multiple disease instruments were available, and Wald ratio was used when only one instrument was available.

To avoid reusing target-protein cis signals as disease instruments, reverse MR was reported for both all distance-pruned disease instruments and a more conservative target-cis-excluded set that removed disease instruments within +/-1 Mb of the target protein cis locus. These analyses were treated as candidate-level screening checks rather than definitive bidirectional MR.

### AF-Mediated Effect Analysis

To explore whether AF may mediate part of the protein-HF association, we used a summary-level product-of-coefficients approach for the pathway genetically predicted protein abundance -> AF -> HF. First, AF instruments were selected from the local AF GWAS at P < 5 x 10^-8 and distance-pruned at 10 Mb. The corresponding HF effects were extracted from HERMES and harmonized to the AF effect allele. The AF -> HF association was estimated using fixed-effect IVW.

For each protein with both protein -> AF and protein -> HF primary MR estimates, the indirect effect was calculated as beta(protein -> AF) x beta(AF -> HF). The standard error was approximated using the delta method: SE_indirect = sqrt(beta_AF_to_HF^2 x SE_protein_to_AF^2 + beta_protein_to_AF^2 x SE_AF_to_HF^2). The direct effect was approximated as the total protein -> HF effect minus the indirect effect, and the mediated proportion was calculated as indirect effect divided by total effect. Because AF and HF are binary outcomes and the AF instrument set was distance-pruned rather than LD-reference clumped, mediation estimates were considered exploratory pathway evidence.

### Target Prioritization

FGF5 and LPA were prioritized using a structured scorecard incorporating primary MR evidence, AF/HF direction concordance, FinnGen replication, formal colocalization, exploratory AF-mediated support, instrument confidence, allele-direction audit, reverse MR findings, clinical tractability, and public drug-development context. Open Targets tractability information and ClinicalTrials.gov records for Lp(a)-lowering therapies were used as translational context. The scorecard was intended as an internal evidence-ranking tool, not as a direct judgment of clinical efficacy or safety.

### Statistical Software and Reproducibility

Data processing, harmonization, MR, colocalization, sensitivity analyses, reverse MR, mediation analysis, target prioritization, and table generation were performed using Python scripts stored in the project `scripts/` directory. The frozen local analysis environment used Python 3.12.4 on Windows 10.0.19045. Core analysis scripts relied on Python standard-library modules for file parsing, statistics, compression, and table generation; `synapseclient` 4.12.0 was used only to assist controlled UKB-PPP pGWAS downloads. No R scripts were used. The software environment is summarized in Supplementary Table S22. The study protocol was not prospectively registered.

### Ethics

This study used publicly available or access-controlled summary-level GWAS and pQTL data and did not involve new individual-level data collection. Ethical approval and informed consent were obtained in the original studies. Additional institutional review is expected to be unnecessary for the present summary-level analysis, subject to local institutional requirements.

## Results

### Exposure Instruments and Outcome Coverage

The primary exposure resource comprised 737 UKB-PPP Olink inflammatory-panel proteins, including 368 proteins from Inflammation and 369 from Inflammation_II. Using the UKB-PPP European pQTL results, P < 5 x 10^-8, and a +/-1 Mb cis window, 4,229 significant pQTL associations were screened, of which 545 were in cis. A total of 540 proteins retained at least one valid cis-pQTL. After lead cis-pQTL selection and MHC exclusion, 529 proteins entered the primary MR analysis.

All lead instruments had F statistics >10. The median F statistic was 788.5 and the minimum was 45.8, suggesting low weak-instrument risk. The F statistics for FGF5 and LPA were 7,582.3 and 5,377.8, respectively.

For the AF outcome, 492 unique instruments were locally matched and 489 remained after allele harmonization. For the HF outcome, 418 unique instruments were matched and 410 remained after harmonization.

### Inflammatory Proteins and AF Risk

In the primary Wald ratio MR analysis, 37 inflammatory-panel proteins were nominally associated with AF risk and five reached FDR significance. The FDR-significant proteins were FGF5, TNFSF12, SPON1, PKLR, and NFATC1. Genetically predicted higher FGF5, SPON1, and NFATC1 levels were associated with higher AF risk, whereas genetically predicted higher TNFSF12 and PKLR levels were associated with lower AF risk. FGF5 was one of the strongest AF signals (OR 1.063, 95% CI 1.041-1.085; P=7.57 x 10^-9; FDR=3.70 x 10^-6). TNFSF12 showed an inverse association (OR 0.896, 95% CI 0.862-0.932; P=5.33 x 10^-8; FDR=1.30 x 10^-5).

### Inflammatory Proteins and HF Risk

For HF, 41 inflammatory-panel proteins were nominally associated with risk and four reached FDR significance. These were CELSR2, LPA, ABO, and APOA2. Genetically predicted higher LPA and ABO levels were associated with higher HF risk, whereas genetically predicted higher CELSR2 and APOA2 levels were associated with lower HF risk. CELSR2 showed the strongest HF association (OR 0.883, 95% CI 0.847-0.920; P=3.27 x 10^-9; FDR=1.34 x 10^-6). LPA was associated with higher HF risk (OR 1.080, 95% CI 1.050-1.111; P=1.10 x 10^-7; FDR=2.26 x 10^-5).

### Shared AF-HF Candidate Proteins

Six proteins overlapped between the nominal AF and HF association sets: CXCL17, FGF5, HEG1, HRG, LPA, and ZBP1. Applying the exploratory shared-candidate rule selected FGF5 and LPA. FGF5 reached FDR significance for AF and nominal significance for HF with concordant risk direction. LPA reached FDR significance for HF and nominal significance for AF with concordant risk direction. Allele-direction auditing showed that the FGF5 C allele and the LPA T allele increased the corresponding protein levels and pointed toward higher risk for both outcomes. Venn and UpSet summaries showed limited overlap between AF and HF nominal signals and no direct overlap at the dual-FDR level, indicating that FGF5 and LPA represent AF-dominant and HF-dominant shared-candidate patterns, respectively.

### Replication and Colocalization

Candidate exact-variant FinnGen replication showed that FGF5 was directionally consistent and significant for AF (OR 1.130, 95% CI 1.100-1.160; P=2.68 x 10^-19), whereas the exact strict-HF variant effect was unavailable or not significant. LPA was directionally consistent and significant for FinnGen AF (OR 1.065, 95% CI 1.035-1.096; P=1.49 x 10^-5) and directionally consistent with nominal significance for FinnGen strict HF (OR 1.030, 95% CI 1.003-1.059; P=0.0318).

Full-panel FinnGen R12 replication was completed for all 529 lead cis-pQTL instruments. Each FinnGen outcome had 496 harmonized instruments. In FinnGen AF, 467 proteins with primary AF estimates were matched; 50 showed same-direction nominal replication and 14 reached same-direction FDR replication. In FinnGen HF, 399 proteins with primary HF estimates were matched; 24 showed same-direction nominal replication and one reached same-direction FDR replication. Among primary FDR signals, FGF5 and SPON1 replicated at FDR level for AF, while LPA replicated nominally for HF. Among FGF5/LPA shared candidates, FGF5-AF and LPA-AF reached same-direction FinnGen FDR replication, LPA-HF reached same-direction nominal replication, and FGF5-HF was directionally consistent but not nominally significant.

UKB/OpenGWAS supplementary replication was completed using `ukb-b-964` and `ukb-d-HEARTFAIL` VCF.GZ files. However, coverage was very low: each outcome harmonized only seven instruments, and neither FGF5, LPA, nor any primary FDR signal was covered. No same-direction nominal or FDR replication signal was observed in this low-coverage supplementary analysis.

Formal coloc.abf colocalization using dense UKB-PPP regional pQTL summary statistics showed strong support for FGF5-AF colocalization. In the chr4:79,261,400-81,261,400 region (GRCh38), 9,454 harmonized overlapping SNPs were available for FGF5-AF, with PP.H4=0.987 and PP.H3=0.013; the lead shared SNP was rs12509595. FGF5-HF had 4,724 overlapping SNPs but did not show strong colocalization (PP.H4=0.054, PP.H3=0.068). LPA-AF had 10,607 overlapping SNPs (PP.H4=0.099, PP.H3=0.405), and LPA-HF had 6,414 overlapping SNPs (PP.H4=0.040, PP.H3=0.960), suggesting distinct-signal or LD complexity rather than a shared single causal variant under the current model.

FinnGen precomputed pQTL-disease colocalization records provided additional supportive evidence for FGF5-AF, including CLPP=0.204 and CLPA=0.557 for a FinnGen Olink pQTL record and CLPP=0.249 and CLPA=0.249 for a UK Biobank PPP Olink 3k record. No comparable FinnGen pQTL-disease colocalization record was detected for LPA in AF or strict HF.

### Sensitivity Analyses and Reverse MR

Because the primary MR design used one lead cis-pQTL per protein, MR-Egger, weighted median, MR-PRESSO, and leave-one-out analyses were not applicable to the primary estimates. Among 899 harmonized primary MR records, no instrument had F <= 10. A total of 141 records involved palindromic alleles, but the FGF5 and LPA candidate instruments were not palindromic. Proxy-Steiger directionality checks supported the protein-to-disease direction for 896 of 899 primary MR records and for all four FGF5/LPA candidate rows. This analysis used observed log-odds scale approximations for binary outcomes and was interpreted as supportive rather than definitive.

Candidate reverse MR used genome-wide significant AF and HF variants from the local primary outcome GWAS files and 10 Mb distance pruning. AF yielded 75 distance-pruned instruments and HF yielded 11. After harmonization with genome-wide UKB-PPP pGWAS data for FGF5 and LPA, AF-related reverse MR included 74 instruments and HF-related reverse MR included nine. In the full distance-pruned set, AF genetic liability was positively associated with FGF5 protein levels (beta=0.081, SE=0.011; P=1.40 x 10^-13), and HF genetic liability was positively associated with LPA protein levels (beta=0.925, SE=0.042; P underflowed to 0). After excluding disease instruments within +/-1 Mb of the target protein cis region, the HF -> LPA signal disappeared (beta=-0.0066, SE=0.0437; P=0.879), suggesting that the unfiltered result was driven by local LPA-region signal reuse. The AF -> FGF5 signal remained nominal after target-cis exclusion (beta=0.024, SE=0.011; P=0.0278), but heterogeneity remained high and LD reference-panel clumping was not performed. Therefore, reverse MR findings were treated as exploratory evidence of possible bidirectionality or pleiotropy rather than definitive reverse causation.

### AF-Mediated Effects

To explore whether AF may mediate part of the inflammatory protein-HF association, we first estimated the genetic association between AF liability and HF risk. AF instruments were selected at P < 5 x 10^-8 from the local AF GWAS and distance-pruned at 10 Mb. Of 75 pruned AF instruments, 73 were harmonized with the HERMES HF GWAS. Fixed-effect IVW supported a positive AF -> HF association (OR 1.265, 95% CI 1.233-1.298; P=3.57 x 10^-71).

For FGF5, the AF-mediated indirect effect was OR 1.0145 (95% CI 1.0093-1.0197; P=3.87 x 10^-8), corresponding to an exploratory mediated proportion of 49.7%. For LPA, the AF-mediated indirect effect was OR 1.0081 (95% CI 1.0022-1.0140; P=0.0068), corresponding to an exploratory mediated proportion of 10.5%. These estimates suggest that AF-mediated pathways may contribute to the observed FGF5/LPA-HF associations, particularly for FGF5. Because the analysis used summary-level binary outcome estimates and distance-pruned AF instruments without a formal LD reference panel, the mediation results were treated as exploratory pathway evidence rather than definitive mechanistic proof.

### Target Prioritization

Target prioritization incorporated primary MR evidence, AF/HF direction concordance, FinnGen replication, formal colocalization, exploratory AF-mediated support, instrument confidence, allele-direction audit, reverse MR findings, and clinical tractability. FGF5 received an updated overall score of 10.0 and was ranked as a Tier 1 high-priority genetic candidate. Its strengths were AF FDR-significant MR, nominal and directionally concordant HF MR, FinnGen AF replication, consistent allele direction, strong FGF5-AF colocalization, and exploratory AF-mediated pathway support. Its limitations were incomplete HF replication, lack of FGF5-HF colocalization, and prior FGF5-AF MR literature. Thus, FGF5 is best described as an AF-dominant, HF-directionally concordant genetic candidate with exploratory AF-mediated HF relevance rather than a fully validated dual-outcome AF-HF target.

LPA received an updated overall score of 8.0 and was ranked as a Tier 2 candidate with high translational tractability but cautious genetic localization. Its strengths were HF FDR-significant MR, nominal and directionally concordant AF MR, same-direction FinnGen AF and strict-HF replication, a smaller but nominally significant AF-mediated indirect effect, and multiple ongoing phase 3 Lp(a)-lowering cardiovascular outcome trials. Its main limitation was lack of formal colocalization support for LPA-AF or LPA-HF, especially the high PP.H3 for LPA-HF, suggesting complex LD or distinct causal signals. LPA should therefore be considered a clinically attractive but genetically complex secondary candidate.

## Discussion

In this cis-pQTL MR study of circulating inflammatory proteins, we integrated AF and HF genetic outcome data to prioritize inflammatory mediators relevant to the AF-HF continuum. The main analysis screened 529 lead cis-pQTL instruments from UKB-PPP Olink inflammatory panels and identified distinct AF- and HF-associated protein signals after false-discovery-rate correction. Applying an exploratory shared-candidate rule that required one outcome to reach FDR significance, the second outcome to reach nominal significance, and both effects to point in the same direction nominated FGF5 and LPA for deeper evaluation. Subsequent evidence integration separated these two candidates into different interpretive classes: FGF5 showed the most coherent genetic evidence, driven by AF, FinnGen replication, consistent allelic direction, proxy-Steiger support, strong FGF5-AF colocalization, and exploratory AF-mediated support; LPA showed strong MR and replication features, especially for HF and AF directionality, and a smaller AF-mediated component, but formal colocalization did not support a single shared causal variant in the current model.

These findings are consistent with the concept that AF and HF share upstream biological pathways but also illustrate why a continuum framework should not be reduced to simple overlap of significant MR results. Inflammation can influence atrial electrophysiology, thrombogenic signaling, fibrotic remodeling, endothelial function, and myocardial repair, all of which may contribute to AF initiation, AF persistence, HF onset, or HF progression [19-21]. However, inflammatory proteins measured in plasma can represent multiple biological states, including immune activation, vascular injury, myocardial stress, metabolic disease, and downstream responses to established cardiovascular pathology. By combining cis-pQTL MR with replication and colocalization, the present study attempts to move beyond biomarker association and rank signals by the strength and localization of genetic support.

FGF5 emerged as the current high-priority genetic candidate in this analysis. The AF association was robust in the primary MR analysis and replicated in FinnGen with the same risk-increasing direction. Formal colocalization using dense UKB-PPP regional pQTL data provided strong support that the FGF5 pQTL and AF GWAS signals share a regional causal variant. This is important because prior inflammation-focused AF MR studies have also highlighted FGF5, meaning that an AF-only claim would have limited novelty [11,22]. The contribution of the present work is therefore not merely to rediscover FGF5 for AF, but to place it within a dual-outcome AF-HF framework and to show that its strongest evidence remains AF-dominant. The HF signal was directionally consistent and nominally significant in the primary analysis, but FinnGen HF support was weaker and FGF5-HF colocalization was not convincing. Thus, FGF5 should be interpreted as an AF-prioritized inflammatory mediator with possible HF relevance rather than as a fully validated shared AF-HF target.

LPA had a different evidence profile. Genetically predicted higher LPA protein levels were associated with higher HF risk in the primary analysis, and the AF signal was directionally concordant. FinnGen supported the same direction for both AF and strict HF, with AF reaching FDR-level replication and HF reaching nominal replication. This makes LPA translationally attractive, especially because Lp(a)-lowering therapies targeting apolipoprotein(a)/LPA are already being tested in phase 3 cardiovascular outcome trials [23-25]. Nevertheless, formal colocalization was not supportive for either AF or HF, and the HF analysis showed a high posterior probability for distinct signals. The LPA locus is genetically complex and strongly linked to lipid, atherosclerotic, thrombotic, and valvular biology. Therefore, the LPA result may reflect a biologically meaningful cardiovascular risk axis, but the present data do not establish a single shared inflammatory-protein causal variant for AF or HF at this locus. In the target-prioritization framework, this justifies treating LPA as a clinically tractable but genetically cautious secondary candidate.

The replication and sensitivity analyses further refine these interpretations. FinnGen full-panel replication provided useful external support for several AF signals and confirmed same-direction evidence for the two shared candidates. In contrast, the UKB/OpenGWAS supplementary replication had very low variant coverage and did not include FGF5 or LPA, limiting its value for candidate confirmation. Because UKB-PPP exposures and UKB-derived outcomes may overlap, these UKB results are best viewed as a low-coverage technical check rather than independent validation. Instrument quality checks did not identify weak instruments among the harmonized primary MR rows, and the candidate instruments were not palindromic. Proxy-Steiger analyses supported the protein-to-disease direction for nearly all primary MR rows and all FGF5/LPA candidate rows, although this approximation used observed log-odds scale information for binary outcomes. Candidate reverse MR was also informative but exploratory. The AF-to-FGF5 signal persisted at nominal significance after target-cis exclusion, suggesting possible bidirectionality or residual pleiotropy that warrants caution. Conversely, the apparent HF-to-LPA reverse signal disappeared after excluding the target cis region, indicating that the unfiltered result was likely driven by local LPA-region signal reuse.

The colocalization results are central to the biological interpretation of this study. A significant cis-pQTL MR estimate can arise from a causal protein effect, horizontal pleiotropy, or linkage disequilibrium between separate causal variants. The strong FGF5-AF PP.H4 result strengthens the case for a localized shared signal, whereas the low PP.H4 results for FGF5-HF, LPA-AF, and LPA-HF argue against over-interpreting MR significance as proof of shared causality. At the same time, lack of colocalization under a single-causal-variant coloc.abf model does not exclude biology. It may reflect multiple causal variants, incomplete variant coverage, ancestry differences, assay-specific protein measurement, or complex regional LD. Fine-mapping approaches that allow multiple causal signals, such as conditional colocalization or SuSiE-based colocalization, would be appropriate next steps, particularly for the LPA region.

This work has several implications. Methodologically, it supports a tiered target-prioritization strategy for proteomic MR studies: primary MR can nominate signals, but replication, allele auditing, directionality checks, colocalization, mediation analysis, reverse MR, and druggability context are needed before target claims are made. Clinically, the findings support inflammatory and vascular-inflammatory pathways as plausible contributors to the AF-HF continuum, but they also show that the evidence may be asymmetric across the two diseases. FGF5 may be most relevant to AF biology, with a substantial exploratory AF-mediated component for HF, whereas LPA may reflect a broader cardiovascular risk pathway with clearer translational infrastructure but less localized colocalization support in this analysis. These distinctions matter because target prioritization for AF-HF prevention should account not only for statistical association, but also for disease specificity, mechanistic localization, pathway evidence, and feasibility of therapeutic modulation.

Several limitations should be acknowledged. First, most exposure and outcome datasets were derived primarily from individuals of European ancestry, which limits generalizability to other populations. Second, the main MR design used one lead cis-pQTL per protein, making robust multi-instrument sensitivity methods such as MR-Egger, weighted median, MR-PRESSO, and leave-one-out unsuitable for the primary estimates. Third, plasma protein abundance may not fully capture protein activity in atrial myocardium, ventricular myocardium, vascular tissue, or immune-cell microenvironments. Fourth, cis-pQTLs are more interpretable than trans-pQTLs for target MR, but they can still act through nearby genes, splice effects, protein isoforms, assay-binding effects, or LD with other functional variants. Fifth, formal colocalization depends on regional variant density, allele harmonization, prior settings, and the assumption structure of the model. Sixth, summary-level datasets may differ in phenotype definitions, case ascertainment, and covariate adjustment, and some replication resources may have sample overlap with the exposure source. Seventh, genetic proxies represent lifelong differences in protein abundance and cannot be assumed to reproduce the timing, magnitude, tissue distribution, or safety profile of pharmacological intervention.

The mediation analysis adds pathway-level support but should remain exploratory. In the product-of-coefficients framework, AF mediated a larger proportion of the FGF5-HF association than of the LPA-HF association, which is consistent with an AF-dominant interpretation for FGF5 and a broader cardiovascular-risk interpretation for LPA. However, mediation estimates were derived from summary-level binary outcome data and a distance-pruned AF instrument set without a formal LD reference panel, so they should not be read as definitive proof of a mechanistic chain. Future work should incorporate multi-ancestry pQTL resources, cardiac tissue expression and single-cell data, fine-mapped colocalization, phenome-wide scans for safety signals, multivariable MR, and experimental studies in atrial and myocardial systems. For now, the evidence supports FGF5 as the highest-priority AF-dominant genetic candidate in this project and LPA as a translationally attractive but genetically complex secondary candidate.

## Tables and Figures

Figure 1. Study workflow: pQTL screening, cis instrument construction, AF/HF MR, replication, colocalization, mediation analysis, and target prioritization.

Figure 2. Primary AF and HF MR volcano plots.

Figure 3. Venn or UpSet plot of AF and HF significant inflammatory proteins.

Figure 4. Forest plot of primary and replication effects for FGF5 and LPA.

Figure 5. Regional colocalization plots for FGF5 and LPA candidate loci.

Figure 6. Target-prioritization heatmap or scorecard.

Table 1. Data sources, sample sizes, ancestry, phenotype definitions, local files, and access links.

Table 2. Lead cis-pQTL instruments.

Table 3. Primary MR results.

Table 4. Replication and sensitivity analyses.

Table 5. Formal colocalization and AF-mediated effect estimates.

Table 6. Target-prioritization scorecard.

## Data Availability

All analyses were conducted using summary-level genetic association data. The primary plasma protein pQTL data were obtained from the UK Biobank Pharma Proteomics Project (UKB-PPP) Olink Explore 3072 resource, available through the UKB-PPP Synapse project (https://www.synapse.org/Synapse:syn51364943) and the Registry of Open Data on AWS (https://registry.opendata.aws/ukbppp/), subject to the applicable access conditions. The atrial fibrillation GWAS summary statistics from Nielsen et al. are available from the University of Michigan Center for Statistical Genetics (http://csg.sph.umich.edu/willer/public/afib2018/nielsen-thorolfsdottir-willer-NG2018-AFib-gwas-summary-statistics.tbl.gz) and the associated repository record (https://doi.org/10.18710/VC2PSH). The heart failure GWAS summary statistics from the HERMES consortium are available from the HERMES/Broad download route (https://personal.broadinstitute.org/ryank/HERMES_Jan2019_HeartFailure_summary_data.txt.zip) and are described by Shah et al. (https://doi.org/10.1038/s41467-019-13690-5). FinnGen R12 replication endpoints are available from FinnGen and Risteys (https://r12.finngen.fi/; https://r12.risteys.finngen.fi/endpoints/I9_AF; https://r12.risteys.finngen.fi/endpoints/I9_HEARTFAIL). UKB/OpenGWAS supplementary replication datasets are available through IEU OpenGWAS (https://gwas.mrcieu.ac.uk/datasets/ukb-b-964/; https://gwas.mrcieu.ac.uk/datasets/ukb-d-HEARTFAIL/).

Derived non-restricted summary tables, figure files, analysis scripts, and file manifests are available at https://github.com/Cristo198/af-hf-inflammation-mr-2026 and archived on Zenodo at https://doi.org/10.5281/zenodo.20433665. Restricted or third-party raw summary-statistics files will not be redistributed; users should obtain these files from the original data providers under the relevant terms of use. No individual-level data were generated or analyzed in this study.

## Code Availability

All custom analysis scripts used for exposure selection, cis-pQTL instrument construction, outcome extraction, allele harmonization, Wald ratio MR, replication, formal coloc.abf-style colocalization, sensitivity checks, reverse MR, mediation analysis, target prioritization, and table generation are available at https://github.com/Cristo198/af-hf-inflammation-mr-2026 and archived on Zenodo at https://doi.org/10.5281/zenodo.20433665. The public repository includes the frozen software environment table, derived non-restricted output tables, and a manifest of required external input files.

## Ethics Statement

The analysis used summary-level genetic data from studies with existing ethical approval and participant consent. No new individual-level data were collected for this study. Additional institutional review should be confirmed according to local policy.

## Author Contributions

Yangfeng Qin and Runchuan Feng jointly conceived the study and designed the analysis plan. Yangfeng Qin performed data processing, statistical analyses, evidence synthesis, and manuscript drafting. Runchuan Feng contributed clinical interpretation, cardiovascular framing, manuscript revision, and correspondence. Both authors approved the final manuscript.

## Funding

This research received no specific grant from any funding agency in the public, commercial, or not-for-profit sectors.

## Conflicts of Interest

The authors declare no competing interests.

## Protocol Registration

This study was not prospectively registered. The analytic workflow, prespecified candidate-selection rule, and reproducibility scripts are provided with the manuscript and supplementary materials.

## References

Numbered citations correspond to the working reference seed library in `references_seed.md`. The current draft uses a conservative numbered reference style suitable for many biomedical journals. Before submission, all references should be exported from Zotero, EndNote, or another reference manager in the selected target journal style and checked for author names, title, journal, year, volume, pages, DOI, PMID, and data-access links.
