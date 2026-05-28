# Target Prioritization Scorecard

Scope: FGF5 and LPA, the two exploratory AF-HF shared candidates selected by the rule of one primary outcome FDR < 0.05, the other primary outcome P < 0.05, and consistent risk direction.

Scoring framework: primary MR evidence (0-4), same-direction AF-HF signal (0-1), FinnGen replication (0-2), formal coloc.abf evidence (0-4), exploratory AF-mediated pathway support (0-1), instrument confidence (0-1), clinical tractability (0-2), and predefined evidence penalty.

| Rank | Protein | Priority class | Overall score | Genetic score | Key interpretation |
|---|---|---|---:|---:|---|
| 1 | FGF5 | Tier 1 - 遗传证据高优先级 | 10.0 | 9.5 | 当前遗传证据优先级最高的候选：以AF为主导，主MR、FinnGen AF复制、等位基因方向和FGF5-AF正式共定位一致；HF证据方向一致但仅为名义显著，且未获正式共定位支持。 |
| 2 | LPA | Tier 2 - 中等优先级/转化可开发性突出 | 8.0 | 6.0 | 次级候选：HF主MR和FinnGen AF/HF复制方向一致，且Lp(a)靶向药物开发基础很强；但本研究正式共定位未支持AF/HF共享单一因果变异。 |

## Evidence Notes

### FGF5

- Primary MR: AF reaches FDR significance and HF is nominally significant. AF: OR 1.063 (1.041-1.085), P=7.57e-09, FDR=3.70e-06; HF: OR 1.029 (1.004-1.055), P=0.0205, FDR=0.4417.
- Replication: FinnGen supports the same risk-increasing direction for AF; strict HF exact-variant replication is unavailable or not significant.
- Formal coloc: Strong formal coloc for AF; HF does not show strong coloc. AF PP.H4=0.987, HF PP.H4=0.054, AF PP.H3=0.013, HF PP.H3=0.068.
- Mediation: Exploratory AF-mediated indirect effect: OR 1.014 (1.009-1.02), P=3.87e-08, mediated proportion 49.7%.
- Clinical tractability: No mature FGF5-directed cardiovascular drug-development signal found in the targeted public scan; prior AF-focused MR literature exists.
- Caution: HF support remains nominal and HF formal coloc is not supportive; FGF5-AF has prior MR literature, reducing novelty for an AF-only claim.

### LPA

- Primary MR: HF reaches FDR significance and AF is nominally significant. AF: OR 1.035 (1.01-1.06), P=0.0062, FDR=0.2167; HF: OR 1.08 (1.05-1.111), P=1.10e-07, FDR=2.26e-05.
- Replication: FinnGen supports the same risk-increasing direction for AF and strict HF.
- Formal coloc: No formal coloc support under the single-causal-variant coloc.abf model. AF PP.H4=0.099, HF PP.H4=0.04, AF PP.H3=0.405, HF PP.H3=0.96.
- Mediation: Exploratory AF-mediated indirect effect: OR 1.008 (1.002-1.014), P=0.0068, mediated proportion 10.5%.
- Clinical tractability: High: multiple Lp(a)-lowering RNA therapies targeting apolipoprotein(a)/LPA are in phase 3 cardiovascular outcome trials.
- Caution: Formal coloc is not supportive for AF/HF and the LPA region shows strong distinct-signal/LD concern, especially for HF.

## External Sources Used For Tractability/Precedence

- Open Targets target-prioritisation/tractability documentation: https://platform-docs.opentargets.org/web-interface/target-prioritisation and https://platform-docs.opentargets.org/target/tractability
- FGF5 AF MR prior literature: https://pubmed.ncbi.nlm.nih.gov/39059473/
- Pelacarsen Lp(a)HORIZON: https://clinicaltrials.gov/study/NCT04023552
- Olpasiran OCEAN(a)-Outcomes: https://clinicaltrials.gov/study/NCT05581303
- Lepodisiran ACCLAIM-Lp(a): https://clinicaltrials.gov/study/NCT06292013

CSV scorecard: `D:\Backup\Documents\New project 2\cardio_no_lab_paper_project\results\prioritization\target_priority_scorecard.csv`
