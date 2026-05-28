# Results段落初稿

说明：以下只更新Results，不包含Introduction或Discussion。正式主结局共定位仍等待FGF5/LPA区域级密集UKB-PPP pQTL summary statistics。

## 3.1 数据源、工具变量和质量控制

本研究从UKB-PPP Olink Explore炎症相关面板中纳入737种蛋白。基于EUR人群pQTL结果、P < 5 × 10^-8阈值和目标基因上下游1 Mb cis窗口，共筛得545条cis-pQTL关联。按每种蛋白保留一个lead cis-pQTL并排除MHC区域后，529种蛋白进入主MR分析。所有lead工具变量F统计量均大于10，F统计量中位数为788.5，最小值为45.8，提示弱工具变量偏倚风险较低。FGF5和LPA的lead工具变量F统计量分别为7582.3和5377.8。由于主分析采用单lead cis-pQTL设计，MR-Egger、weighted median、MR-PRESSO和leave-one-out不适用于主分析估计；本阶段主要通过cis限制、MHC排除、强工具变量筛选、等位基因协调、复制验证和后续共定位控制潜在多效性风险。

心房颤动主结局来自Nielsen等2018年GWAS，心力衰竭主结局来自HERMES GWAS。本地结局提取中，AF结局成功提取492个唯一工具变量，等位基因协调后489个进入MR；HF结局成功提取418个唯一工具变量，等位基因协调后410个进入MR。

## 3.2 循环炎症蛋白与心房颤动风险

在AF主结局中，37种蛋白达到名义显著，5种蛋白达到FDR校正显著。FDR显著蛋白包括：FGF5（OR=1.063, 95%CI 1.041-1.085, P=7.57e-09, FDR=3.70e-06）、TNFSF12（OR=0.896, 95%CI 0.862-0.932, P=5.33e-08, FDR=1.30e-05）、SPON1（OR=1.114, 95%CI 1.065-1.164, P=1.96e-06, FDR=3.19e-04）、PKLR（OR=0.862, 95%CI 0.8-0.928, P=8.84e-05, FDR=0.0108）、NFATC1（OR=1.277, 95%CI 1.119-1.458, P=2.78e-04, FDR=0.0272）。其中，FGF5为最强信号之一，遗传预测FGF5水平升高与AF风险升高相关（OR=1.063, 95%CI 1.041-1.085, P=7.57e-09, FDR=3.70e-06）。

## 3.3 循环炎症蛋白与心力衰竭风险

在HF主结局中，41种蛋白达到名义显著，4种蛋白达到FDR校正显著。FDR显著蛋白包括：CELSR2（OR=0.883, 95%CI 0.847-0.92, P=3.27e-09, FDR=1.34e-06）、LPA（OR=1.08, 95%CI 1.05-1.111, P=1.10e-07, FDR=2.26e-05）、ABO（OR=1.034, 95%CI 1.02-1.049, P=3.18e-06, FDR=4.35e-04）、APOA2（OR=0.643, 95%CI 0.507-0.814, P=2.46e-04, FDR=0.0252）。其中，LPA与HF风险升高显著相关（OR=1.08, 95%CI 1.05-1.111, P=1.10e-07, FDR=2.26e-05）。

## 3.4 AF-HF共享炎症蛋白筛选

AF和HF名义显著信号共有6种蛋白重叠，分别为CXCL17、FGF5、HEG1、HRG、LPA、ZBP1。进一步采用“一个结局FDR < 0.05、另一个结局P < 0.05且方向一致”的探索性共享候选标准后，FGF5和LPA被筛选为AF-HF连续体候选蛋白。FGF5在AF中达到FDR显著（OR=1.063, 95%CI 1.041-1.085, P=7.57e-09, FDR=3.70e-06），在HF中达到名义显著且方向一致（OR=1.029, 95%CI 1.004-1.055, P=0.0205, FDR=0.4417）。LPA在HF中达到FDR显著（OR=1.08, 95%CI 1.05-1.111, P=1.10e-07, FDR=2.26e-05），在AF中达到名义显著且方向一致（OR=1.035, 95%CI 1.01-1.06, P=0.0062, FDR=0.2167）。等位基因方向复核显示，FGF5的C等位基因和LPA的T等位基因分别为蛋白升高等位基因，两个候选蛋白在AF和HF中的主分析方向均指向风险增加。

Venn和UpSet图显示，AF和HF名义显著集合存在有限重叠，而FDR层面的双结局直接重叠为0；FGF5和LPA分别代表AF主导和HF主导的共享候选模式。

## 3.5 FinnGen复制和预计算共定位线索

候选复制使用FinnGen R12 PheWeb公开变异接口。FGF5在FinnGen AF中方向一致并显著（OR=1.13, 95%CI 1.1-1.16, P=2.68e-19），但FinnGen strict HF在该精确变异处无有效beta。LPA在FinnGen AF中方向一致并显著（OR=1.065, 95%CI 1.035-1.096, P=1.49e-05），在FinnGen strict HF中方向一致并达到名义显著（OR=1.03, 95%CI 1.003-1.059, P=0.0318）。

FinnGen R12预计算pQTL-疾病共定位记录为FGF5-AF提供支持性线索：FinnGen Olink记录CLPP=0.204、CLPA=0.557，UK Biobank PPP Olink 3k记录CLPP=0.249、CLPA=0.249。FGF5-HF及LPA-AF/HF暂未检索到对应的FinnGen pQTL-疾病共定位记录。需要强调的是，这些CLPP/CLPA结果属于FinnGen fine-mapping框架下的预计算支持证据，并不等同于本研究主结局Nielsen AF和HERMES HF的正式coloc.abf PP.H4结果。正式主结局共定位仍需获取FGF5和LPA区域级密集pQTL summary statistics后完成。

## 3.6 图表和补充表

本阶段已生成AF/HF主MR火山图、FDR显著结果森林图、FGF5/LPA主分析与复制森林图、AF-HF名义重叠Venn图和MR显著性UpSet图。补充表S1-S10整理了工具变量QC、结局协调、全部MR结果、FDR显著结果、AF/HF重叠蛋白、FGF5/LPA证据矩阵、等位基因方向复核、FinnGen复制和预计算共定位记录。

## 暂不填写的Results小节

中介分析和靶点优先级排序尚未完成，相关Results段落暂不填写。Introduction和Discussion按用户后续专业提示词另行处理。