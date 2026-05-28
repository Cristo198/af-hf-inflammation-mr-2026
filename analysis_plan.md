# 分析计划

## 1. 文献和重复性检查

先检索近3年是否已有完全相同题目。重点检查：

1. 炎症蛋白 -> AF的MR。
2. 炎症蛋白 -> HF的MR。
3. 同时分析AF和HF共同炎症靶点的研究。
4. 是否包含最新HERMES亚型GWAS或最新pQTL。

如果已有类似文章，差异化方向：

- 加入HF亚型。
- 使用cis-pQTL和共定位作为主证据。
- 做AF-HF共享靶点，而不是单疾病列表。
- 加入药物靶点优先级。

## 2. 工具变量筛选

暴露为每个炎症蛋白。

主暴露数据源使用UKB-PPP Olink Explore炎症相关蛋白pQTL。主文主分析优先纳入inflammation和inflammation II面板蛋白；扩展分析纳入经GO、Reactome、KEGG、UniProt或Open Targets注释为免疫-炎症相关的血浆蛋白。deCODE SomaScan pQTL用于跨平台验证，Zhao/SCALLOP 91炎症蛋白pQTL用于经典炎症面板敏感性分析。

筛选标准：

1. genome-wide significance: P < 5e-8。
2. 若工具变量太少，可在探索性分析中使用P < 1e-5，但主文需明确为探索性。
3. LD clumping: r2 < 0.001，窗口10,000 kb。
4. F statistic > 10。
5. 主分析仅用cis-pQTL，cis窗口定义为编码基因上下游1Mb。
6. trans-pQTL不进入主分析，仅作为机制网络或补充分析。
7. 对protein-altering variants进行标记；必要时做剔除敏感性分析。
8. MHC区域(chr6:25.5-34.0Mb)做剔除敏感性分析。

## 3. 数据协调

1. 对齐效应等位基因。
2. 移除palindromic ambiguous SNP。
3. 统一beta、SE、effect allele、other allele、EAF、P值。
4. 记录被剔除SNP和原因。

## 4. MR主分析

根据工具变量数量选择：

- 1个SNP：Wald ratio。
- 2个SNP：IVW fixed effect为主，敏感性有限。
- >=3个SNP：IVW random effects为主。

结果报告：

- OR。
- 95%CI。
- P值。
- FDR校正P值。
- SNP数量。

## 5. 敏感性分析

1. MR-Egger。
2. Weighted median。
3. Weighted mode。
4. MR-PRESSO。
5. Cochran Q异质性。
6. MR-Egger intercept水平多效性。
7. Leave-one-out。
8. Steiger方向性检验。

## 6. 反向MR

对显著结果做：

AF -> 炎症蛋白。

HF -> 炎症蛋白。

目的：判断结果是否可能由疾病倾向反向影响炎症蛋白水平。

## 7. 复制分析

使用独立或半独立结局数据：

1. FinnGen AF/HF。
2. UK Biobank AF/HF。
3. Biobank Japan作为跨ancestry探索性验证。

复制标准：

- 方向一致。
- P < 0.05或FDR后仍显著。
- 无明显异质性/多效性。

## 8. 共定位或SMR/HEIDI

对重点蛋白做局部共定位：

1. 提取pQTL位点上下游500kb或1Mb。
2. 对应疾病GWAS同一区域。
3. 使用`coloc`估计PP.H4。
4. PP.H4 > 0.70可视为支持共享因果变异；>0.80更强。

如缺少完整区域summary statistics，可用SMR/HEIDI作为替代。

## 9. 靶点优先级

每个候选蛋白按以下维度评分：

1. AF显著性。
2. HF显著性。
3. 复制结果。
4. 共定位证据。
5. cis-pQTL证据。
6. 是否已有药物或抗体。
7. 是否有心血管安全性信号。
8. 生物学合理性。

## 10. 图表

建议图表：

1. 研究流程图。
2. AF和HF MR结果火山图。
3. AF-HF共同阳性蛋白Venn图或UpSet图。
4. 重点蛋白森林图。
5. 共定位区域图。
6. 蛋白互作网络图。
7. 靶点优先级热图。
