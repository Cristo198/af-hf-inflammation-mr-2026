# 数据源

## 结局数据

### 心房颤动

优先数据：

- IEU OpenGWAS: `ebi-a-GCST006414`
- 表型：Atrial fibrillation
- 样本量：约1,030,836；病例60,620，对照970,216
- 用途：主分析AF结局

可选复制数据：

- FinnGen: atrial fibrillation and flutter
- UK Biobank ICD10 I48相关表型
- Biobank Japan：东亚人群心律失常/AF相关表型

### 心力衰竭

优先数据：

- IEU OpenGWAS: `ebi-a-GCST009541`
- 来源：HERMES consortium
- 表型：Heart failure
- 样本量：约977,323；病例47,309，对照930,014
- 用途：主分析HF结局

增强数据：

- HERMES更新版HF及亚型summary statistics。
- FinnGen HF相关表型。
- UK Biobank HF相关表型。

## 暴露数据

### 炎症蛋白pQTL

已确认策略：

- 主暴露：UK Biobank Pharma Proteomics Project (UKB-PPP) Olink Explore 3072 pQTL。
- 验证1：deCODE plasma proteomics pQTL，SomaScan平台。
- 验证2：Zhao/SCALLOP 91个Olink Target Inflammation炎症蛋白pQTL。

主分析暴露定义：

优先纳入UKB-PPP中Olink Explore Inflammation和Inflammation_II面板蛋白。已解析得到Inflammation面板368个蛋白、Inflammation_II面板369个蛋白，合计737个主暴露蛋白。

选择理由：

1. UKB-PPP样本量大，覆盖54,219名UK Biobank参与者，测量2,923个独特蛋白。
2. UKB-PPP含明确inflammation/inflammation II面板，适合作为炎症介质发现集。
3. SCALLOP/Zhao 91炎症蛋白与既往AF和HF单病种MR文章高度重叠，不适合作为主暴露。
4. deCODE为泛蛋白组SomaScan资源，适合作为跨平台、独立人群验证。

主暴露来源：

- UKB-PPP / Olink Explore 3072。
- 样本量：54,219。
- 蛋白覆盖：2,941个蛋白分析物，2,923个独特蛋白。
- pQTL概况：14,287个primary genetic associations；约1,954个蛋白具有cis association。
- 数据入口：https://registry.opendata.aws/ukbppp/；https://www.synapse.org/ukbppp

交叉验证来源：

- deCODE plasma proteomics pQTL。
- 样本量：约35,559名冰岛人。
- 平台：SomaScan。
- 覆盖：4,907个aptamer，约4,719个独特蛋白。
- 数据入口：https://www.decode.com/summarydata/

经典炎症面板验证：

- Zhao/SCALLOP Olink Target Inflammation pQTL。
- 样本量：14,824名欧洲ancestry参与者。
- 覆盖：91个血浆炎症相关蛋白。
- pQTL概况：180个pQTL，其中59个cis、121个trans；条件分析后99个cis、128个trans。
- 数据入口：https://www.phpc.cam.ac.uk/ceu/proteins/；GWAS Catalog GCST90274758-GCST90274848。

备选来源：

- GWAS Catalog中已公开的蛋白质水平GWAS。

工具变量建议：

优先使用cis-pQTL，窗口设为编码基因上下游1Mb。这样更接近“药物靶点MR”，也能降低远端水平多效性。trans-pQTL不进入主分析，仅作为补充或机制网络分析。

当前初筛结果：

- UKB-PPP EUR，P < 5e-8。
- 显著pQTL原始关联：4,229条。
- cis-pQTL关联：545条。
- 至少有一个cis-pQTL的蛋白：540个。
- MHC排除后lead cis-pQTL工具变量：529个。
- 主MR输入：`data/processed/ukbppp_inflammation_cis_pqtl_instruments_lead_no_mhc.csv`。

## 辅助数据库

- IEU OpenGWAS：GWAS summary statistics和MR接口。
- GWAS Catalog：检索公开summary statistics。
- FinnGen：疾病结局复制。
- Cardiovascular Disease Knowledge Portal：心血管GWAS和HERMES资源。
- Open Targets：靶点-疾病证据和药物可开发性。
- DrugBank：药物靶点映射。
- STRING/Reactome/KEGG：通路和蛋白互作分析。

## 数据记录字段

每个数据源都要记录：

1. 数据名称。
2. 访问链接。
3. GWAS ID或accession。
4. 样本量。
5. 病例数/对照数。
6. ancestry。
7. 基因组版本。
8. 效应量单位。
9. 下载日期。
10. 使用限制或许可证。
