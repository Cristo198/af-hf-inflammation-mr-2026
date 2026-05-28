# Results段落初稿：敏感性分析、反向MR和复制状态

说明：以下仅为Results中“敏感性分析和反向MR”小节草稿，不包含Introduction或Discussion。

由于本研究主分析采用每种蛋白一个lead cis-pQTL的单工具变量设计，MR-Egger、weighted median、MR-PRESSO和leave-one-out不适用于主估计。可行的敏感性检查显示，899条协调后的主MR记录中无F统计量≤10的弱工具变量；141条记录涉及回文等位基因，但FGF5和LPA候选工具变量均非回文等位基因。基于观测尺度的proxy-Steiger方向性检查显示，896/899条主MR记录的暴露R² proxy大于结局R² proxy，支持蛋白到疾病方向；FGF5/LPA候选4条主结果均支持蛋白到疾病方向。需要说明的是，该方向性检查对二分类结局使用观测log OR尺度近似，不能完全替代liability-scale Steiger检验。

候选反向MR以本地主结局GWAS中P < 5 × 10^-8的AF和HF变异作为疾病遗传工具变量，并在缺乏LD参考面板的情况下按10 Mb距离剪枝。AF主结局获得75个距离剪枝工具变量，HF主结局获得11个距离剪枝工具变量；与FGF5/LPA全基因组UKB-PPP pGWAS协调后，AF相关反向MR纳入74个工具变量，HF相关反向MR纳入9个工具变量。全距离剪枝集显示AF遗传易感性与FGF5水平存在名义显著正向关联（beta=0.081，SE=0.011，P=1.40 × 10^-13），HF遗传易感性与LPA水平亦显示显著正向关联（beta=0.925，SE=0.042，P值在数值计算中下溢为0）。但在剔除目标蛋白cis区域±1 Mb内的疾病工具变量后，HF->LPA反向信号消失（beta=-0.0066，SE=0.0437，P=0.879），提示全量结果主要可能由LPA局部信号驱动；AF->FGF5仍保留名义显著（beta=0.024，SE=0.011，P=0.0278），但异质性仍较高，且本分析未进行LD参考面板严格clumping，因此只能作为潜在双向性或多效性线索，而不能作为确定性反向因果证据。HF->FGF5和AF->LPA未显示反向MR支持。

复制分析方面，FGF5/LPA的FinnGen R12候选精确变异复制已完成：FGF5在FinnGen AF中方向一致并显著，LPA在FinnGen AF和strict HF中方向一致，其中strict HF达到名义显著。进一步的FinnGen R12全量复制也已完成，两个结局各协调496个工具变量。UKB/OpenGWAS补充复制已完成，但对本研究529个lead cis-pQTL覆盖很低，两个结局各仅协调7个工具变量，且未覆盖FGF5/LPA共享候选；同时UKB结局复制还需额外评估与UKB-PPP暴露样本的重叠问题。
