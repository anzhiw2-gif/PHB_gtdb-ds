# T141 ePhaZ dual-profile candidate 06

服务器 dated deploy `20260902_ephaz_dual_profile_candidate_06` 已完成。
本次只构建并校准候选 MCL-PHA HMM，与既有 PHB-focused HMM 做双 profile
分类；未修改正式 registry、run 13 或全库扫描。

- MCL seed：`Q51718`, `Q6UFW4`, `Q84C08`；其中 Q84C08 是明确不降解 PHB 的 MCL-PHA 对照。
- 阈值：`E<=1e-5`；双命中 bitscore margin `10`；CPU `8`。
- 完整 probe 输出：`7 PHB_like`, `2 MCL_like`, `8 no_hit`。
- 五条实验 PHB 阳性为 `5/5 PHB_like`；Q6UFW4 与 Q84C08 均为 `MCL_like`。

Q5Y152 不作为本次 MCL seed；其 UniProt accession 为未审阅注释记录，不能
替代 accession 级实验阳性。HMM 命中仍表示序列层面的家族/功能潜力。

完整哈希见 run 的 `results/output_sha256.tsv` 和部署清单。
