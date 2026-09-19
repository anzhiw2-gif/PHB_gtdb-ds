# PHB-GTDB 项目规则

- 运行必须使用新的 `runs/<run_id>/`；不得覆盖历史 `results/` 或服务器历史目录。
- `run_id` 必须可审计且不含路径遍历；运行目录保存 `logs/`、`inputs/`、`results/` 与 `input_contract.json`。
- GTDB taxonomy、metadata、tree、HMM、源码和环境必须记录路径、版本、大小及 SHA-256；缺失项写 `pending`，不得伪造哈希。
- 服务器只允许执行 dated `deploy/<run_id>/` 中已绑定源码；不得直接运行服务器根目录旧脚本。
- 服务器计算单任务线程上限 **40**；启动前必须按实际空闲动态取 `min(40, 总核数 − 负载 − 10)` 并保留约 10 核余量，不得让服务器超载（总核数以 `nproc`、负载以 `/proc/loadavg` 1 分钟值实测为准）。
- 本项目的 HMM、domain、SignalP、邻域和树结果表示候选同源或功能潜力，不等同于已验证 PHB 降解表型。
- 允许用 DED 家族序列（含 `annotation_only`）训练**发现层 HMM**，层名固定为 `discovery_hmm_uncalibrated` 且每份输出必须带该标注；发现层只用于**召回**，**不得**用于筛选、删除或降级任何候选（实测其对 563 条 x₁ 混淆嫌疑命中 443 条 = 78.69%，无判别力）。
- 发现层不得写入 `pipeline/config/formal_scan_models.tsv`、不得产生 family 判定；它只能对既有中间产物打分，任何 GTDB 全库重扫仍需**单独授权**。
- with-lipase（`DED_hfam_2`）池外命中是已知噪声（2026-09-19 实测 1,206,655 条），必须保留为 **deferred 结构验证层**（`runs/20260919_phaded_with_lipase_deferred_archive_01`）而**禁止删除**；该层不进入任何候选/高可信度计数，后续只能经结构证据（PDB 8YNV Foldseek 或结构预测比对）召回，通过者再走证据层与校准 gate。
- HMM/profile 必须绑定其训练比对与输入集合的 SHA-256；比对工具不可比特复现时（实测 MAFFT 7.525 同一输入两次得到不同比对），必须把**具体比对哈希**写进 manifest 并声明不可位复现，不得声称逐位可复现。
- 保留历史运行残留和失败证据；清理前先确认精确路径与可恢复性。
- 服务器 `runs/` 实测已达 1,016 GB：已被取代的历史 run 可移入归档目录并保留原名、manifest 与 SHA-256 清单，**禁止删除**。
- 代码改动先写失败测试，再实现；完成后运行相关测试、`compileall` 和 `git diff --check`。
- 本地 `git commit` 允许，但每次必须在交接文档记录；**push 到 GitHub 仍需操作者明确授权**。
- 重算、安装、配置变更、删除与推送仍须明确授权；只读审计应保持只读。
- 发布前核对本地、GitHub、服务器 deploy 和 run manifest 的权威关系，避免混用版本。
