# T141 PhaDED 全库 InterProScan 补充任务

日期：2026-09-14  
运行：`20260914_phaded_interpro_full_01`  
状态：InterProScan 结果已生成；封装脚本因相对路径错误退出，结果已恢复到本地 dated run

## 为什么现在启动

上一轮 `20260914_phaded_full_evidence_reconciliation_01` 只是本地只读合并已有结果，服务器不会增加信息，因此没有启动计算。当前全库只有 244 条优先候选完成 InterProScan；为满足“尽可能完整”的序列层面资料需求，本轮对全库候选执行 InterProScan 二次域注释。

## 执行边界

- 仅运行 candidate-only InterProScan；不修改 registry，不进行正式 GTDB 重扫，不涉及湿实验。
- 线程固定为 40；GPU 不使用。
- 输入 FASTA 去除单一末端 `*`；含内部异常字符的 `351` 条记录写入 `interpro_excluded.tsv`，保留在候选台账中，不视为阴性。
- 服务器只执行 `deploy/20260914_phaded_interpro_full_01/` 中绑定的源码和输入。

## 当前运行

- 服务器：`<SERVER_USER>@<SERVER_HOST>`（`grius`）
- InterProScan：`5.78-109.0`
- PID：`2617348`
- 输入：`108,736` 条规范化序列
- 日志：`runs/20260914_phaded_interpro_full_01/logs/`
- 输出：`runs/20260914_phaded_interpro_full_01/results/interpro_full.tsv`

## 实际结果

- InterProScan 子任务已完成，结果原先写入软件安装目录下的相对路径：`${PHB_REMOTE_ROOT}/software/interproscan-5.78-109.0/interproscan-5.78-109.0/runs/20260914_phaded_interpro_full_01/results/interpro_full.tsv`。
- 结果已复制回本地 `results/interpro_full.tsv`：`671,556` 行、`108,735` 个 accession、`126,215,214` bytes，SHA-256 `6ef482ede8d1ab4bb93840468d5cdd4e883116618a5ac2455f63001c15e6b437`。
- 与原始 `109,087` 条候选相比缺失 `352` 条：`351` 条在工具输入阶段因内部异常字符排除，另 `1` 条无 InterPro 报告；两类都不是生物学阴性。
- 封装脚本错误日志保留在 `logs/interproscan.stderr.log`；没有重跑同一长任务。

## 后续

修复了本地脚本的相对路径问题，并增加回归测试。下一步是验证/解析已恢复的 TSV，再与全库 PhaDED/Pfam 台账合并并重新计算证据等级；未命中或未返回的记录标记为 `not_tested/unresolved`。
