# PHB-GTDB 项目规则

- 运行必须使用新的 `runs/<run_id>/`；不得覆盖历史 `results/` 或服务器历史目录。
- `run_id` 必须可审计且不含路径遍历；运行目录保存 `logs/`、`inputs/`、`results/` 与 `input_contract.json`。
- GTDB taxonomy、metadata、tree、HMM、源码和环境必须记录路径、版本、大小及 SHA-256；缺失项写 `pending`，不得伪造哈希。
- 服务器只允许执行 dated `deploy/<run_id>/` 中已绑定源码；不得直接运行服务器根目录旧脚本。
- 本项目的 HMM、domain、SignalP、邻域和树结果表示候选同源或功能潜力，不等同于已验证 PHB 降解表型。
- 任何重算、安装、配置变更、删除、提交或推送都必须得到明确授权；只读审计应保持只读。
- 保留历史运行残留和失败证据；清理前先确认精确路径与可恢复性。
- 代码改动先写失败测试，再实现；完成后运行相关测试、`compileall` 和 `git diff --check`。
- 发布前核对本地、GitHub、服务器 deploy 和 run manifest 的权威关系，避免混用版本。
