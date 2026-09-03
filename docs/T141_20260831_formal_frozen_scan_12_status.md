# T141 Formal Frozen Scan 12

> **历史启动快照**：本页记录 run-12 启动时的状态。后续完成情况以 [CURRENT_STATUS_20260902.md](CURRENT_STATUS_20260902.md) 和 run-13 文档为准。

## Launch

- Run ID: `20260831_formal_frozen_scan_12`
- Deploy: `deploy/20260831_formal_frozen_scan_12/`
- Server run: `${PHB_REMOTE_ROOT}/PHB_gtdb-ds/runs/20260831_formal_frozen_scan_12/`
- Status at launch: running
- Launch mode: detached `setsid`, one HMMER process per model/shard

## Inputs and models

The run uses the frozen registry and the same 100 GTDB protein shards as the completed preflight. The ten frozen HMMs are copied into the run directory and bound by `input_contract.json`. Thresholds remain `E=1e-5`; `OH` retains `min_cov=0.6`.

## HMMER limit handling

Before each shard is searched, `filter_hmmsearch_shard.py` removes only records longer than 100,000 aa, the HMMER comparison-pipeline limit. Each excluded accession and observed length is appended to `results/overlength_exclusions.tsv`. These exclusions are tool-limit records, not biological negatives or phenotype calls.

## Evidence boundary

The run is not complete until all 1,000 model/shard tasks finish, failed task log is empty, `hits_all.tsv` is nonempty, and the final scan manifest is written. HMM/domain hits remain candidate homology or function-potential evidence and do not establish PHB degradation phenotype.
