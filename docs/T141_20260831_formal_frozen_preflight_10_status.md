# T141 Formal Frozen Screen Preflight

## Run

- Run ID: `20260831_formal_frozen_preflight_10`
- Deploy: `deploy/20260831_formal_frozen_preflight_10/`
- Server run: `${PHB_REMOTE_ROOT}/PHB_gtdb-ds/runs/20260831_formal_frozen_preflight_10/`
- Status: `planned_not_run`
- Executed: 2026-08-31 (Asia/Shanghai; manifest UTC `2026-08-30T16:34:28Z`)

## Verified

- 10 registry models copied and SHA-256 verified.
- Frozen thresholds are `E=1e-5`; `OH` retains `min_cov=0.6`; other models retain registry coverage values.
- `iPhaZ` resolves from `frozen_data_root` and has the calibrated SHA-256.
- 100 protein shards were enumerated and hashed; the complete ledger is in the run `inputs/shards.tsv`.
- GTDB taxonomy, metadata, and tree were verified in `input_contract.json`.
- T141 deploy script passed `bash -n`; environment reports HMMER 3.4.

## Boundary

No `hmmsearch` command was executed. The run contains no `tblout`, `domtblout`, `hits_all.tsv`, or scan manifest. HMM/domain matches remain candidate homology or function-potential evidence and are not phenotype proof.

The 268G formal scan remains pending explicit authorization. If authorized, it must use this dated deploy and run a new run ID; it must not reuse or overwrite this preflight directory.
