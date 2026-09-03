# PHB_gtdb-ds Current Status (2026-09-02; reviewed 2026-09-03)

Repository snapshot: GitHub `main` (reviewed 2026-09-03).

This file supersedes stale running-state statements in older audit snapshots. Historical reports and runs are retained for provenance.

## Executive status

- Core Scheme A candidate screening is complete.
- The formal frozen scan `20260901_formal_frozen_scan_13` on T141 is **completed**.
- ePhaZ production models remain `ePhaZ_curated_core` and `ePhaZ_broad_discovery` with `E<=1e-5` and no generic coverage gate.
- MCL-PHA subfamily profiles are **candidate-only** and are not registered for production scanning.
- Biological interpretation remains **candidate homology/function potential**, not experimentally verified PHB degradation.

## Formal scan 13 acceptance

Evidence: `${PHB_REMOTE_ROOT}/PHB_gtdb-ds/runs/20260901_formal_frozen_scan_13/`.

- `scan_manifest.json`: `status=completed`
- 60 threads; 1000/1000 tasks completed; 297 tasks reused from the parent run
- `logs/failed_tasks.tsv`: header only (no failed task rows)
- `results/hits_all.tsv`: 6,743,198 lines including the header
- `results/overlength_exclusions.tsv`: 4 records, all HMMER `>100000 aa` tool-limit exclusions
- Input contract: `status=verified`, GTDB R232 taxonomy/metadata/tree and model hashes recorded

The four overlength records are tool exclusions, not biological negatives.

## ePhaZ decision

Retain the current production registry unchanged. The MCL candidate work recovered the two previously missed divergent positives and showed that a candidate coverage gate of `0.8` separates the bounded challenge panel, but the divergent profiles are singleton models without independent validation. Do not add them to `pipeline/config/formal_scan_models.tsv` yet.

## Family interpretation

- `iPhaZ`: small-panel separation is good; retain ePhaZ/iPhaZ competition auditing.
- `OH`: retain `min_cov=0.6` to suppress nylon-hydrolase false positives; keep `Q4W8C9` as a boundary case.
- `ArchPhaZ_hydrolase`: usable candidate family.
- `ArchPhaZ_patatin`: broad patatin/lipid-enzyme fold; require neighborhood evidence and do not call it PHB-specific.
- `PhaJ`, `phasin`, and `PhaC`: contextual/auxiliary markers, not standalone degradation calls.
- `BdhA`: background metabolism family, excluded from core depolymerase counts.

## Follow-up work

1. Completed run-13 downstream tier processing in `20260902_formal_scan13_tier_processing_02`; strict four-family union is 38,741 genomes. See `docs/T141_20260902_formal_scan13_tier_processing_02_status.md`.
2. Historical run-12 snapshots are marked as launch-time records; keep them for provenance and use this file for current status.
3. Keep tree/HGT work paused until the relevant inputs are rebuilt and registered.
4. Prepare publication outputs only after scan-13 downstream results are reviewed.

## Governance decision

The ePhaZ repair loop is closed for now. Further small-panel threshold/LOO runs require a new independent evidence set and a concrete failure mode; otherwise they are not scheduled.
