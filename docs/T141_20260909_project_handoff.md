# PHB_gtdb-ds project handoff (2026-09-09)

## Purpose and authority

This is the single handoff for the work completed through 2026-09-09. It is
written for a new session or a new operator. Read it before making any change.

The project produces **computational candidate homology / functional-potential
evidence**, not experimentally verified PHB or PHA degradation phenotypes.
This boundary applies to HMM, coverage, SignalP6, domain, motif, neighborhood,
taxonomy, co-occurrence, and tree evidence.

Read in this order:

1. `AGENTS.md` - live execution and authorization rules.
2. This handoff - current workstream, decisions, exact locations, and next action.
3. `docs/STATUS.md` - prior repository-wide run-13 authority map. It was last
   updated on 2026-09-03 and therefore does not describe the uncommitted
   2026-09-04 to 2026-09-07 work summarized here.
4. The status records listed in [Authoritative records](#authoritative-records).
5. Each run's `input_contract.json` and result manifest before reusing its data.

When an older report conflicts with a dated accepted run manifest, use the
accepted run manifest and its status record. Never overwrite historical runs or
retroactively modify their contracts.

## Current repository and worktree state

- Workspace: `<REPO_ROOT>`
- Branch: `main`
- Local `HEAD` inspected on 2026-09-09:
  `2e44aafe2c91263447f37bac625f97c56ddafc05`
- Remote configured as: `https://github.com/anzhiw2-gif/PHB_gtdb-ds`
- No `git fetch`, commit, or push was performed during this handoff. Do not
  infer the remote state from this document.
- The 2026-09-04 to 2026-09-07 research scripts, tests, status records, and
  figure outputs are currently **untracked local files**. They are real local
  work products but are not represented by `HEAD` or guaranteed to be on GitHub.
  Preserve them. Commit/push requires a separate explicit authorization.
- The worktree also contains pre-existing untracked archives, oddly named files,
  and LibreOffice uninstall logs. They are outside this handoff's scientific
  scope. Do not delete, move, or normalize them without exact-path confirmation
  and explicit authorization.

## Project decisions that are now fixed

### 1. Main line: only `ePhaZ` and `iPhaZ`

The current Scheme A reporting line retains only two main candidate families:

| Family | Meaning in this project | What it is not |
|---|---|---|
| `ePhaZ` | extracellular PHA-depolymerase-like candidate family, represented by the `ePhaZ_curated_core` strict layer plus a separate broad discovery layer | a PHB-substrate-specific or experimentally confirmed enzyme call |
| `iPhaZ` | intracellular PHA mobilization/depolymerase-like candidate family | evidence that it is extracellular or that a genome has a measured PHB phenotype |

`ePhaZ` and `iPhaZ` were **not assigned from GTDB species labels**. Their seeds
were collected separately from literature/UniProt EC and functional annotation,
plus reference organisms; separate HMMs were then built. Candidate layering
combines family HMM evidence, sequence characteristics, SignalP6 scope, and an
ePhaZ/iPhaZ competition audit. Taxonomy is a downstream description only.

Keep the following layers separate:

| Main-family layer | Proteins | Genomes | Interpretation |
|---|---:|---:|---|
| `ePhaZ_curated_core_secreted` | 21,856 | 17,223 | curated-core candidate with a SignalP6 secretion-related prediction |
| `ePhaZ_curated_core_nonsecreted` | 16,818 | 13,213 | curated-core candidate without that prediction; not a lower-quality phenotype call |
| `ePhaZ_tier2_review` | 31,990 | 25,763 | review candidate, not equivalent to curated core |
| `ePhaZ_competition_review` | 18 | 18 | unresolved ePhaZ/iPhaZ competition case |
| `iPhaZ_tier1` | 32,926 | 25,920 | high-confidence computational intracellular candidate layer |
| `iPhaZ_tier2_review` | 5,479 | 5,153 | review candidate |

The accepted deduplicated main-line table has 109,087 proteins and 76,141
`genome x family` records. The deduplication key is `(family, accession)`;
the accepted output has zero duplicate keys. In the no-signal competition audit,
16,809 are `ePhaZ_like`, 9 are `iPhaZ_like`, and 18 are ambiguous.

**Do not call every main-line object `Tier1`.** In current presentation language,
use `ePhaZ curated-core` and `iPhaZ tier1` exactly. The earlier HMM-rescore
workflow used a curated HMM threshold of `E < 1e-20`, but that numeric rule must
not be generalized as the sole definition of all ePhaZ/iPhaZ figure objects.

### 2. Other families are not part of the current main line

- `ePhaZ_broad_discovery` is a discovery layer and stays outside strict/main
  counts. It can include remote or MCL-PHA-related homologous space.
- `OH` remains governed by `E <= 1e-5` and `min_cov = 0.6`; it is not to be
  forced into the two-family main report.
- `ArchPhaZ_hydrolase` is not an archaeal-main-line family. Its high-coverage
  layer (`coverage >= 0.8`) was entirely bacterial in the audited sample.
- `ArchPhaZ_patatin`, `PhaJ`, `BdhA`, `PhaC`, and `phasin` are broad-fold,
  contextual, or background evidence. None is a standalone depolymerase-positive
  call, and none may be counted as an `ePhaZ`/`iPhaZ` main candidate.

### 3. The archaeal PhaZh1-like workstream is an independent exploratory branch

It must not be merged into the ePhaZ/iPhaZ main statistics or visual workflow.
It is a bounded candidate-only investigation based on a Haloferax PhaZh1-like
patatin route, not evidence for archaeal PHB-depolymerase expansion.

| Archaeal layer | Records | Meaning |
|---|---:|---|
| `PhaZh1_like_high` | 197 records / 173 genomes | review-required candidate-only priority set |
| `PhaZh1_like_review` | 7 | review-required candidate-only |
| `archaeal_patatin_exploratory` | 2,027 | exploratory only |
| `non_PhaZh1_patatin` | 76 | excluded from the branch result |

For the 197 high records, 188 have an N-terminal `GxSxG` motif and 9 do not.
Neighborhood labels are `BdhA` (192), `BdhA;PhaJ` (4), and `PhaJ` (1). These
are context annotations, not activity evidence. The 2026-09-06 registered-HMM
domain check found `ArchPhaZ_patatin` in all 204 high+review inputs, but also
in 5 negative controls. It therefore cannot establish specificity or upgrade a
candidate to a PHB phenotype. Independent Pfam/InterPro annotation is pending.

### 4. Formal scan and historical comparators

- Frozen formal scan `20260901_formal_frozen_scan_13` is complete: 1,000/1,000
  tasks; 6,740,900 registry-accepted rows; four HMMER `>100000 aa` exclusions
  are tool limits, not biological negatives.
- Its strict four-family tier union is 38,741 genomes. This number is a prior
  multi-family reporting result, **not** the current ePhaZ/iPhaZ-only headline.
- Historical Scheme A `44,814` genomes is definition-specific and must never be
  combined with run-13 raw, strict, or current two-family counts.
- No formal scan or registry change was launched for any 2026-09-04 to
  2026-09-07 work. `pipeline/config/formal_scan_models.tsv` remained unchanged.

## Completed work and exact evidence

### A. ePhaZ/iPhaZ evidence layering

- Accepted run: `runs/20260905_ephaz_iphaz_layering_01/`
- Status: `completed_candidate_only`
- Accepted result directory: `results/layers_retry_02/`
- Retained failure evidence: `results/layers_retry_01/`. Do not delete it. The
  first implementation emitted overlapping tier1/tier2 accessions; retry 02
  fixed this with `(family, accession)` deduplication and layer precedence:
  `competition_review > curated/tier1 > tier2 > broad/pending`.
- Status record: `docs/T141_20260905_ephaz_iphaz_layering_01_status.md`

### B. Overall PHA table, conservative PHB focus, and archaeal full audit

- Accepted run: `runs/20260906_overall_pha_phb_archaea_audit_04/`
- Status: `completed_candidate_only`
- Contract: `input_contract.json`
- Result manifest: `results/analysis_manifest.json`
- Main table: `results/overall_pha_phb.tsv`
- Archaeal table: `results/archaeal_phazh1_full_audit.tsv`
- All 109,087 main-line records retain `PHB_candidate_review`; `phb_supported=0`.
  The dual-profile panel had no accession-level match to the GTDB protein IDs,
  so there is no justified PHB-specific upgrade. This is not evidence of absence
  of PHB-active enzymes.
- Earlier `audit_01` through `audit_03` are intermediate/error evidence. The
  `audit_04` run is the accepted version after label compatibility, review-layer,
  and nonsecreted-matching corrections.

### C. Archaeal PhaZh1-like domain and local manual review

- Accepted run: `runs/20260906_archaea_phazh1_domain_annotation_03/`
- Server run/deploy: `${PHB_REMOTE_ROOT}/PHB_gtdb-ds/runs/20260906_archaea_phazh1_domain_annotation_03/`
  and `${PHB_REMOTE_ROOT}/PHB_gtdb-ds/deploy/20260906_archaea_phazh1_domain_annotation_03/`
- Status: `completed_candidate_only`; HMMER 3.4; 18 HMM searches using
  `i-Evalue <= 1e-5`.
- Completion manifest: `results/completion_manifest.json`
- Domain table: `results/domain_annotation/domain_annotation.tsv`
- Control calibration: `results/domain_annotation/positive_negative_calibration.tsv`
- Cross-phylum manual subset: `results/cross_phylum_manual_review.tsv`
- Two earlier attempts are retained failure evidence: `_01` failed on BOM-header
  validation, `_02` on a `phasin` SHA-256 check. Both failed closed before HMMER.
- The manual subset contains 9 records: Halobacteriota review 2, Thermoproteota
  high 2, and Thermoproteota review 5. Its neighborhood labels are all `BdhA`.
  No PhaC/phasin conflict hit was recorded. It is a small audit sample only.

### D. 2026-09-07 progress-report figure suite

- Completed visualization-only run:
  `runs/20260907_progress_report_01/`
- Contract: `input_contract.json` (`status=completed_candidate_only`)
- Renderer: `pipeline/scripts/plot_progress_report_0907.py`
- Output folder: `runs/20260907_progress_report_01/results/0907汇报/`
- Figure manifest: `runs/20260907_progress_report_01/results/figure_manifest.json`
- Source-data folder:
  `runs/20260907_progress_report_01/results/source_data/`
- Five main figures are exported as PDF, PNG, SVG, and TIFF:

| Figure | File stem | What it reports |
|---|---|---|
| 1 | `figure_1_workflow_and_funnel` | GTDB R232 to two-family workflow and data funnel |
| 2 | `figure_2_family_hierarchy` | initial versus strict/main family result hierarchy |
| 3 | `figure_3_signalp_scope` | SignalP6 scope for ePhaZ curated-core only |
| 4 | `figure_4_phylum_distribution` | phylum and genus distributions, stacked by ePhaZ/iPhaZ |
| 5 | `figure_5_genome_representation_and_boundary` | per-genome family representation and competition boundary |

Figure 1 funnel values are 70,682 initial ePhaZ candidates, 38,692 ePhaZ
curated-core candidates, 38,405 initial iPhaZ candidates, 32,926 iPhaZ tier1
candidates, 109,087 deduplicated proteins, and 76,141 `genome x family`
records. It starts from GTDB Release 11 R232 (199,923 genomes), then
Pyrodigal+FAA/GFF QC, DIAMOND blastp (6,769,772), HMMER 3.4 (6,740,900), and
ePhaZ/iPhaZ evidence layering.

> **⚠️ 更正（2026-09-17，经全库实测核实）：本段漏斗描述有一处事实错误。**
> 冻结扫描 `20260901_formal_frozen_scan_13` **没有 DIAMOND 步骤**——该 run 全目录 `grep -ril diamond` **零命中**；其 `results/scan_manifest.json` 记录的是**纯 hmmsearch**：输入 = 100 个预筛分片（`inputs/scan_shards/`，267 GB）、模型 = 10 个 HMM（`inputs/formal_scan_models.tsv`）、任务 = 1,000（100 分片 × 10 模型，`task_completed=1000`）。
> `6,769,772` 不是 DIAMOND 输出行数，而是项目根 `results/logs/rerun_candidates.log` 第 3 行记录的一次**早先 900-tbl 聚合**的命中行数（"聚合完成: 900 个 tbl, 6769772 行命中"）；`6,740,900` 为 registry-accepted 行数（`hits_all.tsv` 原始数据行实测为 6,743,198，见 `scan_manifest.json` 的 `hits_all_sha256`）。正确的漏斗应为：**GTDB R232 (199,923 基因组) → Pyrodigal/FAA/GFF QC 分片 → 10 模型 hmmsearch（1,000 任务）→ ePhaZ/iPhaZ 分层**。
> 证据来源：`runs/20260917_phaded_cys_targeted_recall_01/results/recall_summary.json`（`incidents.upstream_doc_error`）与 `docs/T141_20260917_phaded_cys_targeted_recall_status.md` §3 结论 4。本历史文档原文保留，仅加此更正标注。

Figure 3 used **only** the 38,692 ePhaZ curated-core sequences. SignalP6 results
are SP 16,545; LIPO 5,170; TAT 122; TATLIPO 19; OTHER 16,836. Thus 56.5% have
a recognized secretion-related class and 43.5% are `OTHER`. `OTHER` is not a
secretion route, and localization is not substrate specificity or degradation
activity. SignalP6 also has a `Pilin` class, but it did not occur in this data.

Figure 4's unit is **unique `genome x family` assignment**, not unique genomes
collapsed across families. Before publication or external distribution, update
the x-axis/title language accordingly. Its GTDB placeholder genus `Palsa-465`
should either be labelled explicitly as a placeholder or replaced by the next
formal genus; do not present it as a standard named genus. The panel is a
taxonomy distribution of computational candidates, not a list of experimentally
verified PHA/PHB-degrading organisms.

Figure 5 reports 34,109 genomes with one main family and 9,825 with both.

## Authoritative records

| Question | First record to read |
|---|---|
| project rules and authorization | `AGENTS.md` |
| frozen scan 13 QA | `docs/T141_20260903_run13_qa_02_status.md` |
| tier-processing counts and strict/broad split | `docs/T141_20260902_formal_scan13_tier_processing_02_status.md` |
| family-specific HMM interpretation | `docs/T141_20260903_family_specific_hmm_standards.md` |
| OH / ArchPhaZ cross-phylum boundary | `docs/T141_20260904_oh_arch_crossphylum_audit_status.md` |
| accepted main-family layering | `docs/T141_20260905_ephaz_iphaz_layering_01_status.md` |
| accepted overall PHA/PHB-focus and archaeal audit | `docs/T141_20260906_overall_pha_phb_archaea_audit_04_status.md` |
| archaeal domain audit and retained failures | `docs/T141_20260906_archaea_phazh1_domain_annotation_03_status.md` |
| detailed previous discussion and decisions | `docs/T141_20260906_conversation_process_and_decision_log.md` |
| current visual report | `runs/20260907_progress_report_01/input_contract.json` and `results/figure_manifest.json` |

The primary literature anchors already cited in the decision log are Knoll et
al. 2009 (PMCID `PMC2666664`), Jendrossek and Handrick 2002, Jendrossek et al.
1995 (DOI `10.1128/jb.177.3.596-607.1995`), Papaneophytou et al. 2009 (DOI
`10.1007/s00253-008-1842-2`), and Liu et al. 2015 for archaeal PhaZh1 (DOI
`10.1128/AEM.04269-14`). Do not infer a subtype taxonomy solely from these
anchors; establish a versioned evidence table before creating any classifier.

## Pending work and the only safe next sequence

### Priority 1: define defensible ePhaZ/iPhaZ subtype terminology

The requested next scientific task is a literature-grounded subtype
classification for **only** `ePhaZ` and `iPhaZ`. It is not yet done. The main
open issue is that broad ePhaZ homology contains PHA/MCL-PHA-related space and
should not be relabelled as PHB specific merely by HMM score or SignalP class.

Before creating a classifier or running a new scan:

1. Create a new dated run, for example
   `runs/20260909_ephaz_iphaz_subtype_evidence_01/`, with `inputs/`, `logs/`,
   `results/`, and `input_contract.json`. Do not modify the accepted tables.
2. Freeze an accession-level evidence table for literature/reference sequences:
   accession, exact source database/version/retrieval date, organism, substrate
   tested, experimental assay, cellular context, EC/function text, DOI/PMID,
   positive/negative/control status, and proposed subtype. Record `pending`
   rather than filling unknown fields by inference.
3. Define only subtypes that can be operationalized from that evidence table.
   Keep extracellular/secretory support, intracellular mobilization, MCL-PHA
   challenge, broad/ambiguous, and negative/control categories distinct.
4. Build explicit positive and negative/control panels before any model
   threshold or subtype assignment. Test leakage between training and held-out
   references. A literature label is not interchangeable with a GTDB candidate
   label.
5. Write failing tests first for the input contract, accession uniqueness,
   required evidence fields, exact allowed labels, and output invariants; only
   then implement scripts. Run relevant tests, `compileall`, and `git diff --check`.
6. If a server computation becomes necessary, stop until it is explicitly
   authorized. Bind the exact source snapshot to
   `deploy/<new_run_id>/` and execute only that deploy on T141.

This priority needs a user decision/authorization before it starts because it
creates a new run and may lead to research/network access or server execution.

### Priority 2: editorial repair of the current report figures

No rescan is necessary. When authorized to revise the visual report, make a
new figure-only run rather than overwriting `20260907_progress_report_01`.
The required repairs are:

1. Figure 4: change the x-axis/title to `Tier1 genome-family assignments` (or
   Chinese equivalent), because the input is not a cross-family unique-genome
   count.
2. Figure 4: mark `Palsa-465` as a GTDB placeholder or replace it with the
   next ranked formally named genus.
3. Keep the caption caveat: phylum/genus enrichment reflects GTDB inclusion,
   HMM detection, and annotation bias as well as homology distribution; it is
   not an experimentally validated degraders list.

### Priority 3: optional archaeal deepening, isolated from the main line

Independent Pfam/InterPro evidence for the 197 high + 7 review records is
pending. It may be worthwhile, but no registered independent Pfam/InterPro
database/tool is currently available on T141. Do not present the existing
registered-HMM result as independent domain validation. Any follow-up requires
a new input contract with the candidate FASTA, positive controls, negative
controls, tool/database versions, database checksums, and the planned
candidate-only interpretation. Installation, configuration changes, and server
execution require explicit authorization.

## Prohibited shortcuts

- Do not rerun the formal GTDB scan, modify `formal_scan_models.tsv`, install
  tools, change server configuration, delete files, commit, or push without the
  user's explicit authorization.
- Do not execute a server-root legacy script. Use a dated deploy bound to the
  new run only.
- Do not mix raw accepted hits, strict tier results, historical Scheme A counts,
  ePhaZ/iPhaZ report counts, and genome-family assignments.
- Do not merge archaeal PhaZh1-like records into the ePhaZ/iPhaZ main line.
- Do not turn `SignalP6`, motif, domain, neighborhood, taxonomy, or HMM hits
  into an experimentally verified PHB degradation claim.
- Do not delete failed retries, failed domain runs, overlength exclusions, or
  other negative evidence; they are provenance.

## Read-only resume checklist

Run these commands before any new action to refresh volatile state:

```powershell
Get-Content -Raw AGENTS.md
Get-Content -Raw docs\T141_20260909_project_handoff.md
git status --short
git branch --show-current
git rev-parse HEAD
Get-Content -Raw runs\20260907_progress_report_01\input_contract.json
Get-Content -Raw runs\20260907_progress_report_01\results\figure_manifest.json
```

To verify the accepted local report inputs without modifying them:

```powershell
Get-FileHash -Algorithm SHA256 runs\20260907_progress_report_01\inputs\no_signal_competition.tsv
Get-FileHash -Algorithm SHA256 runs\20260907_progress_report_01\inputs\overall_pha_phb.tsv
Get-FileHash -Algorithm SHA256 runs\20260907_progress_report_01\inputs\tier1_genome_family.tsv
Get-FileHash -Algorithm SHA256 pipeline\scripts\plot_progress_report_0907.py
```

Expected SHA-256 values are respectively
`53e43ebc0c667a714cdc6ba0f3018d894c44852a38b0e7894494dbe0b2613390`,
`2bdbd706a951d3985d422190326f995c2547cd9c87bfebbd29bf48196f624080`,
`35e1f69bca364a0fbc693b0b1abb92a45a765b610bbfe20057be17b0c1964e59`, and
`0dcce982f8c2ca0e42f4b4b9f782ae024632bb7cf53fde46e9c965534ecd469c`.

## Handoff status

- Code/runtime facts: `verified-current` for local files inspected on 2026-09-09.
- Main ePhaZ/iPhaZ report: `completed_candidate_only`.
- Archaeal branch: `completed_candidate_only` for its bounded HMM/manual audit;
  independent Pfam/InterPro remains `pending`.
- New ePhaZ/iPhaZ subtype classification: `pending`, not started.
- Figure 4 terminology/placeholder correction: `pending` editorial revision.
- Formal rescan, server-side work, installation, commit, push, and deletion:
  `out-of-scope unless explicitly authorized`.

