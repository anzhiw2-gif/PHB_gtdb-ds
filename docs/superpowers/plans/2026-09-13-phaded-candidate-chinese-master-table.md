# PhaDED Candidate Chinese Master Table Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Produce an auditable Chinese summary table for the 19 manually reviewed PhaDED candidates.

**Architecture:** Copy no raw sequences and do not rerun classification. Bind the manual-review table, its 14-record priority table, and the completed full-metadata literature bridge as immutable inputs. Publish an aligned TSV and readable Markdown table in a new dated run.

**Tech Stack:** UTF-8 TSV, Markdown, SHA-256 input/output contract, existing frozen run artifacts.

## Global Constraints

- Do not promote homology, Pfam/InterPro, Foldseek, structure, or phylogeny evidence to PHB phenotype evidence.
- Exact literature-to-GTDB bridge is required for direct literature support.
- Preserve the 3 architecture holds and 2 source/structure holds.
- No wet-lab work, GTDB-wide rescan, formal registry change, deploy, commit, push, or deletion.

### Task 1: Create the Chinese candidate master table

**Files:**
- Create: `runs/20260913_phaded_candidate_chinese_master_table_01/input_contract.json`
- Create: `runs/20260913_phaded_candidate_chinese_master_table_01/results/phaded_candidate_chinese_master_table.tsv`
- Create: `docs/T141_20260913_phaded_candidate_chinese_master_table_01.md`

- [x] Bind the three frozen source artifacts with path, size, and SHA-256.
- [x] For each of the 19 candidate IDs, preserve PhaDED classification, nearest reference family, evidence, review decision, retention status, and exact literature-bridge status.
- [x] State explicitly that no individual record among these 19 has an exact literature bridge; do not treat this as a biological negative.
- [x] Add a group summary and glossary for the 14 retained and 5 hold records.
- [x] Verify JSON parsing, TSV row count and headers, file hashes, and `git diff --check`.
