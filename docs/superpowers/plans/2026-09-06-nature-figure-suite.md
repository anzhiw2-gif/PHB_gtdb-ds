# Nature Figure Suite Implementation Plan

> **For agentic workers:** Execute this plan inline in the current session with Python/Matplotlib.

**Goal:** Generate six main and two supplementary publication-style figures from accepted project outputs, with traceable source data and candidate-only evidence boundaries.

**Architecture:** A single Python renderer reads the accepted run tables and recorded manifests, writes normalized source-data TSVs, and exports one directory per figure. Each figure is independently reproducible and receives SVG, PDF, 600-dpi TIFF/PNG, and a QA record.

**Tech Stack:** Python 3.12, pandas, matplotlib, numpy; SVG text remains editable and PDF uses TrueType fonts.

## Global Constraints

- Use only accepted project outputs; do not fabricate counts or mix historical Scheme-A and current run-13 definitions.
- Keep HMM, domain, neighborhood, taxonomy and motif evidence candidate-only; no phenotype confirmation language.
- Create a new run under `runs/20260906_nature_figure_suite_01/` with `logs/`, `inputs/`, `results/`, and `input_contract.json`.
- Do not modify `pipeline/config/formal_scan_models.tsv`; do not run formal scans; do not commit or push.
- Use Python/Matplotlib exclusively for drawing, previewing, exporting and QA.
- Export editable SVG, PDF, 600-dpi TIFF, and PNG preview for every figure.

## Figure Set

1. Workflow and completion status.
2. Overall PHA candidate subtype distribution.
3. Core family candidate scale and genome coverage.
4. Phylum-level family distribution.
5. Archaeal PhaZh1-like audit and motif distribution.
6. Species/genome distribution and manual cross-phylum review.
7. Supplement: positive/negative HMM calibration.
8. Supplement: neighborhood-marker association.

## Verification

- Source-data counts match plotted values.
- Every quantitative panel declares `n` and candidate-only interpretation in metadata.
- SVG contains editable text nodes; files are nonempty and dimensions are recorded.
- Run `python -m pytest pipeline/tests -v`, `python -m compileall -q pipeline/scripts pipeline/tests`, and `git diff --check`.
