# Scheme A Figure Suite Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Generate six English, double-column scientific figures that present the current Scheme A candidate-homolog workflow, seed provenance, classification distribution, SignalP evidence, and genomic-neighborhood context.

**Architecture:** A Python data-preparation script will convert immutable project tables and a downloaded current cluster snapshot into compact, figure-specific TSV source files plus provenance metadata. A separate Python-only matplotlib renderer will consume only those TSV files and produce editable SVG, PDF, 600 dpi TIFF, and PNG previews. A QA script will confirm expected source-data totals, required SVG labels, file dimensions, and exports.

**Tech Stack:** Python 3.11+, pandas, numpy, matplotlib, csv, json, hashlib, unittest, PowerShell/OpenSSH for the current server snapshot.

## Global Constraints

- Use Python only for data preparation, plotting, previews, exports, and visual QA.
- Use English figure text and double-column width (183 mm) for every figure.
- Use only candidate-homolog, signal-prediction, taxonomy, and neighborhood-context language; do not claim enzyme activity, phenotype, pathway completeness, co-transcription, or causation.
- Keep Figure 5 as a current-progress snapshot without a visible remediation badge; its provenance metadata must identify its source snapshot.
- Exclude phylogenetic trees from this figure suite.
- Treat counts as database census-style values; do not draw p-values, confidence intervals, or replicate error bars.
- Export editable SVG/PDF plus 600 dpi TIFF and PNG previews.
- Do not commit, push, install packages, delete existing results, or overwrite legacy figures without explicit user authorization.

---

## File Structure

- Create: `pipeline/scripts/build_scheme_a_figure_data.py`
  - Extracts compact source TSVs and `provenance.json` from seed, tier, SignalP, and cluster tables.
- Create: `pipeline/scripts/plot_scheme_a_figures.py`
  - Renders Figures 0-5 from the compact source TSVs using matplotlib only.
- Create: `pipeline/scripts/check_scheme_a_figures.py`
  - Validates source-data invariants, export presence, SVG editability, and required text labels.
- Create: `pipeline/tests/test_scheme_a_figure_data.py`
  - Unit tests pure extraction and transformation functions using temporary fixture TSVs.
- Create: `pipeline/tests/test_scheme_a_figure_exports.py`
  - Runs the renderer against fixture source data and checks all four export formats.
- Create: `pipeline/tests/test_scheme_a_figure_qa.py`
  - Exercises validation failures independently of renderer fixture tests.
- Create: `results/figures/scheme_a/source_data/`
  - Machine-readable figure inputs generated at runtime; do not hand-edit.
- Create: `results/figures/scheme_a/`
  - Runtime figure exports and PNG previews; do not overwrite existing legacy figures.

## Data Contracts

`build_scheme_a_figure_data.py` accepts these paths:

```text
--seed-manifest pipeline/seeds/seeds_manifest.tsv
--seed-stats pipeline/seeds/seeds_stats.json
--tier-table results/tables/tier1_genome_family.tsv
--tier-phylum results/tables/tier1_phylum_distribution.tsv
--tier-fasta-dir data/screen/tiers
--signal-summary results/tables/signalp_ePhaZ_analysis.tsv
--signal-phylum results/tables/ePhaZ_signalpeptide_phylum.tsv
--cluster-context <current-server-snapshot>/cluster_context.tsv
--cluster-audit <current-server-snapshot>/cluster_locus_audit.tsv
--cluster-summary <current-server-snapshot>/cluster_summary.tsv
--output-dir results/figures/scheme_a/source_data
```

It produces:

```text
figure0_seed_family.tsv
figure0_seed_evidence.tsv
figure0_seed_split.tsv
figure1_funnel.tsv
figure1_genome_coverage.tsv
figure2_core_scale.tsv
figure2_union.tsv
figure3_phylum_matrix.tsv
figure3_phylum_totals.tsv
figure4_signal_classes.tsv
figure4_signal_phyla.tsv
figure5_neighborhood_rate.tsv
figure5_patatin_context.tsv
provenance.json
```

All table writers use UTF-8, tab delimiters, deterministic sort order, and a header row. `provenance.json` records input path labels, SHA-256, row count, extraction time, and the current cluster snapshot identifier.

### Task 1: Build Deterministic Figure Source Data

**Files:**
- Create: `pipeline/scripts/build_scheme_a_figure_data.py`
- Create: `pipeline/tests/test_scheme_a_figure_data.py`

**Interfaces:**
- Consumes: the data-contract paths above.
- Produces: `build_figure_data(args) -> dict[str, pathlib.Path]`, one compact TSV per figure, and `provenance.json`.
- Pure helpers: `count_fasta_headers(path) -> int`, `core_union(tier_df) -> int`, `extract_seed_tables(seed_df) -> tuple[pd.DataFrame, ...]`, `neighborhood_rates(context_df, audit_df) -> pd.DataFrame`.

- [x] **Step 1: Write failing source-data tests**

```python
def test_core_union_excludes_patatin_and_counts_unique_genomes(tmp_path):
    tier = pd.DataFrame([
        {"genome": "G1", "family": "ePhaZ"},
        {"genome": "G1", "family": "iPhaZ"},
        {"genome": "G2", "family": "OH"},
        {"genome": "G3", "family": "ArchPhaZ_patatin"},
    ])
    assert figure_data.core_union(tier) == 2


def test_neighborhood_rate_uses_unique_hit_loci_as_denominator():
    audit = pd.DataFrame([
        {"genome": "G1", "locus": "L1", "family": "ePhaZ", "status": "analyzed"},
        {"genome": "G1", "locus": "L2", "family": "ePhaZ", "status": "analyzed"},
    ])
    context = pd.DataFrame([
        {"genome": "G1", "hit_locus": "L1", "hit_family": "ePhaZ", "marker_family": "BdhA"},
        {"genome": "G1", "hit_locus": "L1", "hit_family": "ePhaZ", "marker_family": "BdhA"},
    ])
    result = figure_data.neighborhood_rates(context, audit)
    assert result.loc[0, "candidate_loci"] == 2
    assert result.loc[0, "supported_loci"] == 1
    assert result.loc[0, "support_rate"] == 0.5
```

- [x] **Step 2: Run tests and verify they fail**

Run:

```powershell
python -m unittest pipeline.tests.test_scheme_a_figure_data -v
```

Expected: failure because `build_scheme_a_figure_data.py` and the tested helpers do not exist.

- [x] **Step 3: Implement source-data helpers and CLI**

```python
CORE_FAMILIES = ("ePhaZ", "iPhaZ", "OH", "ArchPhaZ_hydrolase")


def core_union(tier_df: pd.DataFrame) -> int:
    return tier_df.loc[tier_df["family"].isin(CORE_FAMILIES), "genome"].nunique()


def neighborhood_rates(context_df: pd.DataFrame, audit_df: pd.DataFrame) -> pd.DataFrame:
    analyzed = audit_df.loc[audit_df["status"] == "analyzed", ["genome", "locus", "family"]]
    denominators = analyzed.groupby("family").size().rename("candidate_loci")
    supported = (
        context_df[["genome", "hit_locus", "hit_family", "marker_family"]]
        .drop_duplicates()
        .groupby(["hit_family", "marker_family"])
        .size()
        .rename("supported_loci")
        .reset_index()
    )
    supported["candidate_loci"] = supported["hit_family"].map(denominators)
    supported["support_rate"] = supported["supported_loci"] / supported["candidate_loci"]
    return supported.sort_values(["hit_family", "marker_family"], kind="stable")
```

Implement SHA-256 streaming with 1 MiB chunks. Reject missing, empty, malformed, or incomplete `cluster_locus_audit.tsv` inputs. For Figure 1, keep explicit fields `stage`, `value`, `unit`, and `note` so hit rows are never visually compared as if they were unique proteins.

- [x] **Step 4: Run tests and verify they pass**

Run:

```powershell
python -m unittest pipeline.tests.test_scheme_a_figure_data -v
```

Expected: all source-data tests pass.

- [x] **Step 5: Generate compact source data from current project tables**

First retrieve the current server cluster snapshot into a new, dated local staging directory without overwriting legacy local cluster tables. Then run:

```powershell
python pipeline/scripts/build_scheme_a_figure_data.py `
  --seed-manifest pipeline/seeds/seeds_manifest.tsv `
  --seed-stats pipeline/seeds/seeds_stats.json `
  --tier-table results/tables/tier1_genome_family.tsv `
  --tier-phylum results/tables/tier1_phylum_distribution.tsv `
  --tier-fasta-dir data/screen/tiers `
  --signal-summary results/tables/signalp_ePhaZ_analysis.tsv `
  --signal-phylum results/tables/ePhaZ_signalpeptide_phylum.tsv `
  --cluster-context results/figures/scheme_a/staging/cluster_context.tsv `
  --cluster-audit results/figures/scheme_a/staging/cluster_locus_audit.tsv `
  --cluster-summary results/figures/scheme_a/staging/cluster_summary.tsv `
  --output-dir results/figures/scheme_a/source_data
```

Expected: eleven TSVs and one `provenance.json`, all non-empty.

### Task 2: Render The Six Publication Figures

**Files:**
- Create: `pipeline/scripts/plot_scheme_a_figures.py`
- Create: `pipeline/tests/test_scheme_a_figure_exports.py`

**Interfaces:**
- Consumes: `results/figures/scheme_a/source_data/*.tsv`.
- Produces: `figure_0_seed_library`, `figure_1_workflow_funnel`, `figure_2_core_scale`, `figure_3_phylum_distribution`, `figure_4_ephaz_signal`, and `figure_5_neighborhood_context`, each in SVG, PDF, TIFF, and PNG.
- Public renderers: `render_figure_0(data_dir, output_dir) -> list[Path]` through `render_figure_5(data_dir, output_dir) -> list[Path]`.

- [x] **Step 1: Write failing export tests**

```python
def test_renderer_exports_editable_svg_and_all_requested_formats(tmp_path):
    write_minimal_valid_source_data(tmp_path / "source_data")
    outputs = figures.render_figure_2(tmp_path / "source_data", tmp_path / "out")
    assert {path.suffix for path in outputs} == {".svg", ".pdf", ".tiff", ".png"}
    svg = (tmp_path / "out" / "figure_2_core_scale.svg").read_text(encoding="utf-8")
    assert "<text" in svg
    assert "44,814" in svg
```

- [x] **Step 2: Run tests and verify they fail**

Run:

```powershell
python -m unittest pipeline.tests.test_scheme_a_figure_qa -v
```

Expected: failure because the renderer module and figure function do not exist.

- [x] **Step 3: Implement shared style and all renderers**

At the top of `plot_scheme_a_figures.py`, before figure creation:

```python
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
```

Use 183 mm width (`7.205` inches), a white background, lower-case panel labels,
and stable colors: ePhaZ `#0F4D92`, iPhaZ `#42949E`, OH `#B64342`,
ArchPhaZ_hydrolase `#9A4D8E`, and neutral gray `#767676`.

Figure-specific layout:

```text
Figure 0: three aligned panels: family review status, train/validation split,
          evidence provenance.
Figure 1: left 55% flow boxes and arrows; right 45% horizontal funnel bars.
Figure 2: sequence-scale bars, genome-scale bars, and one union-coverage panel.
Figure 3: top-15 phylum total bars above/alongside a log10(count + 1) core-family matrix.
Figure 4: SignalP-class horizontal bars and leading-phylum stacked bars.
Figure 5: family-by-marker support-rate heatmap plus patatin locus/genome context bars.
```

All Figure 5 labels say `genomic neighborhood support (+/-10 kb)` and do not
use functional verbs. Add no visible remediation badge. Save each figure through
a shared `save_exports(fig, stem)` helper that creates SVG, PDF, TIFF at 600 dpi,
and PNG at 300 dpi.

- [x] **Step 4: Run export tests and verify they pass**

Run:

```powershell
python -m unittest pipeline.tests.test_scheme_a_figure_exports -v
```

Expected: all fixture export tests pass and the SVG includes editable `<text>` nodes.

- [x] **Step 5: Render the approved suite**

Run:

```powershell
python pipeline/scripts/plot_scheme_a_figures.py `
  --data-dir results/figures/scheme_a/source_data `
  --output-dir results/figures/scheme_a
```

Expected: six figure stems, each exported in SVG, PDF, TIFF, and PNG.

### Task 3: Validate Figure Values And Visual Deliverables

**Files:**
- Create: `pipeline/scripts/check_scheme_a_figures.py`
- Modify: `docs/STATUS.md` only if all figure values and source-data provenance checks pass.

**Interfaces:**
- Consumes: compact source TSVs, `provenance.json`, and the 24 export files.
- Produces: `results/figures/scheme_a/qa_report.json` and a concise console summary.

- [x] **Step 1: Write the failing QA test**

```python
def test_qa_rejects_svg_without_editable_text(tmp_path):
    svg = tmp_path / "broken.svg"
    svg.write_text("<svg><path /></svg>", encoding="utf-8")
    with self.assertRaises(figure_qa.FigureQAError):
        figure_qa.check_svg_text(svg)
```

- [x] **Step 2: Run the test and verify it fails**

Run:

```powershell
python -m unittest pipeline.tests.test_scheme_a_figure_exports -v
```

Expected: failure because `check_scheme_a_figures.py` is absent.

- [x] **Step 3: Implement QA checks**

```python
def check_svg_text(path: Path) -> None:
    if "<text" not in path.read_text(encoding="utf-8"):
        raise FigureQAError(f"SVG has no editable text nodes: {path}")


def check_required_exports(output_dir: Path, stem: str) -> None:
    for suffix in (".svg", ".pdf", ".tiff", ".png"):
        path = output_dir / f"{stem}{suffix}"
        if not path.is_file() or path.stat().st_size == 0:
            raise FigureQAError(f"missing or empty export: {path}")
```

Also verify: Figure 2 union equals 44,814; Figure 2 denominator equals 199,923;
Figure 3 excludes `ArchPhaZ_patatin`; Figure 5 marker set equals
`PhaC,PhaE,PhaJ,BdhA,phasin,PHA_gran_rgn`; PNG dimensions are nonzero; every
source-data file listed in `provenance.json` still hashes to its recorded value.

- [x] **Step 4: Run QA tests and verify they pass**

Run:

```powershell
python -m unittest pipeline.tests.test_scheme_a_figure_data pipeline.tests.test_scheme_a_figure_exports pipeline.tests.test_scheme_a_figure_qa -v
```

Expected: all figure tests pass.

- [x] **Step 5: Run final QA and inspect Python-rendered previews**

Run:

```powershell
python pipeline/scripts/check_scheme_a_figures.py `
  --data-dir results/figures/scheme_a/source_data `
  --output-dir results/figures/scheme_a `
  --report results/figures/scheme_a/qa_report.json
```

Inspect the six PNG previews using the Python-produced raster exports. Confirm
at final size: no clipping or overlap, legible labels, consistent family colors,
no rainbow colormap, no unsupported functional claim, and readable panel letters.

Expected: `qa_report.json` records `status: pass` for all six figures. Update
`docs/STATUS.md` with a link to the figure suite only after this final pass.

### Approved Figure 1 Amendment: Separate Genome And Protein Information

**Files:**
- Modify: `pipeline/scripts/build_scheme_a_figure_data.py`
- Modify: `pipeline/scripts/plot_scheme_a_figures.py`
- Modify: `pipeline/tests/test_scheme_a_figure_data.py`
- Modify: `pipeline/tests/test_scheme_a_figure_exports.py`

**Interface:** `figure1_genome_coverage.tsv` contains `stage`, `value`, `unit`,
and `note` rows for all GTDB representatives and the four-core-family tier1
genome union. Figure 1 uses it only for the genome-level panel; its existing
funnel remains gene/protein-level.

- [x] **Step 1: Approve Figure 1 genome/protein separation**
- [x] **Step 2: Add failing source-data and renderer tests**
- [x] **Step 3: Generate the genome-level source data and revise the renderer**
- [x] **Step 4: Regenerate Figure 1 and run the complete QA suite**

## Plan Self-Review

- Spec coverage: Tasks 1-3 cover all six figures, compact source data, provenance,
  editable exports, numerical invariants, and visual QA.
- Scope: workflow/funnel, seed provenance, family scale, phylum distribution,
  SignalP distribution, and neighborhood context are separate panels but share one
  data and rendering architecture; no phylogeny or new biological analysis is added.
- Ambiguity resolved: Figure 5 uses current snapshot tables, has no visible repair
  callout, and labels its support metric as genomic proximity only.
- Commit policy: deliberately omitted because the user has not authorized commits.
