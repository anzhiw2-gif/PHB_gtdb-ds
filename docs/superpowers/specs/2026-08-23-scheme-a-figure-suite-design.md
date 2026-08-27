# Scheme A Figure Suite Design

## Purpose

Create an English, double-column figure suite that communicates the current
GTDB R232 PHB/PHA candidate-homolog workflow and distribution results. The
figures show sequence, taxonomy, and genomic-context evidence only. They do
not claim enzyme activity, PHB/PHA degradation phenotype, complete pathways,
or co-transcription.

## Shared Contract

- Backend: Python only (matplotlib).
- Primary exports: editable SVG and PDF, plus 600 dpi TIFF and PNG preview.
- Typography: English, sans-serif, readable at approximately 183 mm width.
- Scope: current Scheme A outputs. No phylogenetic-tree figure.
- Evidence language: `candidate homologs`, `putative`, `signal-peptide prediction`,
  and `genomic neighborhood context`.
- Colors: one restrained palette across the suite; family color assignments are
  stable across all figures.
- Source data: each quantitative panel gets a TSV/CSV derived directly from the
  recorded project tables and an adjacent machine-readable provenance file.

## Figure 0: Curated Seed Library And Evidence Provenance

Core claim: the reference library is family-resolved and records review status,
training/validation assignment, and evidence provenance.

- Panel a: seed counts per family, split into reviewed and unreviewed entries.
- Panel b: train versus validation allocation.
- Panel c: evidence-source composition (PMID, DOI, or both) per family.
- Inputs: seed manifest and seed statistics.

## Figure 1: Workflow And Screening Funnel

Core claim: a fail-closed workflow links a genome-level census to a separate
gene/protein-level filtering record; the two count types are not conflated.

- Left panel: GTDB R232 through protein prediction, nine-family HMM screening,
  family-aware coverage/arbitration, sequence validation, and tier rescoring.
- Upper-right panel: 199,923 GTDB representative genomes and the 44,814-genome
  union with at least one tier1 candidate from the four core families.
- Lower-right panel: hit rows, unique protein candidates, validated candidates,
  and strict tier1 sequences for the four core families.
- Genome and gene/protein counts are separate panels with explicit units;
  broad patatin is not used as the core-funnel endpoint.

## Figure 2: Core Candidate Family Scale

Core claim: the four core candidate families differ substantially in sequence
and genome-level representation, while their genome union spans 44,814 GTDB
representatives.

- Panel a: tier1 sequence counts by family.
- Panel b: tier1 genome counts by family.
- Panel c: four-family genome union relative to all 199,923 GTDB representatives.
- Inputs: tier1 genome-family table and tier1 FASTA counts.

## Figure 3: Phylum-Level Distribution

Core claim: core candidate homologs occur across multiple GTDB phyla, with
family-specific concentration patterns.

- Panel a: total core-candidate genomes in leading phyla.
- Panel b: phylum-by-family matrix on a log scale.
- Broad patatin is excluded from the core matrix.
- Inputs: tier1 phylum distribution and tier1 genome-family table.

## Figure 4: ePhaZ Evidence-Stratified Distribution

Core claim: ePhaZ candidate homologs contain distinct signal-peptide prediction
classes and these classes have different phylum-level representation.

- Panel a: SignalP class composition of tier1 ePhaZ candidates.
- Panel b: leading phyla among candidates with a predicted signal peptide.
- Inputs: SignalP summary and phylum tables.

## Figure 5: Genomic Neighborhood Context

Core claim: current candidate loci show family-specific local proximity to the
six available PHB/PHA-associated marker families.

- Panel a: candidate-family by marker-family neighborhood support rate within
  plus/minus 10 kb.
- Panel b: ArchPhaZ_patatin loci and genomes with at least one nearby marker.
- Inputs: current cluster context, locus audit, and summary tables from the
  Scheme A run snapshot.
- The figure presents current progress without a visible remediation callout.
  Provenance files retain the exact run and source-table identifiers.

## Quality Checks

- Every figure includes a panel label and readable text at final width.
- Quantitative source data are reproduced from the raw TSV/FASTA inputs and
  cross-checked against displayed labels.
- SVG text remains editable; PDF uses embedded TrueType fonts.
- No statistical significance test or uncertainty interval is shown because
  these are census-style database counts, not replicate experiments.
- Final visual QA checks clipping, overlaps, color contrast, text editability,
  and agreement between labels and source data.
