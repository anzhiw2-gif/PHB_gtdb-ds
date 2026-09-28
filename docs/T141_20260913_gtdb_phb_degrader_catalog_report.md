# GTDB PHB Degrader Evidence Catalog

Date: 2026-09-13  
Run: `runs/20260913_gtdb_phb_degrader_literature_reconciliation_01/`

## Main conclusion

Using the current Tier1 GTDB proxy table and accession-bound public records,
three proteins have direct PHB degradation evidence in the literature:

- `O05527`, *Delftia acidovorans* YM1609, PMID `9406404`, DOI
  `10.1128/aem.63.12.4844-4852.1997`.
- `P12625`, *Alcaligenes faecalis* T1 / current *Ralstonia pickettii*
  record, PMID `2644188`, DOI `10.1128/jb.171.1.184-189.1989`.
- `Q0K9H3`, *Cupriavidus necator* H16, genome `NC_008313.1`, PMID
  `12775684`, DOI `10.1128/jb.185.12.3485-3490.2003`.

None of these three literature records was found as an exact accession in the
current Tier1 GTDB proxy table. Therefore the GTDB genome catalog has zero
Tier A genomes in this run. This is an unresolved bridge state, not evidence
that the organisms lack PHB degradation.

## GTDB genome-level catalog

The genome evidence ledger contains 62,666 unique GTDB genomes from 109,087
candidate proteins:

| Tier | Count | Interpretation |
|---|---:|---|
| A | 0 | exact literature-supported GTDB genome bridge |
| B | 50,758 | genome-supported computational candidate; no exact PHB phenotype record |
| C | 11,903 | sequence-only or unresolved candidate evidence |
| H | 5 | architecture/source/structure hold |

The underlying unit is one unique GTDB genome. Protein count, genome x
PhaDED-superfamily count (`54,673`), and genome x PhaDED-family count
(`58,290`) are separate denominators and must not be combined.

## High-information candidates

The 19 manually reviewed records remain candidate-only:

- 3 dPHASCL2 architecture conflicts remain `hold_architecture_conflict`.
- 2 source/gene-model or low-structure records remain `hold_gene_model_or_structure`.
- 14 records remain bounded dPHASCL1-like candidates; their review order is in
  `results/experimental_priority_14.tsv` from the manual-review run.

No Foldseek score, pLDDT, MRCA support, InterPro hit, or PhaDED family label is
treated as direct PHB degradation evidence.

## Evidence gaps

- The current bridge input is a Tier1 genome-family proxy, not a complete GTDB
  metadata/protein accession table. Exact literature-to-GTDB matching therefore
  remains incomplete.
- The three direct literature records require assembly/protein accession
  resolution against a full GTDB release table before any GTDB Tier A claim.
- Two annotation-only UniProt records are retained as annotation-only and do
  not support a PHB phenotype claim.
- Missing structures, missing Foldseek pairs, unresolved taxonomy, and
  incomplete gene models remain pending/hold states, not biological negatives.

## Methods and provenance

The analysis reuses frozen PhaDED assignment, feature, taxonomy, and manual
review outputs. Literature sources were queried from UniProtKB REST and Europe
PMC REST API v6.9 on 2026-09-13; accession-level hashes and citations are in
`inputs/phb_literature_ledger.tsv`. The GTDB proxy is recorded as `GTDB R232;
Tier1 summary proxy` with SHA-256 in `results/gtdb_bridge_manifest.json`.

No wet-lab work was planned or executed. No HMM, formal registry, historical
run, server deploy, or GTDB-wide rescan was modified or started.
