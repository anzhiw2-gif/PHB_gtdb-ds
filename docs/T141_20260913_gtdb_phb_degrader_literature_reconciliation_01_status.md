# T141 GTDB PHB-degrader literature reconciliation

Date: 2026-09-13  
Run: `runs/20260913_gtdb_phb_degrader_literature_reconciliation_01/`  
Status: `completed_candidate_only_genome_catalog`

## Completed layers

- Literature ledger: 5 accession-level rows, including 3 direct PHB
  experimental-positive records and 2 annotation-only records.
- GTDB bridge: 5 rows retained; all are unresolved against the current Tier1
  GTDB proxy, so no Tier A genome was assigned.
- Genome evidence ledger: 62,666 unique genomes from 109,087 proteins,
  classified as B=50,758, C=11,903, H=5, A=0.
- High-information adjudication: 19 unique candidates, with 3 architecture
  holds, 2 source/structure holds, and 14 retained candidate-only records.
- No-wet-lab scope: no experiment, assay result, or phenotype promotion was
  performed or implied.

## Interpretation

The catalog answers two different questions separately:

1. Which organisms/proteins have direct PHB degradation evidence in public
   literature? Three protein-level records are currently supported.
2. Which GTDB genomes contain computational PhaDED/ePhaZ/iPhaZ candidates?
   50,758 are Tier B computationally supported candidates in this proxy-based
   genome summary, while 11,903 are Tier C and 5 remain Tier H holds.

The current data do not support a claim that any specific GTDB Tier B genome
is an experimentally verified PHB-degrading bacterium. Exact literature-to-GTDB
assembly/protein bridging remains the next information-retrieval limitation.

## Provenance and boundaries

Source database versions, retrieval dates, accession identifiers, input/output
hashes, and the unchanged registry hash are bound in
`runs/20260913_gtdb_phb_degrader_literature_reconciliation_01/input_contract.json`.
The formal registry, historical runs, HMMs, server state, and full GTDB scan
were not modified or started.
