# GTDB PHB degrader catalog: full-metadata amendment

Date: 2026-09-13  
Run: `runs/20260913_gtdb_phb_degrader_literature_reconciliation_02/`

This amendment supersedes the prior proxy-only bridge count for literature
records. The computational candidate universe and its B/C/H denominators are
unchanged; only the literature-to-GTDB accession bridge was improved.

## Tier A: exact literature-supported GTDB record

| GTDB representative | literature accession | organism/strain | evidence | citation |
|---|---|---|---|---|
| `RS_GCF_000009285.1` | `NC_008313.1`; `Q0K9H3` | *Cupriavidus necator* H16 | direct PHB/PhaZ2 hydrolysis experiment | PMID `12775684`; DOI `10.1128/jb.185.12.3485-3490.2003` |

This is one exact assembly-bound Tier A record. It is not a claim that every
*C. necator* genome or every homolog is experimentally validated.

## Tier B/C/H computational catalog

The frozen candidate universe remains 109,087 proteins across 62,666 genomes:

| unit | count | interpretation |
|---|---:|---|
| genome | 62,666 | unique GTDB genome denominator |
| genome x superfamily | 54,673 | separate denominator |
| genome x family | 58,290 | separate denominator |
| Tier B | 50,758 | computationally supported candidate |
| Tier C | 11,903 | sequence-only/unresolved candidate |
| Tier H | 5 | architecture/source/structure hold |

These B/C/H values remain candidate potential, not phenotype calls. Foldseek,
Pfam/InterPro, SignalP, structures, neighborhoods, and tree placement are
supporting evidence only.

## Bridge outcomes for the five literature records

| record | result | boundary |
|---|---|---|
| `phb-uniprot-q0k9h3` | exact `NC_008313.1` -> `RS_GCF_000009285.1` | Tier A |
| `phb-uniprot-o05527` | 129 species-context rows for *Delftia acidovorans* | taxonomy-only; not Tier A |
| `phb-uniprot-p12625` | no exact T1 assembly/strain bridge | unresolved hold |
| `annotation-uniprot-a0a0f6ygp7` | no direct assay and no exact bridge | annotation-only hold |
| `annotation-uniprot-a0a1c4yai9` | conflicting annotation and no exact bridge | annotation-only hold |

## Missing-evidence table

| gap | status | consequence |
|---|---|---|
| Exact `Delftia` YM1609 strain assembly | unresolved | species context cannot transfer phenotype |
| Exact *Alcaligenes* T1 assembly | unresolved | taxonomy/renaming discrepancy remains |
| Protein accession-to-GTDB gene mapping for Q0K9H3 | pending | assembly bridge is exact; protein locus bridge is not asserted |
| Full per-genome Pfam/InterPro/SignalP/Foldseek/tree join | retained in prior candidate runs | supports potential only; not needed to establish this exact Tier A bridge |
| Wet-lab confirmation for GTDB records | `not_planned_by_scope` | absence is not a negative result |

## Provenance

The remote GTDB R232 metadata source, taxonomy/tree hashes, local compact
metadata subset, bridge output hash, and mapping policy are recorded in the
amendment `input_contract.json` and `results/gtdb_bridge_manifest.json`.
Retrieval used read-only SSH from `<SERVER_USER>@<SERVER_HOST>`; no server script was
executed, and no historical result, registry, HMM, deploy, or run was changed.
