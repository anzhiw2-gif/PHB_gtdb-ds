# T141 GTDB PHB degrader literature reconciliation amendment

Date: 2026-09-13  
Run: `runs/20260913_gtdb_phb_degrader_literature_reconciliation_02/`  
Status: `completed_full_metadata_accession_reconciliation`

## Amendment result

The prior run used a Tier1 summary proxy and therefore could not resolve exact
assembly accessions. This amendment used the server-bound GTDB R232 bac120
metadata file (`${PHB_REMOTE_ROOT}/GTDB/metadata/bac120_metadata_r232.tsv.gz`)
through a read-only SSH retrieval and retained a compact, hashed subset.

`Q0K9H3` (*Cupriavidus necator* H16) now has an exact assembly bridge from
`NC_008313.1` to GTDB representative `RS_GCF_000009285.1`. Because the ledger
row is `experimental_positive`, this is one Tier A direct-literature-supported
GTDB genome record.

`O05527` (*Delftia acidovorans* YM1609) has only species-level taxonomy context
in R232 and remains sequence-only context. `P12625` remains unresolved because
the cited T1 strain was not exactly identified in the metadata subset. The two
annotation-only rows remain unresolved and are not phenotype evidence.

## Counts

| bridge class | rows |
|---|---:|
| exact accession / Tier A | 1 |
| taxonomy-only context | 129 |
| unresolved hold | 3 |
| total bridge rows | 133 |

The 129 taxonomy-only rows are one-to-many species-context matches and must not
be interpreted as 129 experimentally validated degraders. The exact bridge is
the only Tier A result in this amendment.

## Provenance and boundaries

Remote resource sizes and SHA-256 values, local subset hashes, commands, and
the mapping policy are recorded in
`runs/20260913_gtdb_phb_degrader_literature_reconciliation_02/input_contract.json`
and `results/gtdb_bridge_manifest.json`. No wet-lab work, formal registry edit,
GTDB-wide rescan, server script execution, deploy, deletion, commit, or push was
performed.

The Tier A label is limited to the exact accession-bound literature record; all
PhaDED/domain/structure/tree evidence remains computational support and does
not independently prove a PHB degradation phenotype.
