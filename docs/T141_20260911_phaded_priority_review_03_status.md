# T141 PhaDED priority review integrity correction status

Date: 2026-09-11
Run: `runs/20260911_phaded_priority_review_03/`

## Completed

- Corrected sequence-integrity classification in `annotate_phaded_features.py`:
  - one terminal `*` after standard amino acids -> `terminal_stop_only`;
  - internal `*`, multiple terminal stops, or another non-IUPAC residue -> `invalid_internal_character`.
- Updated `prioritize_phaded_candidates.py` so `terminal_stop_only` receives the same score as `valid`, while internal invalid residues retain a `-4` penalty.
- Rebuilt the 109,087-record corrected feature annotation:
  - `terminal_stop_only`: 105,748
  - `valid`: 1,388
  - `invalid_internal_character`: 351
  - `possible_N_truncation`: 1,600
- Rebuilt the bounded review queue and InterPro priority selection from the frozen priority source.
- Old and new accession order and sets are identical:
  - review queue: 71,671 records;
  - InterPro priority candidates: 244 records;
  - no accession additions, removals, or reordering.
- No HMMER, InterProScan, GTDB rescan, formal registry edit, or server execution was performed.

An early intermediate directory, `results/interpro_priority_v3/`, used the pre-Pfam assignment table and produced 71,618 queue records / 161 InterPro candidates. It is retained as audit residue but is superseded and must not be used. The authoritative rebuilt selection is `results/interpro_priority_v4/`.

The first corrected feature annotation is also retained as superseded residue. Version 2 distinguishes terminal-stop sequences that do not begin with `M` as `possible_N_truncation`; this avoids treating likely N-terminal truncations as fully intact proteins.

## Provenance correction

The historical InterPro run retains its original contract unchanged. A provenance amendment records that the declared hash for `interpro_summary_manifest.json` was stale:

- declared: `DEEC616D0AB136E909062A54F0A7F584E543DEF4B5D7EA531464EA0D90129BD4`
- actual: `816D43F106557845354F18BBE6F7E2FE1635C7256E42C37EE7F3464F0194BBD3`

See `runs/20260911_phaded_priority_review_03/results/interpro_summary_manifest_hash_amendment.json`.

## Interpretation boundary

The correction changes quality metadata and provenance only. It does not change the PhaDED assignment, Pfam architecture state, InterProScan result, or candidate-only phenotype boundary.

## Next action

Proceed to sequence-specific review in this order:

1. Four explicit extracellular PHA-related candidates.
2. Twelve PhaC-versus-PhaZ competition candidates.
3. Eight motif/domain or assignment-gap exceptions.
4. Eighteen `ePhaZ_competition_review` records.

Only after those reviews should structure prediction, phylogenetic placement, or experimental shortlisting be performed. The four-seed ePhaZ repair remains a separate calibration gate before any replacement model or formal rescan.
