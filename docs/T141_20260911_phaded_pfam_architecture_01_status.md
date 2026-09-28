# T141 PhaDED Pfam architecture evidence status

Date: 2026-09-11  
Run: `runs/20260911_phaded_pfam_architecture_01/`  
Scope: independent Pfam architecture evidence for the frozen PhaDED candidate mapping.

## Inputs and execution

- Candidate FASTA: 109,087 accessions, SHA-256 `64895bc74218eff0d798c32a78c9f655c071396eac9acf7b9d31316db0a8b4cc`.
- Pfam-A database: `${PHB_REMOTE_ROOT}/GTDB/pfam/Pfam-A.hmm`, size 2,246,909,846 bytes, SHA-256 `4b0da6399b97d2b23329de82190b2e5705bde1dea4cb2d028a52712d5966739e`.
- Tool: `${PHB_REMOTE_ROOT}/miniconda3/envs/phb_gtdb/bin/hmmscan`, HMMER 3.4 (Aug 2023), `--cut_ga --cpu 40`.
- Reference-first procedure: full Pfam-A scan of the PhaDED reference FASTA, followed by a targeted scan using 23 Pfam models present in the reference-derived fingerprint.

## Results

- Reference scan: 1,343 Pfam domain hits.
- Reference fingerprint: 6 superfamilies and 23 targeted Pfam models. Superfamilies without at least two shared reference hits remain without a fingerprint.
- Candidate targeted scan: 98,419 candidates with at least one targeted Pfam hit.

| Pfam architecture state | Count |
| --- | ---: |
| `partial` | 66,632 |
| `conflicting` | 83 |
| `pending` | 42,372 |

The `partial` state means that at least one Pfam accession is present in the reference fingerprint of the candidate's assigned PhaDED superfamily. It is not a catalytic-domain or activity proof. `conflicting` means that a Pfam accession associated with another reference superfamily was observed and the candidate requires review; the PhaDED assignment is retained. `pending` means that the targeted fingerprint did not resolve the architecture or that no superfamily fingerprint was available; it is not a biological negative.

## Candidate review boundary

The final table combines the original SignalP/motif feature fields with Pfam accessions and coordinates. It does not overwrite the original PhaDED mapping. The 83 conflicting records should be reviewed first, followed by the high-information ePhaZ competition and curated-core strata. InterProScan 5.78-109.0 remains available for a smaller priority set, but was not used as a silent substitute for this Pfam run.

All evidence remains candidate-only. Pfam, SignalP, motif, and PhaDED profile evidence cannot establish that an organism degrades PHB/PHA.

