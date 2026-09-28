# PhaDED intermediate-structure burden review

Date: 2026-09-12  
Run: `runs/20260912_phaded_intermediate_structure_review_01/`  
Status: `completed_candidate_only`

## Scope and ranking

All nine `intermediate_structure_manual_review` dPHASCL1-like candidates were
audited. Ranking uses the predeclared rule in the input contract: lowest
outside-domain low-pLDDT fraction, longest outside-domain low run, then lowest
within-domain pLDDT burden and highest domain mean. The first four are the
priority manual-review subset; the other five remain in the complete audit
table.

## Priority four

| Rank | Candidate | PFAM interval | Domain pLDDT | Outside-domain pLDDT | Longest outside low run | PAE domain internal | PAE N-term/domain |
|---:|---|---:|---:|---:|---:|---:|---:|
| 1 | `GCA_963791955.1|CAWBUO010000005.1_63` | 132-294 | 83.765 | 68.660 | 50 | 6.301 | 11.620 |
| 2 | `GCA_023389175.1|JAEWEW010000389.1_3` | 157-345 | 92.099 | 64.959 | 108 | 4.058 | 27.447 |
| 3 | `GCA_051861205.1|JBHYDC010000002.1_261` | 162-356 | 90.744 | 66.084 | 105 | 4.749 | 28.065 |
| 4 | `GCA_031386065.1|JANNBZ010000062.1_14` | 107-314 | 89.017 | 65.039 | 90 | 6.050 | 28.047 |

## Interpretation

The PFAM cores are generally more reliable than the flanking regions. Rank 1
has the weakest domain confidence and should be inspected first for a genuine
domain-boundary problem. Ranks 2 and 3 have strong cores but long low-confidence
regions outside the core, which may be linker/SBD or unresolved flexible
sequence; this is not evidence of absence. Rank 4 combines a relatively weak
domain and a 90-residue N-terminal low-confidence run.

The PAE values support compact internal cores (approximately 4-6 Å mean
internal PAE) but higher N-terminal/domain separation (approximately 11.6-28.1
Å), so termini should not be interpreted as rigidly attached to the catalytic
core. The candidate sequences remain `candidate-only`; architecture and
structure do not establish degradation activity.

## Next action

Inspect these four PDBs and their Pfam/InterPro coordinates in a molecular
viewer, documenting whether low-confidence regions correspond to signal
peptides, linkers, substrate-binding domains, or unresolved insertions. Then
review the remaining five intermediate structures. Foldseek/TM-align remains
deferred until reference structures are bound and the comparison contract is
approved.

