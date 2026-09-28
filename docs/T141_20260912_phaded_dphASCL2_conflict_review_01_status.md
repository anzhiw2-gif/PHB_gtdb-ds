# PhaDED dPHASCL2 architecture-conflict structure review

Date: 2026-09-12  
Run: `runs/20260912_phaded_dphASCL2_conflict_review_01/`  
Status: `completed_candidate_only`

## Scope

This run audits the best-ranked ColabFold PDB for each of the three explicit
dPHASCL2-like candidates. It maps the observed PF00756 interval onto the PDB
residue numbering and reports pLDDT separately for the PF00756 interval, an
operational N-terminal window (residues 1-30), and the remainder of the chain.
The N-terminal window is not treated as a SignalP cleavage prediction.

## Findings

All three PDBs have 100% residue coverage. The PF00756 intervals are entirely
high-confidence and contain no residue with pLDDT <70:

| Candidate | PF00756 interval | Domain pLDDT | N-term 1-30 pLDDT | Low-confidence runs |
|---|---:|---:|---:|---|
| `GCA_963818845.1|CAWFTO010000003.1_420` | 64-207 | 97.247 | 25.642 | 1-58; 290-304; 394-396 |
| `GCF_000759345.1|NZ_JPZL01000003.1_113` | 39-141 | 97.005 | 44.029 | 1-28; 178-181; 352-354 |
| `GCF_050269715.1|NZ_JBNNIW010000004.1_132` | 34-141 | 98.071 | 59.284 | 1-21; 173-176; 178; 347 |

The low N-terminal segments are compatible with signal-peptide or flexible
leader behavior, but this audit does not infer cleavage sites. The predicted
fold is therefore structurally reliable for a compact PF00756-containing core;
it is not a failed prediction.

## Architecture decision

The three candidates remain `P0_architecture_conflict_hold`. Their high-quality
PF00756 cores do not supply the missing dPHASCL2 reference fingerprint
(`PF00041`, `PF00326`, `PF00734`, `PF10503`). InterPro's explicit PHA-related
records and the dPHASCL2 profile placement are compatible with a PHA-related
esterase, but they cannot distinguish a dPHASCL2 family assignment from a
generic/alternative esterase architecture when the Pfam evidence is
contradictory. High pLDDT, SignalP class, or phylogenetic proximity does not
override this conflict.

For the second and third candidates the PhaDED family gaps are 0.8 and 0.1,
respectively, which further supports retaining ambiguity. The first has a
larger gap but the same missing type-2 architecture.

## Next action

Do not promote these records, alter `ePhaZ_curated_core`, or repair the
four-seed model. The next bounded comparison may use Foldseek/TM-align against
explicit dPHASCL1 and dPHASCL2 reference structures, but any structural
similarity result remains candidate evidence and cannot by itself resolve the
Pfam architecture conflict or establish PHB/PHA degradation phenotype.

