# T141 ePhaZ MCL-PHA cross-genus candidate 07

## Scope

This is a bounded candidate run only. It adds three experimentally supported,
sequence-complete MCL-PHA positives from *Bdellovibrio*, *Burkholderia* and
*Streptomyces*. No formal registry, run 13, or genome-wide result was changed.

## Evidence

| Candidate | Sequence binding | Experimental evidence | Completeness |
|---|---|---|---|
| `Q6MH49` | UniProtKB Q6MH49; Bd3709; 271 aa | Recombinant extracellular mcl-PHA depolymerase; PHB inactive; PMID 22706067 | complete protein record |
| `WHU94860.1` | CP120731.1 CDS 314007..315101(+); P4G95_16805; 364 aa | Tn-seq, RhaCAST and GC-FID gene-to-function evidence for extracellular MCL-PHA degradation; PMID 41930961 | complete annotated CDS |
| `AZSS01000334.1:12616-13485(-)` | WGS contig AZSS01000334.1, strand-aware coordinate | Cloned/purified PhaZSex2 with MCL-PHA/PHACOS degradation; PMID 26156240 | 870 nt = 290 codons, one terminal stop, 289-aa protein |

Raw records, source URLs, coordinate translation and SHA-256 values are under
`runs/20260902_ephaz_mcl_crossgenus_candidate_07/`.

## Results

The six-sequence MCL candidate profile recovers all three new positives in the
same-profile probe. The independent clean LOO run is
`results/loo_calibration_v4/loo_summary.tsv`:

- `Q51718`, `Q6UFW4`, `Q6MH49`: held-out hit at `E <= 1e-5`.
- `WHU94860.1`, `AZSS01000334.1:12616-13485(-)`: held-out miss at every tested threshold.
- Three intracellular non-ePhaZ negatives: zero hits.

## Decision

Do **not** request formal scan authorization. The single MCL profile does not
generalize across the divergent lipase-associated and Streptomyces branches.
The next candidate-only step is to collect additional accession-level,
experimentally confirmed examples for those branches and evaluate separate
subfamily profiles with independent LOO and cross-family challenge panels.

HMM hits remain evidence of sequence-family or functional potential, not proof
of a validated PHB degradation phenotype.

## Provenance note

An earlier two-source `scp` synchronization left a copy at the server project
root named `20260902_ephaz_mcl_crossgenus_candidate_07`. It was not used for
execution and has been retained as residue; the authoritative paths are the
dated `runs/` and `deploy/` directories above.
