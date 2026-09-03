# Public repository guide

This repository contains source code, curated seed records, HMM profiles, compact
reports, and reproducibility contracts. Raw GTDB proteomes, bulk HMMER output,
server logs, `runs/`, and `deploy/` directories are machine-local and are not
part of the public release surface.

## Configure a private execution environment

Set these variables in the shell or job scheduler rather than committing them:

```bash
export PHB_REPO_ROOT=/path/to/PHB_gtdb-ds
export PHB_SOURCE_ROOT=/path/to/PHB_gtdb-ds
export PHB_PYTHON=/path/to/python
export PHB_GTDB_ROOT=/path/to/gtdb_genomes_reps_r232/database
export PHB_GTDB_TAXONOMY=/path/to/bac120_taxonomy_r232.tsv
export PHB_GTDB_METADATA=/path/to/bac120_metadata_r232.tsv.gz
export PHB_GTDB_TREE=/path/to/bac120_r232.tree
export PHB_CONDA_SH=/path/to/conda.sh
```

The synchronization helper also requires explicit `SYNC_SERVER` and
`SYNC_REMOTE` values. No account name, host address, hostname, or personal
directory is assumed by the repository.

## Interpretation boundary

HMM, domain, SignalP, neighborhood, and tree results identify candidate
homology or functional potential. They do not by themselves establish an
experimentally verified PHB or MCL-PHA degradation phenotype. Family-specific
thresholds and the current run authority are recorded in `docs/STATUS.md` and
`docs/CURRENT_STATUS_20260902.md`.

Before publishing a new run, record the input paths, versions, sizes, and
SHA-256 values in its `input_contract.json`; use placeholders or environment
variables in public documentation.
