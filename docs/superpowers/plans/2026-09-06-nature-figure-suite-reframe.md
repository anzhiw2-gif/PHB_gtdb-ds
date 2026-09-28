# Nature figure suite reframe (2026-09-06)

## Scientific contract

The main figure suite describes the bacterial/GTDB candidate census as two primary
families only: `ePhaZ` and `iPhaZ`. All labels are computational candidate evidence;
none are phenotype confirmations. The archaeal `PhaZh1-like` audit is an independent
supplementary branch and is never counted in the main-family totals.

## Evidence boundaries

- `ePhaZ` SignalP6 panels use only the 38,692 curated-core tier1 ePhaZ proteins.
- ePhaZ `OTHER` competition review uses only the 16,836 curated-core ePhaZ proteins
  predicted as `OTHER`.
- Main family/phylum/species distributions use the bound `tier1_genome_family.tsv`
  taxonomy table and filter to `ePhaZ`/`iPhaZ`; they do not borrow archaeal taxonomy.
- Domain/motif/neighborhood and positive/negative controls are shown only in the
  archaeal supplementary branch (204 domain-annotated records; 197 high + 7 review).

## Figure map

1. Workflow and main-data funnel; archaeal audit shown only as a detached side branch.
2. Two-family hierarchy and evidence layers at protein and genome level.
3. ePhaZ SignalP6 scope and signal classes, with iPhaZ explicitly out of scope.
4. ePhaZ/iPhaZ tier1 genome distribution by GTDB phylum.
5. Species/genus representation and ePhaZ/iPhaZ genome co-occurrence.
6. ePhaZ `OTHER` competition audit and conservative PHB candidate boundary.
7. Supplementary archaeal PhaZh1-like audit strata and domain/motif evidence.
8. Supplementary archaeal neighborhood markers plus positive/negative controls.

The renderer uses Python/Matplotlib only, writes source-data TSVs, and exports SVG,
PDF, TIFF and PNG under a new dated run without overwriting historical runs.
