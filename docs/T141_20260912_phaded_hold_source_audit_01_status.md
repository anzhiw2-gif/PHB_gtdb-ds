# T141 hold-candidate source audit

Date: 2026-09-12  
Run: `runs/20260912_phaded_hold_source_audit_01/`  
Status: `completed_candidate_only_pending_external_source`

## Finding

The local project contains the two translated proteins and their candidate
mapping, but not the source contigs or GFF3/GBFF annotation needed to verify
ORF boundaries.  The available neighborhood audit only records both loci as
`analyzed`; it does not contain nucleotide coordinates or strand information.

### `GCA_030946875.1|JALYXY010000192.1_6`

- Cleaned protein length: 209 aa; original translated record: 210 aa with a
  terminal `*`.
- Pfam spans begin at residue 6 and structure pLDDT is high (88.04), but the
  corrected feature annotation marks `possible_N_truncation`.
- Keep as `hold_gene_model_or_structure`.  An upstream CDS/contig check is
  required before structural comparison or any shortlist decision.

### `GCF_006442545.2|NZ_CP083952.1_5059`

- Cleaned protein length: 202 aa; original translated record: 203 aa with a
  terminal `*`.
- PF10503 covers residues 2-100, while the best structure has pLDDT 64.09,
  49.0% residues below 70, and an 85-residue longest low-confidence run.
- The same genome contains a separate 384-aa ePhaZ-like record,
  `GCF_006442545.2|NZ_CP083952.1_5067`.  This is evidence for a competing
  neighboring/alternative hit, not permission to merge the two proteins.
- Keep as `hold_low_structure_gene_model` until the source annotation is
  available.

## Required external evidence

For each candidate, obtain the exact matching contig FASTA and GFF3/GBFF,
record accession/version, retrieval date, file sizes and SHA-256, then verify:

1. CDS start, end, strand and translation against the protein accession;
2. whether an upstream start codon or alternate CDS explains the possible
   truncation;
3. whether `5059` and `5067` are separate ORFs, overlapping ORFs, or annotation
   artifacts;
4. whether the predicted protein sequence changes before repeating structure
   prediction.

Until those files are bound, missing source annotation remains a provenance gap,
not a biological negative.  No curated model, positive training set, or formal
GTDB scan is changed.

