#!/usr/bin/env python3
"""Curate the frozen 723-row PhaDED reference ledger with the v2 evidence columns.

Task F13 of ``docs/superpowers/plans/2026-09-28-phaded-evidence-model-redesign-followup.md``.

What this script is allowed to do
---------------------------------
It is **append-only**.  It reads

* the frozen reference ledger (723 rows, 26 columns),
* the 30-row primary-source amendment produced by
  ``runs/20260920_phaded_three_gaps_01`` (one row per ``experimental_positive``),
* the frozen reference FASTA (sequence integrity) and
* the 15-column evidence contract in ``pipeline/config``,

and writes a **new** ledger carrying the same rows plus the 15 v2 evidence columns.
Every original column keeps its original value byte for byte; the frozen ledger is
never opened for writing and its SHA-256 is re-checked before and after the run.

How a row is graded
-------------------
The grading rules come from the design spec (Task 4) and are applied strictly --
when the frozen evidence does not prove a grade, the row is downgraded:

``E3``
    purified-protein direct substrate conversion, product or kinetics evidence.
``E2``
    knockout/complementation, cellular polymer change or other strong in-vivo
    causal evidence, i.e. a cloned gene product with a demonstrated activity.
``E1``
    expression, correlative phenotype or another indirect experimental link.
``A``
    database annotation, sequence inference or an inherited historical label.

The amendment classifies every positive by ``evidence_type`` and this script maps
that documented classification onto the grade ladder:

===================  =====  ==================
``evidence_type``    grade  directness
===================  =====  ==================
``purified_activity``  E3   direct
``cloned_activity``    E2   direct
``sequence_only``       A   annotation
``same_gene_duplicate`` A   annotation
===================  =====  ==================

A ``cloned_activity`` row is never promoted to ``E3`` here: the frozen curation
inputs document a cloning/activity claim, not purified-protein kinetics.  The one
``purified_activity`` row (``DED_hfam_70_0016``, purified and crystallised enzyme)
is the only row whose evidence_type documents purification, so it is the only
``E3``.  No row is *promoted* above its documented evidence_type, and no grade is
invented: anything the frozen evidence does not fix is written as ``pending`` and
listed in ``pending_value_classes``.

How independence groups are derived
-----------------------------------
Two positives are *the same piece of evidence* when

* they name the same gene -- the documented ``seed_gi`` of the Knoll 2009 Table 1
  seed, cross-checked between the frozen ledger's ``notes`` column
  (``Knoll 2009 Table 1 seed GI <n>``) and the amendment's ``seed_gi`` column, or
* they come from the same primary study -- the amendment's ``primary_pmid`` /
  ``primary_doi`` for indexed papers and a ``CITATION:`` anchor for the documented
  non-indexed ones (Briese 1994, Kobayashi 1999, GenBank submissions).

The independence group is therefore the connected component of those two
documented relations (union-find over the rows in ledger order), and the group id
records which relation bound it: ``IG_GI_<seed_gi>`` when all members are one
gene, otherwise ``IG_<study>``.  ``functional_positive_count`` is then counted by
unique ``independence_group`` by the frozen validator, never by row count.

``discovery_training_eligible`` and ``sequence_integrity``
---------------------------------------------------------
``sequence_integrity`` is evaluated from the bound FASTA sequence with the rule
the project already uses in ``build_phaded_evidence_tiers.normalized_sequence_integrity``
and ``annotate_phaded_features``, mapped onto the evidence vocabulary:

* clean 20-residue alphabet, initiator Met (an optional single trailing ``*`` is
  stripped) -> ``complete``;
* clean alphabet but no initiator Met -> ``partial`` (possible N-terminal
  truncation or mature form; completeness cannot be established);
* any ambiguous/non-standard residue or an internal stop -> ``unresolved``.

``discovery_training_eligible`` is ``true`` only for a ``complete`` sequence with
a consistent family call; every other row is ``false`` (fail-closed).
``functional_calibration_eligible`` is ``true`` only for a direct ``E2``/``E3``
row with a documented independence group and a ``complete`` sequence.

Acceptance gate
---------------
The script re-runs the frozen validator itself in ``strict`` evidence-column mode
over the ledger it just wrote and refuses to finish if the gate does not pass.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import re
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


# --- frozen contracts -----------------------------------------------------------

#: The 26 frozen ledger columns, in file order.  The curated ledger keeps them in
#: this order and with these exact values.
FROZEN_LEDGER_COLUMNS = (
    "reference_id", "accession", "sequence_sha256", "sequence_source", "source_database",
    "source_database_version", "retrieval_date", "organism", "taxonomy_id",
    "phaded_superfamily", "phaded_family_id", "reported_localization", "substrate_class",
    "substrate_detail", "experimental_assay", "experimental_result", "evidence_status",
    "positive_negative_control", "catalytic_residues", "lipase_box_or_ahsmg",
    "oxyanion_hole_evidence", "domain_architecture", "primary_doi", "pmid", "pmcid", "notes",
)

#: The 17 columns of the 30-row primary-source amendment.
PROVENANCE_COLUMNS = (
    "reference_id", "accession", "organism", "phaded_superfamily", "phaded_family_id",
    "reported_localization", "substrate_class", "seed_gi", "ledger_pmid", "ledger_doi",
    "primary_pmid", "primary_doi", "primary_title", "evidence_type", "confidence",
    "source_url", "note",
)

EVIDENCE_COLUMNS_CONFIG = (
    Path(__file__).resolve().parents[1] / "config" / "phaded_reference_evidence_columns.tsv"
)
VALIDATOR_PATH = Path(__file__).resolve().with_name("validate_phaded_reference_ledger.py")

EVIDENCE_STATUS_ANNOTATION = "annotation_only"
EVIDENCE_STATUS_POSITIVE = "experimental_positive"

#: "nothing here to assess" declared value for the free-text evidence columns, used
#: exactly as the validator's own grade-defaults table uses it.
NOT_ASSESSED = "not_assessed"
#: The AGENTS.md sentinel for "no value exists / not curated"; never inferred.
PENDING = "pending"

AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")
_SEED_GI_PATTERN = re.compile(r"seed GI (\d+)")
_VERSION_SUFFIX = re.compile(r"^(?P<base>[A-Za-z0-9]+)\.(?P<version>\d+)$")
GENE_ID_PREFIX = "gi:"

#: Sequences whose integrity cannot be established: excluded from discovery training.
INTEGRITY_VOCABULARY = ("complete", "partial", "truncated", "unresolved")

EVIDENCE_COLUMNS_FALLBACK = (
    "accession_version", "gene_id", "study_id", "experiment_unit_id", "independence_group",
    "experimental_evidence_grade", "assay_directness", "experimental_substrate_class",
    "experimental_system", "negative_control_description", "catalytic_mutant_evidence",
    "localization_evidence", "sequence_integrity", "discovery_training_eligible",
    "functional_calibration_eligible",
)


# --- the auditable curation table ------------------------------------------------
#
# One entry per `experimental_positive` row of the frozen ledger.  Everything here
# is copied from the frozen ledger's `notes` column or from the 30-row amendment
# (evidence_type / confidence / primary_pmid / primary_doi / primary_title / note);
# nothing is inferred from outside the frozen evidence.  ``study_key`` is the
# primary source of the row: an indexed paper is anchored on ``PMID:<pmid>`` and a
# documented non-indexed source on ``CITATION:<...>``; the value becomes the row's
# ``study_id`` and, together with ``seed_gi``, the binding used for the
# independence groups.
#
# Per-evidence-type defaults keep the table small; only the documented row-specific
# facts are spelled out.  ``grade_basis`` is the one-line justification that the
# curation report prints for every non-E3 row.

EVIDENCE_TYPE_RULES: dict[str, dict[str, str]] = {
    "purified_activity": {
        "grade": "E3",
        "assay_directness": "direct",
        "experimental_system": "purified_enzyme_in_vitro",
        "reason": (
            "amendment evidence_type=purified_activity: purified-protein direct substrate "
            "conversion/product evidence -> E3"
        ),
    },
    "cloned_activity": {
        "grade": "E2",
        "assay_directness": "direct",
        "experimental_system": "cloned_gene_product_activity_assay",
        "reason": (
            "amendment evidence_type=cloned_activity: a cloned gene product with a demonstrated "
            "activity is direct in-vivo/cellular causal evidence, but no purified-protein "
            "conversion or kinetics are documented in the frozen evidence -> E2, not E3"
        ),
    },
    "sequence_only": {
        "grade": "A",
        "assay_directness": "annotation",
        "experimental_system": "sequence_record_only_no_assay",
        "reason": (
            "amendment evidence_type=sequence_only, confidence=none: the record carries a "
            "database/sequence annotation only, with no experimental assay -> A"
        ),
    },
    "same_gene_duplicate": {
        "grade": "A",
        "assay_directness": "annotation",
        "experimental_system": "same_gene_duplicate_record_no_independent_assay",
        "reason": (
            "amendment evidence_type=same_gene_duplicate: the record is a second database entry "
            "of an already counted gene and adds no independent assay -> A"
        ),
    },
}

DEFAULT_NEGATIVE_CONTROL = NOT_ASSESSED
DEFAULT_CATALYTIC_MUTANT = NOT_ASSESSED
#: An inherited DED/GenBank localization label is a report, never experimental
#: localization evidence (AGENTS.md: transport-signal prediction != localization truth).
DEFAULT_LOCALIZATION_EVIDENCE = "inherited_database_family_label"

POSITIVE_CURATION: dict[str, dict[str, str]] = {
    "DED_hfam_65_0001": {
        "study_key": "PMID:11114905",
        "study_basis": (
            "amendment primary_pmid 11114905 = Saegusa 2001, 'Cloning of an intracellular "
            "poly[D(-)-3-hydroxybutyrate] depolymerase gene from Ralstonia eutropha H16 and "
            "characterization of the gene product' (the gene's characterization paper)"
        ),
        "experimental_system": "cloned_gene_product_characterization_native_host",
        "catalytic_mutant_evidence": (
            "reported_in_a_separate_study_PMID:16233560_not_the_primary_characterization"
        ),
        "grade_basis": (
            "evidence_type=cloned_activity, confidence=high: the gene product was characterized "
            "after cloning, but no purified-protein conversion/kinetics are documented in the "
            "frozen evidence -> E2; the C183 catalytic-mutant result belongs to a different "
            "paper (PMID 16233560) and is not counted towards this row's grade"
        ),
        "note_excerpt": (
            "seed_gi 3641686 = BAA33394.1 (same gene as CAJ92291.1, PhaZ1/PhaZa1). Gene product "
            "characterized. The C183 catalytic-residue mutagenesis is a later paper (PMID "
            "16233560), not the primary characterization."
        ),
    },
    "DED_hfam_2_0013": {
        "study_key": "PMID:16936025",
        "study_basis": (
            "amendment primary_pmid 16936025 = Tseng 2006, 'Identification and characterization "
            "of the Bacillus thuringiensis phaZ gene, encoding new intracellular "
            "poly-3-hydroxybutyrate depolymerase'"
        ),
        "negative_control_description": (
            "S102A catalytic mutant abolishes activity (same primary study)"
        ),
        "catalytic_mutant_evidence": "S102A_abolishes_activity_primary_study",
        "grade_basis": (
            "evidence_type=cloned_activity, confidence=high with a documented catalytic-mutant "
            "(S102A) result, but no purified-protein conversion/kinetics in the frozen evidence "
            "-> E2, not E3; the phaZ identity comes from Tseng 2006, not from the dead RefSeq "
            "record ZP_00743157.1"
        ),
        "note_excerpt": (
            "seed_gi 75763431 = ZP_00743157.1 (RefSeq, status=dead, annotated '3-oxoadipate "
            "enol-lactonase'). The phaZ identity comes from Tseng 2006, NOT from the GenBank "
            "record (which has no PUBMED link). S102A abolishes activity."
        ),
    },
    "DED_hfam_4_0001": {
        "study_key": "PMID:1989978",
        "study_basis": (
            "amendment primary_pmid 1989978 = Huisman 1991, 'Metabolism of poly(3-hydroxyalkanoates) "
            "(PHAs) by Pseudomonas oleovorans ... function of the encoded proteins in the synthesis "
            "and degradation of PHA'"
        ),
        "experimental_system": "in_vivo_complementation_of_a_PHA_degradation_mutant",
        "negative_control_description": (
            "phaZ mutant host does not degrade PHA (SwissProt ECO:0000269 complementation target)"
        ),
        "grade_basis": (
            "evidence_type=cloned_activity with SwissProt-curated [FUNCTION] 'complements a mutant "
            "that does not degrade PHA' (ECO:0000269): in-vivo causal complementation evidence -> E2"
        ),
        "note_excerpt": (
            "SwissProt curated [FUNCTION] 'complements a mutant that does not degrade PHA' "
            "(evidence ECO:0000269, direct assay). Strain GPo1 (Ectopseudomonas oleovorans)."
        ),
    },
    "DED_hfam_4_0002": {
        "study_key": "CITATION:GenBank_AY113181_unpublished_PHA_synthase_reference",
        "study_basis": (
            "amendment primary_title/note: the GenBank record links no characterization paper; its "
            "reference is an unpublished PHA-synthase paper, so the record's own provenance is the "
            "study anchor"
        ),
        "grade_basis": (
            "evidence_type=sequence_only, confidence=none: no activity evidence on the record and "
            "no linked characterization paper -> A; this row is one of two duplicate records of the "
            "same gene (seed GI 21689574) and can never add an independent positive"
        ),
        "note_excerpt": (
            "No characterization paper linked on the GenBank record (reference is about PHA "
            "synthase, unpublished). Product is phaZ (PHA depolymerase) but no activity evidence. "
            "Same gene as AAN70570.1 (identical 283 aa)."
        ),
    },
    "DED_hfam_4_0033": {
        "study_key": "PMID:12534463",
        "study_basis": (
            "amendment primary_pmid 12534463 = Nelson 2002, Pseudomonas putida KT2440 complete "
            "genome paper (the only publication linked to this gene)"
        ),
        "grade_basis": (
            "evidence_type=sequence_only, confidence=none: genome annotation only (locus PP_5004, "
            "phaZ) with no activity characterization -> A; same gene as DED_hfam_4_0002"
        ),
        "note_excerpt": (
            "Genome annotation (locus PP_5004, phaZ); no activity characterization. Same gene as "
            "AAM63408.1. Not an independent positive."
        ),
    },
    "DED_hfam_52_0001": {
        "study_key": "PMID:2644188",
        "study_basis": (
            "amendment primary_pmid 2644188 = Saito 1989, 'Cloning, nucleotide sequence, and "
            "expression in Escherichia coli of the gene for poly(3-hydroxybutyrate) depolymerase "
            "from Alcaligenes faecalis'"
        ),
        "grade_basis": (
            "evidence_type=cloned_activity, confidence=high: the cloned gene was expressed and its "
            "product degrades water-insoluble and water-soluble PHB to monomer (SwissProt "
            "[FUNCTION]); purification/kinetics are not documented in the frozen evidence -> E2, not E3"
        ),
        "note_excerpt": (
            "Ralstonia pickettii T1 (= Alcaligenes faecalis T1). SwissProt [FUNCTION] degrades "
            "water-insoluble and water-soluble PHB to monomeric 3-hydroxybutyrate."
        ),
    },
    "DED_hfam_52_0002": {
        "study_key": "PMID:9177489",
        "study_basis": (
            "amendment primary_pmid 9177489 = Kita 1997, 'Cloning of poly(3-hydroxybutyrate) "
            "depolymerase from a marine bacterium, Alcaligenes faecalis AE122, and characterization "
            "of its gene product'"
        ),
        "grade_basis": (
            "evidence_type=cloned_activity, confidence=high: cloned gene product characterized, but "
            "no purified-protein conversion/kinetics documented in the frozen evidence -> E2; strain "
            "AE122 is a distinct gene/source (enzyme properties also in Kita 1995, PMID 7646009)"
        ),
        "note_excerpt": "Strain AE122. Enzyme properties also in Kita et al. 1995 (PMID 7646009).",
    },
    "DED_hfam_52_0003": {
        "study_key": "PMID:10742216",
        "study_basis": (
            "amendment primary_pmid 10742216 = Schober 2000, 'Poly(3-hydroxyvalerate) depolymerase "
            "of Pseudomonas lemoignei'"
        ),
        "grade_basis": (
            "evidence_type=cloned_activity, confidence=high for a distinct gene (PHV depolymerase); "
            "no purified-protein conversion/kinetics documented in the frozen evidence -> E2, not E3"
        ),
        "note_excerpt": "PHV depolymerase (distinct gene).",
    },
    "DED_hfam_52_0004": {
        "study_key": "PMID:7836292",
        "study_basis": (
            "amendment primary_pmid 7836292 = Jendrossek 1994, 'Biochemical and molecular "
            "characterization of the Pseudomonas lemoignei polyhydroxyalkanoate depolymerase system'"
        ),
        "grade_basis": (
            "evidence_type=cloned_activity, confidence=high (depolymerase A of the depolymerase A/B "
            "system); same primary paper as DED_hfam_52_0007 -> one independence group, and no "
            "purified-protein conversion/kinetics documented in the frozen evidence -> E2, not E3"
        ),
        "note_excerpt": (
            "Same paper as AAA65705.1 (hfam_52_0007); depolymerase A/B system."
        ),
    },
    "DED_hfam_52_0005": {
        "study_key": "CITATION:GenBank_BAA82057.1_Zhang_Saito_Ralstonia_pickettii_A1",
        "study_basis": (
            "amendment note: GenBank AUTHORS Zhang K, Saito T; the record has no PUBMED link and was "
            "not found in PubMed/Europe PMC (older, non-indexed journal), so the record itself is the "
            "documented primary source"
        ),
        "grade_basis": (
            "evidence_type=cloned_activity, confidence=medium: a cloning/activity claim for a distinct "
            "strain (A1), so E2 rather than E1; the bibliographic gap (no PubMed-indexed primary "
            "citation) is recorded in study_id as a CITATION: anchor and is NOT converted into a grade "
            "change"
        ),
        "note_excerpt": (
            "GenBank AUTHORS Zhang K, Saito T; no PUBMED on record and not found in PubMed/Europe PMC "
            "(older journal). Distinct strain A1 (distinct gene from P12625/T1 and BAA04986/K1)."
        ),
    },
    "DED_hfam_52_0006": {
        "study_key": "CITATION:GenBank_BAA04986.1_Yukawa_Uchida_Ralstonia_pickettii_K1",
        "study_basis": (
            "amendment note: GenBank AUTHORS Yukawa H, Uchida Y, ...; no PMID found, so the record "
            "itself is the documented primary source"
        ),
        "grade_basis": (
            "evidence_type=cloned_activity, confidence=medium: a cloning/activity claim for a distinct "
            "strain (K1), so E2; the missing PubMed-indexed citation is recorded as a CITATION: study "
            "anchor, not as a grade change"
        ),
        "note_excerpt": (
            "GenBank AUTHORS Yukawa H, Uchida Y, ...; no PMID found. Distinct strain K1. Cloning "
            "context = biodegradability-monitoring reporter."
        ),
    },
    "DED_hfam_52_0007": {
        "study_key": "PMID:7836292",
        "study_basis": (
            "amendment primary_pmid 7836292 = Jendrossek 1994 (same primary paper as DED_hfam_52_0004)"
        ),
        "grade_basis": (
            "evidence_type=cloned_activity, confidence=high (depolymerase A precursor): a different "
            "gene from DED_hfam_52_0004 but the SAME primary paper -> same independence group, so it "
            "adds no second independent positive (amendment note: 'NOT an independent paper from row 9')"
        ),
        "note_excerpt": (
            "Depolymerase A precursor. Same paper as AAA65703.1 (hfam_52_0004) -> NOT an independent "
            "paper from row 9."
        ),
    },
    "DED_hfam_53_0001": {
        "study_key": "PMID:9872779",
        "study_basis": (
            "amendment primary_pmid 9872779 = Ohura 1999, 'Cloning and characterization of the "
            "polyhydroxybutyrate depolymerase gene of Pseudomonas stutzeri'"
        ),
        "grade_basis": (
            "evidence_type=cloned_activity, confidence=high: the gene product was characterized after "
            "cloning, but no purified-protein conversion/kinetics documented in the frozen evidence "
            "-> E2, not E3"
        ),
        "note_excerpt": "Same gene as ACG63775.1 (hfam_53_0008).",
    },
    "DED_hfam_53_0008": {
        "study_key": "CITATION:GenBank_2008_Yoon_Rho_Baek_direct_submission",
        "study_basis": (
            "amendment primary_title/note: a GenBank direct submission (Yoon, Rho & Baek, 2008) "
            "re-sequencing the same gene; it has no independent paper"
        ),
        "grade_basis": (
            "evidence_type=same_gene_duplicate, confidence=none: the record is a second database entry "
            "of the gene already characterized by DED_hfam_53_0001 (same seed GI 75538924) -> A, no "
            "independent positive"
        ),
        "note_excerpt": (
            "GenBank re-sequencing of the P. stutzeri PHB depolymerase (same gene as O82950). No "
            "independent paper. NOT an independent positive."
        ),
    },
    "DED_hfam_55_0001": {
        "study_key": "PMID:8269961",
        "study_basis": (
            "amendment primary_pmid 8269961 = Jendrossek 1993, 'Cloning and characterization of the "
            "poly(hydroxyalkanoic acid) depolymerase gene locus of Pseudomonas lemoignei'"
        ),
        "grade_basis": (
            "evidence_type=cloned_activity, confidence=high (PHA depolymerase C / PhaZ C): the gene "
            "locus was cloned and characterized, but no purified-protein conversion/kinetics in the "
            "frozen evidence -> E2, not E3"
        ),
        "note_excerpt": "PHA depolymerase C (PhaZ C).",
    },
    "DED_hfam_55_0002": {
        "study_key": "CITATION:Briese_1994_J_Environ_Polym_Degrad_2_75",
        "study_basis": (
            "amendment primary_title: Briese, Schmidt & Jendrossek 1994, 'Pseudomonas lemoignei has "
            "five poly(hydroxyalkanoic acid) (PHA) depolymerase genes: a comparative study of the "
            "enzymes', J Environ Polym Degrad 2:75 -- journal not PubMed-indexed (no PMID)"
        ),
        "grade_basis": (
            "evidence_type=cloned_activity, confidence=medium (depolymerase D): a comparative enzyme "
            "characterization, so E2; the same non-indexed paper also reports depolymerase B "
            "(DED_hfam_55_0003), which therefore shares this row's independence group"
        ),
        "note_excerpt": (
            "Depolymerase D. Journal not PubMed-indexed (no PMID). Same paper as AAB17150.1 (dep B)."
        ),
    },
    "DED_hfam_55_0003": {
        "study_key": "CITATION:Briese_1994_J_Environ_Polym_Degrad_2_75",
        "study_basis": (
            "amendment primary_title: Briese, Schmidt & Jendrossek 1994 (same non-indexed paper as "
            "DED_hfam_55_0002)"
        ),
        "grade_basis": (
            "evidence_type=cloned_activity, confidence=medium (depolymerase B): a different gene from "
            "DED_hfam_55_0002 but the SAME primary paper -> same independence group, so it is not a "
            "second independent positive"
        ),
        "note_excerpt": "Depolymerase B. Same paper as AAB48166.1 (dep D).",
    },
    "DED_hfam_58_0001": {
        "study_key": "PMID:17064368",
        "study_basis": (
            "amendment primary_pmid 17064368 = Takaku 2006, 'Isolation of a Gram-positive "
            "poly(3-hydroxybutyrate)-degrading bacterium and cloning of its depolymerase'"
        ),
        "grade_basis": (
            "evidence_type=cloned_activity, confidence=high: the depolymerase was cloned from "
            "Priestia (Bacillus) megaterium N-18-25-9, but no purified-protein conversion/kinetics in "
            "the frozen evidence -> E2, not E3; distinct from the intracellular B. megaterium "
            "depolymerase of PMID 19561190"
        ),
        "note_excerpt": (
            "Priestia megaterium (Bacillus megaterium) N-18-25-9. Distinct from the intracellular "
            "B. megaterium depolymerase of Chen et al. 2009 (PMID 19561190)."
        ),
    },
    "DED_hfam_70_0001": {
        "study_key": "CITATION:Kobayashi_1999_J_Environ_Polym_Degrad_7_281",
        "study_basis": (
            "amendment primary_title: Kobayashi 1999, 'Biochemical and genetical characterization of "
            "an extracellular poly(3-hydroxybutyrate) depolymerase from Acidovorax sp. strain TP4', "
            "J Environ Polym Degrad 7:281 -- journal not PubMed-indexed (no PMID)"
        ),
        "grade_basis": (
            "evidence_type=cloned_activity, confidence=medium: a biochemical and genetical "
            "characterization of the strain TP4 depolymerase, so E2; the missing PubMed-indexed "
            "citation is recorded as a CITATION: study anchor, not as a grade change"
        ),
        "note_excerpt": "Journal not PubMed-indexed (no PMID).",
    },
    "DED_hfam_70_0002": {
        "study_key": "PMID:16232882",
        "study_basis": (
            "amendment primary_pmid 16232882 = Takeda 2006, 'Cloning and expression of the gene "
            "encoding thermostable poly(3-hydroxybutyrate) depolymerase'"
        ),
        "grade_basis": (
            "evidence_type=cloned_activity, confidence=high (thermostable PHB depolymerase): no "
            "purified-protein conversion/kinetics documented in the frozen evidence -> E2, not E3"
        ),
        "note_excerpt": "Thermostable PHB depolymerase.",
    },
    "DED_hfam_70_0003": {
        "study_key": "PMID:7606660",
        "study_basis": (
            "amendment primary_pmid 7606660 = Jendrossek 1995, 'Characterization of the extracellular "
            "poly(3-hydroxybutyrate) depolymerase of Comamonas sp.'"
        ),
        "grade_basis": (
            "evidence_type=cloned_activity, confidence=high (strain DSM 6781): no purified-protein "
            "conversion/kinetics documented in the frozen evidence -> E2, not E3"
        ),
        "note_excerpt": "Strain DSM 6781.",
    },
    "DED_hfam_70_0004": {
        "study_key": "PMID:9406404",
        "study_basis": (
            "amendment primary_pmid 9406404 = Kasuya 1997, 'Biochemical and molecular characterization "
            "of the poly(3-hydroxybutyrate) depolymerase of Delftia acidovorans'"
        ),
        "grade_basis": (
            "evidence_type=cloned_activity, confidence=high (SwissProt curated): no purified-protein "
            "conversion/kinetics documented in the frozen evidence -> E2, not E3"
        ),
        "note_excerpt": "SwissProt curated.",
    },
    "DED_hfam_70_0005": {
        "study_key": "PMID:15340791",
        "study_basis": (
            "amendment primary_pmid 15340791 = Romen 2004, 'Thermotolerant poly(3-hydroxybutyrate)-"
            "degrading bacteria from hot compost and characterization of the PHB depolymerase of "
            "Schlegelella sp. KB1a'"
        ),
        "grade_basis": (
            "evidence_type=cloned_activity, confidence=high (thermotolerant KB1a depolymerase): no "
            "purified-protein conversion/kinetics documented in the frozen evidence -> E2, not E3"
        ),
        "note_excerpt": "Thermotolerant.",
    },
    "DED_hfam_70_0006": {
        "study_key": "PMID:8810505",
        "study_basis": (
            "amendment primary_pmid 8810505 = Klingbeil 1996, 'Taxonomic identification of Streptomyces "
            "exfoliatus K10 and characterization of its poly(3-hydroxybutyrate) depolymerase gene'"
        ),
        "grade_basis": (
            "evidence_type=cloned_activity, confidence=high (strain K10): no purified-protein "
            "conversion/kinetics documented in the frozen evidence -> E2, not E3"
        ),
        "note_excerpt": "Strain K10.",
    },
    "DED_hfam_70_0016": {
        "study_key": "PMID:16405909",
        "study_basis": (
            "amendment primary_pmid 16405909 = Hisano 2006, 'The crystal structure of "
            "polyhydroxybutyrate depolymerase from Penicillium funiculosum provides insights into the "
            "recognition and degradation of biopolyesters'"
        ),
        "grade_basis": (
            "evidence_type=purified_activity, confidence=high: purified and crystallised enzyme with "
            "documented substrate recognition/degradation -> E3"
        ),
        "note_excerpt": (
            "seed_gi 88192747 = PDB 2D80 chain A (Talaromyces funiculosus = teleomorph of Penicillium "
            "funiculosum). Ledger accession BAG32152.1 is the protein, seed GI is the structure chain "
            "-> same enzyme. Purified + crystallized."
        ),
    },
    "DED_hfam_7_0001": {
        "study_key": "PMID:11457823",
        "study_basis": (
            "amendment primary_pmid 11457823 = Handrick 2001, 'A new type of thermoalkalophilic "
            "hydrolase of Paucimonas lemoignei with high specificity for amorphous polyesters of short "
            "chain-length hydroxyalkanoic acids'"
        ),
        "grade_basis": (
            "evidence_type=cloned_activity, confidence=high (PhaZ7, AHSMG motif): no purified-protein "
            "conversion/kinetics documented in the frozen evidence -> E2, not E3"
        ),
        "note_excerpt": "PhaZ7 (AHSMG motif).",
    },
    "DED_hfam_8_0001": {
        "study_key": "PMID:15995648",
        "study_basis": (
            "amendment primary_pmid 15995648 = Kim 2005, 'Molecular characterization of extracellular "
            "medium-chain-length poly(3-hydroxyalkanoate) depolymerase genes from Pseudomonas "
            "alcaligenes strains'"
        ),
        "grade_basis": (
            "evidence_type=cloned_activity, confidence=high (strain M4-7): no purified-protein "
            "conversion/kinetics documented in the frozen evidence -> E2, not E3; the bound sequence "
            "contains an ambiguous residue, so the row cannot be used for functional calibration"
        ),
        "note_excerpt": "Strain M4-7. Same paper as AAO73963.1 (LB19).",
    },
    "DED_hfam_8_0002": {
        "study_key": "PMID:15995648",
        "study_basis": (
            "amendment primary_pmid 15995648 = Kim 2005 (same primary paper as DED_hfam_8_0001)"
        ),
        "grade_basis": (
            "evidence_type=cloned_activity, confidence=high (strain LB19): a different gene from "
            "DED_hfam_8_0001 but the SAME primary paper -> same independence group, so it is not a "
            "second independent positive"
        ),
        "note_excerpt": "Strain LB19. Same paper as AAQ72538.1 (M4-7).",
    },
    "DED_hfam_8_0003": {
        "study_key": "PMID:7961472",
        "study_basis": (
            "amendment primary_pmid 7961472 = Schirmer & Jendrossek 1994, 'Molecular characterization "
            "of the extracellular poly(3-hydroxyoctanoic acid) [P(3HO)] depolymerase gene of "
            "Pseudomonas fluorescens GK13 and of its gene product'"
        ),
        "grade_basis": (
            "evidence_type=cloned_activity, confidence=high (SwissProt 'Evidence at protein level'): "
            "no purified-protein conversion/kinetics documented in the frozen evidence -> E2, not E3"
        ),
        "note_excerpt": "SwissProt 'Evidence at protein level'.",
    },
    "DED_hfam_3_0001": {
        "study_key": "PMID:15489436",
        "study_basis": (
            "amendment primary_pmid 15489436 = Handrick 2004, \"The 'intracellular' "
            "poly(3-hydroxybutyrate) (PHB) depolymerase of Rhodospirillum rubrum is a periplasm-located "
            "protein with specificity for native PHB ...\""
        ),
        "localization_evidence": (
            "experimental_periplasmic_localization_demonstrated_in_the_primary_study"
        ),
        "grade_basis": (
            "evidence_type=cloned_activity, confidence=high: the primary study establishes the "
            "periplasmic localization and native-PHB specificity of the gene product; no "
            "purified-protein conversion/kinetics documented in the frozen evidence -> E2, not E3"
        ),
        "note_excerpt": (
            "This is the paper that establishes periplasmic localization (correcting the "
            "'intracellular' label)."
        ),
    },
}


# --- small helpers --------------------------------------------------------------


def _text(value: Any) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if value is None:
        return ""
    return str(value).strip()


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _slug(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_")


def _regular_file(path: Path, label: str) -> Path:
    if not path.is_file():
        raise ValueError(f"{label} is missing or not a regular file: {path}")
    return path


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def read_tsv(path: Path, label: str) -> tuple[list[str], list[dict[str, str]]]:
    """Read a TSV with a checked header; a missing file is a refused input."""
    _regular_file(path, label)
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        rows = [dict(row) for row in reader]
    if not fields:
        raise ValueError(f"{label} has no header: {path}")
    return fields, rows


def read_fasta(path: Path) -> dict[str, tuple[str, str]]:
    """Read ``>reference_id|accession`` records into ``{reference_id: (accession, sequence)}``."""
    _regular_file(path, "reference FASTA")
    records: dict[str, tuple[str, str]] = {}
    current: str | None = None
    accession: str = ""
    chunks: list[str] = []
    for raw in path.read_text(encoding="ascii").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current is not None:
                records[current] = (accession, "".join(chunks))
            header = line[1:].strip()
            tokens = header.split("|")
            if len(tokens) != 2 or not tokens[0].strip() or not tokens[1].strip():
                raise ValueError(f"reference FASTA header must be '>reference_id|accession': {header!r}")
            current, accession, chunks = tokens[0].strip(), tokens[1].strip(), []
        else:
            if current is None:
                raise ValueError("reference FASTA sequence appears before any header")
            chunks.append(line)
    if current is not None:
        records[current] = (accession, "".join(chunks))
    if not records:
        raise ValueError(f"reference FASTA has no records: {path}")
    return records


# --- validator reuse ------------------------------------------------------------


_VALIDATOR_CACHE: dict[str, Any] = {}


def load_validator(path: str | Path | None = None):
    """Import the frozen validator module.

    The validator owns the evidence-column contract, the vocabulary and the
    independence counting, so the curation reuses it instead of restating those
    rules.  A missing validator is a refused input, not a silent fallback.
    """
    resolved = Path(path) if path is not None else VALIDATOR_PATH
    key = str(resolved)
    if key not in _VALIDATOR_CACHE:
        _regular_file(resolved, "validator")
        spec = importlib.util.spec_from_file_location("validate_phaded_reference_ledger", resolved)
        if spec is None or spec.loader is None:
            raise ValueError(f"cannot import validator: {resolved}")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        _VALIDATOR_CACHE[key] = module
    return _VALIDATOR_CACHE[key]


def load_evidence_columns(path: str | Path | None = None) -> tuple[str, ...]:
    """Read and cross-check the evidence-column contract via the validator."""
    validator = load_validator()
    config = Path(path) if path is not None else EVIDENCE_COLUMNS_CONFIG
    try:
        return tuple(validator.load_evidence_columns(config))
    except ValueError as error:
        raise ValueError(f"evidence columns config is invalid: {error}") from error


# --- evidence helpers -----------------------------------------------------------


def parse_seed_gi(notes: Any) -> str:
    """Return the Knoll 2009 Table 1 seed GI recorded in a frozen ledger ``notes`` cell."""
    match = _SEED_GI_PATTERN.search(_text(notes))
    return match.group(1) if match else ""


def sequence_integrity(sequence: Any) -> str:
    """Map the project's FASTA-integrity rule onto the evidence vocabulary.

    Mirrors ``build_phaded_evidence_tiers.normalized_sequence_integrity`` /
    ``annotate_phaded_features`` (clean alphabet + initiator Met == valid) and maps
    ``valid``/``terminal_stop_only`` -> ``complete``, ``possible_N_truncation`` ->
    ``partial``, ``invalid_internal_character`` -> ``unresolved``.
    """
    seq = _text(sequence).upper()
    if not seq:
        return "unresolved"
    if seq.endswith("*") and seq.count("*") == 1:
        seq = seq[:-1]
        if seq and all(residue in AMINO_ACIDS for residue in seq):
            return "complete" if seq.startswith("M") else "partial"
    if any(residue not in AMINO_ACIDS for residue in seq):
        return "unresolved"
    if not seq.startswith("M"):
        return "partial"
    return "complete"


def family_call_consistent(row: Mapping[str, Any]) -> bool:
    """Whether the frozen family call is a source-bound PhaDED family (validator vocabulary)."""
    validator = load_validator()
    superfamily = _text(row.get("phaded_superfamily"))
    family = _text(row.get("phaded_family_id"))
    if superfamily not in validator.ALLOWED_SUPERFAMILIES:
        return False
    if not family or family.lower() in {"unresolved", "unassigned"}:
        return False
    return True


def grade_for_evidence_type(evidence_type: Any) -> str:
    """Documented grade for one amendment ``evidence_type`` (fail-closed on unknown values)."""
    key = _text(evidence_type)
    if key not in EVIDENCE_TYPE_RULES:
        raise ValueError(
            f"undocumented amendment evidence_type {key!r}; documented values: "
            + ",".join(sorted(EVIDENCE_TYPE_RULES))
        )
    return EVIDENCE_TYPE_RULES[key]["grade"]


def accession_version_for(ledger_accession: Any, provenance_accession: Any = "") -> tuple[str, str]:
    """Return ``(accession_version, basis)`` without inventing a version.

    The frozen ``accession`` value is authoritative.  When the amendment records the
    same accession with a version suffix the amendment value is used and the basis
    says so; a genuine accession disagreement is refused by the caller.
    """
    ledger_value = _text(ledger_accession)
    provenance_value = _text(provenance_accession)
    if provenance_value and provenance_value != ledger_value:
        ledger_match = _VERSION_SUFFIX.match(ledger_value)
        provenance_match = _VERSION_SUFFIX.match(provenance_value)
        ledger_base = ledger_match.group("base") if ledger_match else ledger_value
        provenance_base = provenance_match.group("base") if provenance_match else provenance_value
        if provenance_base and provenance_base == ledger_base:
            return provenance_value, "amendment_versioned_accession"
    return ledger_value, "frozen_accession_verbatim"


def experimental_substrate_class(row: Mapping[str, Any]) -> str:
    """Compose the frozen substrate fields into one experimental substrate token."""
    substrate = _text(row.get("substrate_class"))
    detail = _text(row.get("substrate_detail"))
    if not substrate:
        return PENDING
    token = _slug(f"{substrate}_{detail}") if detail else _slug(substrate)
    return token or PENDING


# --- independence groups --------------------------------------------------------


def independence_groups(members: Sequence[Mapping[str, str]]) -> dict[str, dict[str, Any]]:
    """Union positives by documented same-gene (seed GI) and same-primary-study bindings.

    ``members`` carry ``reference_id``, ``seed_gi`` and ``study_key``.  Returns
    ``{reference_id: {"group_id":..., "members": [...], "basis": [...], "binding": ...}}``.
    """
    parent = {member["reference_id"]: member["reference_id"] for member in members}

    def find(item: str) -> str:
        while parent[item] != item:
            parent[item] = parent[parent[item]]
            item = parent[item]
        return item

    def union(left: str, right: str) -> None:
        left_root, right_root = find(left), find(right)
        if left_root != right_root:
            # Deterministic: the first row in ledger order wins the root.
            parent[max(left_root, right_root)] = min(left_root, right_root)

    bindings: dict[str, list[str]] = {}
    for label, selector in (("seed_gi", "seed_gi"), ("study_key", "study_key")):
        buckets: dict[str, list[Mapping[str, str]]] = {}
        for member in members:
            value = _text(member.get(selector))
            if value:
                buckets.setdefault(value, []).append(member)
        for value, group in buckets.items():
            if len(group) < 2:
                continue
            references = sorted(member["reference_id"] for member in group)
            for member in group[1:]:
                union(group[0]["reference_id"], member["reference_id"])
            bindings.setdefault(label, []).append(
                f"same_{label}={value} binds {', '.join(references)}"
            )

    components: dict[str, list[Mapping[str, str]]] = {}
    for member in members:
        components.setdefault(find(member["reference_id"]), []).append(member)

    result: dict[str, dict[str, Any]] = {}
    for component in components.values():
        members_sorted = sorted(component, key=lambda member: member["reference_id"])
        seed_gis = {_text(member.get("seed_gi")) for member in members_sorted}
        gene_bound = len(members_sorted) > 1 and len(seed_gis) == 1 and "" not in seed_gis
        if gene_bound:
            group_id = f"IG_GI_{next(iter(seed_gis))}"
            binding = "same_gene_seed_gi"
        else:
            anchors = sorted(_text(member.get("study_key")) for member in members_sorted)
            indexed = [anchor for anchor in anchors if anchor.startswith("PMID:")]
            anchor = (indexed or anchors)[0]
            group_id = f"IG_{_slug(anchor)}"
            binding = "same_primary_study" if len(members_sorted) > 1 else "single_primary_study"
        reason = [
            binding,
        ]
        if gene_bound:
            reason.append(f"shared seed_gi={next(iter(seed_gis))}")
        else:
            reason.append("shared primary study anchor=" + _text(members_sorted[0].get("study_key")))
        for member in members_sorted:
            result[member["reference_id"]] = {
                "group_id": group_id,
                "binding": binding,
                "members": [item["reference_id"] for item in members_sorted],
                "seed_gis": sorted(seed_gis),
                "basis": "; ".join(reason),
            }
    relevant = []
    for label, entries in sorted(bindings.items()):
        relevant.extend(entries)
    for entry in result.values():
        entry["documented_bindings"] = list(relevant)
    return result


# --- curation -------------------------------------------------------------------


def assert_append_only(
    original_fields: Sequence[str],
    original_rows: Sequence[Mapping[str, Any]],
    curated_rows: Sequence[Mapping[str, Any]],
) -> None:
    """Refuse any change to an original column value (append-only guarantee)."""
    if len(original_rows) != len(curated_rows):
        raise ValueError(
            "append-only violation: row count changed from "
            f"{len(original_rows)} to {len(curated_rows)}"
        )
    for index, (original, curated) in enumerate(zip(original_rows, curated_rows), start=2):
        label = _text(original.get("reference_id")) or f"row {index}"
        if _text(curated.get("reference_id")) != _text(original.get("reference_id")):
            raise ValueError(
                f"append-only violation at row {index}: reference_id changed from "
                f"{_text(original.get('reference_id'))!r} to {_text(curated.get('reference_id'))!r}"
            )
        for field in original_fields:
            # Raw comparison on purpose: even a whitespace-only change is a violation.
            if curated.get(field) != original.get(field):
                raise ValueError(
                    f"append-only violation at row {index} ({label}): column {field} changed from "
                    f"{original.get(field)!r} to {curated.get(field)!r}"
                )


def build_curation(
    ledger_fields: Sequence[str],
    ledger_rows: Sequence[Mapping[str, str]],
    provenance_fields: Sequence[str],
    provenance_rows: Sequence[Mapping[str, str]],
    sequences: Mapping[str, tuple[str, str]],
    evidence_columns: Sequence[str],
) -> dict[str, Any]:
    """Build the curated rows, the derivation table and the pending-value classes."""
    missing_ledger = [column for column in FROZEN_LEDGER_COLUMNS if column not in ledger_fields]
    if missing_ledger:
        raise ValueError("frozen ledger is missing required columns: " + ",".join(missing_ledger))
    missing_provenance = [column for column in PROVENANCE_COLUMNS if column not in provenance_fields]
    if missing_provenance:
        raise ValueError(
            "provenance amendment is missing required columns: " + ",".join(missing_provenance)
        )
    if not provenance_rows:
        raise ValueError("provenance amendment has no data rows")

    provenance_by_reference: dict[str, Mapping[str, str]] = {}
    for row in provenance_rows:
        reference = _text(row.get("reference_id"))
        if not reference:
            raise ValueError("provenance amendment has a row without reference_id")
        if reference in provenance_by_reference:
            raise ValueError(f"provenance amendment has duplicate reference_id: {reference}")
        provenance_by_reference[reference] = row

    positives = [row for row in ledger_rows if _text(row.get("evidence_status")) == EVIDENCE_STATUS_POSITIVE]
    positive_references = {_text(row.get("reference_id")) for row in positives}
    unsupported = sorted(set(provenance_by_reference) - positive_references)
    if unsupported:
        raise ValueError(
            "provenance amendment documents rows that are not experimental positives in the frozen "
            "ledger (refusing to guess what they belong to): " + ",".join(unsupported)
        )

    members: list[dict[str, str]] = []
    curation_detail: dict[str, dict[str, Any]] = {}
    for row in positives:
        reference = _text(row.get("reference_id"))
        accession = _text(row.get("accession"))
        provenance = provenance_by_reference.get(reference)
        if provenance is None:
            raise ValueError(
                f"experimental positive {reference} ({accession}) has no row in the provenance "
                "amendment; its evidence cannot be graded or independence-grouped and will not be "
                "fabricated"
            )
        if _text(provenance.get("accession")) and _text(provenance.get("accession")) != accession:
            ledger_match = _VERSION_SUFFIX.match(accession)
            provenance_match = _VERSION_SUFFIX.match(_text(provenance.get("accession")))
            ledger_base = ledger_match.group("base") if ledger_match else accession
            provenance_base = provenance_match.group("base") if provenance_match else _text(provenance.get("accession"))
            if ledger_base != provenance_base:
                raise ValueError(
                    f"provenance amendment accession disagrees with the frozen ledger for {reference}: "
                    f"ledger={accession!r} amendment={_text(provenance.get('accession'))!r}"
                )
        for column in ("phaded_superfamily", "phaded_family_id", "organism"):
            if _text(provenance.get(column)) and _text(provenance.get(column)) != _text(row.get(column)):
                raise ValueError(
                    f"provenance amendment {column} disagrees with the frozen ledger for {reference}: "
                    f"ledger={_text(row.get(column))!r} amendment={_text(provenance.get(column))!r}"
                )
        seed_gi = _text(provenance.get("seed_gi"))
        note_seed_gi = parse_seed_gi(row.get("notes"))
        if note_seed_gi and seed_gi and note_seed_gi != seed_gi:
            raise ValueError(
                f"{reference}: seed GI disagrees between the frozen ledger notes ({note_seed_gi}) and "
                f"the provenance amendment ({seed_gi})"
            )
        if not seed_gi:
            seed_gi = note_seed_gi
        if not seed_gi:
            raise ValueError(
                f"{reference}: no documented seed GI in either the frozen ledger notes or the "
                "provenance amendment; the gene identity (and therefore independence) cannot be "
                "derived and will not be invented"
            )
        evidence_type = _text(provenance.get("evidence_type"))
        grade = grade_for_evidence_type(evidence_type)
        curation = POSITIVE_CURATION.get(reference)
        if curation is None:
            raise ValueError(
                f"{reference}: no curation entry in POSITIVE_CURATION for evidence_type "
                f"{evidence_type!r}; add the documented study anchor and grade basis before curating"
            )
        study_key = _text(curation.get("study_key"))
        if not study_key:
            raise ValueError(f"{reference}: POSITIVE_CURATION entry has no study_key")
        primary_pmid = _text(provenance.get("primary_pmid"))
        primary_doi = _text(provenance.get("primary_doi"))
        primary_title = _text(provenance.get("primary_title"))
        if study_key.startswith("PMID:"):
            declared_pmid = study_key.split(":", 1)[1]
            if not primary_pmid:
                raise ValueError(
                    f"{reference}: the curation table anchors an indexed study {study_key} but the "
                    "amendment records no primary_pmid"
                )
            if primary_pmid != declared_pmid:
                raise ValueError(
                    f"{reference}: study anchor disagrees with the amendment: table={study_key} "
                    f"amendment primary_pmid={primary_pmid}"
                )
            study_id_basis = f"amendment primary_pmid={primary_pmid}"
            if primary_doi:
                study_id_basis += f" (primary_doi={primary_doi})"
        elif study_key.startswith("CITATION:"):
            if primary_pmid:
                raise ValueError(
                    f"{reference}: the amendment records primary_pmid={primary_pmid}, so the "
                    f"non-indexed citation anchor {study_key} is wrong"
                )
            study_id_basis = (
                "amendment primary_pmid and primary_doi are blank; documented non-indexed source: "
                + (primary_title or "see amendment note")
            )
        else:
            raise ValueError(
                f"{reference}: undocumented study anchor form {study_key!r}; documented forms are "
                "PMID:<pmid> and CITATION:<source>"
            )
        members.append({"reference_id": reference, "seed_gi": seed_gi, "study_key": study_key})
        curation_detail[reference] = {
            "provenance": provenance,
            "curation": curation,
            "grade": grade,
            "evidence_type": evidence_type,
            "seed_gi": seed_gi,
            "seed_gi_source": (
                "ledger_note_and_amendment" if note_seed_gi and seed_gi == note_seed_gi
                else ("amendment_only" if seed_gi else "ledger_note_only")
            ),
            "study_key": study_key,
            "study_id_basis": study_id_basis,
        }

    groups = independence_groups(members)

    curated_rows: list[dict[str, str]] = []
    derivation_rows: list[dict[str, str]] = []
    pending_counts: Counter[tuple[str, str]] = Counter()
    pending_reasons: dict[tuple[str, str], str] = {}

    def record_pending(column: str, row_class: str, reason: str) -> None:
        pending_counts[(column, row_class)] += 1
        pending_reasons[(column, row_class)] = reason

    annotation_reason = (
        "annotation_only row: the frozen ledger records no experimental study, experiment unit, "
        "assay or substrate (its `experimental_assay`/`experimental_result` are the frozen `pending` "
        "sentinel and the Knoll 2009 Table 1 seed designation is a database label, not a study); the "
        "AGENTS.md `pending` sentinel is written instead of a fabricated unit"
    )
    annotation_gene_reason = (
        "annotation_only row: the frozen ledger records no per-row gene identifier; a gene id is not "
        "invented"
    )
    frozen_pending_reason = (
        "written by the frozen ledger itself and preserved unchanged by the append-only rule (an "
        "original column value is never altered by curation)"
    )

    for row in ledger_rows:
        reference = _text(row.get("reference_id"))
        status = _text(row.get("evidence_status"))
        curated = {column: row.get(column, "") for column in ledger_fields}
        for column in evidence_columns:
            curated[column] = ""

        integrity = sequence_integrity(sequences[reference][1])
        if integrity not in INTEGRITY_VOCABULARY:
            raise ValueError(f"{reference}: produced an invalid sequence_integrity {integrity!r}")
        consistent = family_call_consistent(row)
        discovery = "true" if (integrity == "complete" and consistent) else "false"

        if status == EVIDENCE_STATUS_POSITIVE:
            detail = curation_detail[reference]
            provenance = detail["provenance"]
            curation = detail["curation"]
            grade = detail["grade"]
            rules = EVIDENCE_TYPE_RULES[detail["evidence_type"]]
            group = groups[reference]
            accession_version, accession_version_basis = accession_version_for(
                row.get("accession"), provenance.get("accession")
            )
            substrate = experimental_substrate_class(row)
            calibration = (
                "true"
                if (grade in {"E2", "E3"} and integrity == "complete" and group["group_id"])
                else "false"
            )
            curated.update({
                "accession_version": accession_version,
                "gene_id": GENE_ID_PREFIX + detail["seed_gi"],
                "study_id": detail["study_key"],
                "experiment_unit_id": f"EU_{group['group_id']}__gi{detail['seed_gi']}",
                "independence_group": group["group_id"],
                "experimental_evidence_grade": grade,
                "assay_directness": rules["assay_directness"],
                "experimental_substrate_class": substrate,
                "experimental_system": _text(curation.get("experimental_system", rules["experimental_system"])),
                "negative_control_description": _text(
                    curation.get("negative_control_description", DEFAULT_NEGATIVE_CONTROL)
                ),
                "catalytic_mutant_evidence": _text(
                    curation.get("catalytic_mutant_evidence", DEFAULT_CATALYTIC_MUTANT)
                ),
                "localization_evidence": _text(
                    curation.get("localization_evidence", DEFAULT_LOCALIZATION_EVIDENCE)
                ),
                "sequence_integrity": integrity,
                "discovery_training_eligible": discovery,
                "functional_calibration_eligible": calibration,
            })
            counts_as_positive = (
                grade in {"E2", "E3"} and rules["assay_directness"] == "direct" and bool(group["group_id"])
            )
            derivation_rows.append({
                "reference_id": reference,
                "accession": _text(row.get("accession")),
                "accession_version": accession_version,
                "accession_version_basis": accession_version_basis,
                "gene_id": curated["gene_id"],
                "seed_gi": detail["seed_gi"],
                "seed_gi_source": detail["seed_gi_source"],
                "evidence_type": detail["evidence_type"],
                "confidence": _text(provenance.get("confidence")),
                "primary_pmid": _text(provenance.get("primary_pmid")),
                "primary_doi": _text(provenance.get("primary_doi")),
                "study_key": detail["study_key"],
                "study_id": curated["study_id"],
                "study_id_basis": detail["study_id_basis"],
                "study_basis": _text(curation.get("study_basis")),
                "experiment_unit_id": curated["experiment_unit_id"],
                "independence_group": group["group_id"],
                "independence_binding": group["binding"],
                "independence_members": ",".join(group["members"]),
                "independence_basis": group["basis"],
                "independence_documented_bindings": " | ".join(group["documented_bindings"]),
                "grade": grade,
                "grade_basis": _text(curation.get("grade_basis", rules["reason"])),
                "assay_directness": curated["assay_directness"],
                "experimental_substrate_class": substrate,
                "experimental_substrate_class_basis": (
                    "frozen substrate_class+substrate_detail of the primary study row"
                    if grade in {"E2", "E3"}
                    else "frozen substrate_class+substrate_detail (inherited DED label, no assay)"
                ),
                "experimental_system": curated["experimental_system"],
                "negative_control_description": curated["negative_control_description"],
                "catalytic_mutant_evidence": curated["catalytic_mutant_evidence"],
                "localization_evidence": curated["localization_evidence"],
                "sequence_integrity": integrity,
                "discovery_training_eligible": discovery,
                "functional_calibration_eligible": calibration,
                "counts_as_functional_positive": "true" if counts_as_positive else "false",
                "amendment_note": _text(provenance.get("note")),
            })
        else:
            curated.update({
                "accession_version": accession_version_for(row.get("accession"))[0],
                "gene_id": PENDING,
                "study_id": PENDING,
                "experiment_unit_id": PENDING,
                "independence_group": PENDING,
                "experimental_evidence_grade": "A",
                "assay_directness": "annotation",
                "experimental_substrate_class": PENDING,
                "experimental_system": PENDING,
                "negative_control_description": DEFAULT_NEGATIVE_CONTROL,
                "catalytic_mutant_evidence": DEFAULT_CATALYTIC_MUTANT,
                "localization_evidence": DEFAULT_LOCALIZATION_EVIDENCE,
                "sequence_integrity": integrity,
                "discovery_training_eligible": discovery,
                "functional_calibration_eligible": "false",
            })
            for column, reason in (
                ("gene_id", annotation_gene_reason),
                ("study_id", annotation_reason),
                ("experiment_unit_id", annotation_reason),
                ("independence_group", annotation_reason),
                ("experimental_substrate_class", annotation_reason),
                ("experimental_system", annotation_reason),
            ):
                record_pending(column, EVIDENCE_STATUS_ANNOTATION, reason)
            for column in ("experimental_assay", "experimental_result"):
                if _text(row.get(column)).lower() == PENDING:
                    record_pending(column, "frozen_original_column", frozen_pending_reason)

        curated_rows.append(curated)

    assert_append_only(ledger_fields, ledger_rows, curated_rows)

    pending_classes = [
        {
            "column": column,
            "row_class": row_class,
            "count": count,
            "reason": pending_reasons[(column, row_class)],
        }
        for (column, row_class), count in sorted(pending_counts.items())
    ]
    positives_with_pending = [
        {
            "reference_id": _text(row.get("reference_id")),
            "column": column,
        }
        for row in curated_rows
        if _text(row.get("evidence_status")) == EVIDENCE_STATUS_POSITIVE
        for column in evidence_columns
        if _text(row.get(column)).lower() == PENDING
    ]
    return {
        "curated_rows": curated_rows,
        "derivation_rows": derivation_rows,
        "groups": groups,
        "pending_value_classes": pending_classes,
        "positives_with_pending_values": positives_with_pending,
    }


# --- report ---------------------------------------------------------------------


def _markdown_table(header: Sequence[str], rows: Iterable[Sequence[Any]]) -> str:
    lines = ["| " + " | ".join(header) + " |", "|" + "|".join("---" for _ in header) + "|"]
    for row in rows:
        lines.append("| " + " | ".join(_text(cell).replace("|", "\\|") for cell in row) + " |")
    return "\n".join(lines)


def build_report(
    *,
    run_id: str,
    generated_at: str,
    payload: Mapping[str, Any],
    inputs: Mapping[str, Mapping[str, Any]],
) -> str:
    curated_rows = payload["curated_rows"]
    derivation_rows = payload["derivation_rows"]
    gate = payload["strict_gate"]
    counts = payload["counts"]
    family_counts = payload["family_functional_positive_counts"]
    genus_counts = payload["family_distinct_genus_counts"]
    family_rows = payload["family_row_counts"]
    positive_rows_per_family = payload["family_positive_row_counts"]

    positives = [row for row in curated_rows if _text(row.get("evidence_status")) == EVIDENCE_STATUS_POSITIVE]
    detail_by_reference = {row["reference_id"]: row for row in derivation_rows}
    grade_counts = Counter(row["experimental_evidence_grade"] for row in positives)

    lines: list[str] = []
    lines.append("# PhaDED reference-ledger evidence curation (F13)")
    lines.append("")
    lines.append(f"- run: `{run_id}`")
    lines.append(f"- generated (UTC): `{generated_at}`")
    lines.append(
        "- scope: the frozen 723-row reference ledger (30 `experimental_positive` + 693 "
        "`annotation_only`), curated with the 15 v2 evidence columns"
    )
    lines.append(
        "- authority: the frozen ledger is **append-only**; this run writes a new ledger and never "
        "modifies the frozen one"
    )
    lines.append("")
    lines.append("## 1. Inputs")
    lines.append("")
    lines.append(_markdown_table(
        ["input", "path", "size", "sha256"],
        [
            [name, entry["path"], entry["size"], entry["sha256"]]
            for name, entry in sorted(inputs.items())
        ],
    ))
    lines.append("")
    lines.append("## 2. Method")
    lines.append("")
    lines.append(
        "Grades are applied strictly from the amendment's documented `evidence_type` "
        "(`purified_activity` -> E3, `cloned_activity` -> E2, `sequence_only` / "
        "`same_gene_duplicate` -> A). A `cloned_activity` row is never promoted to E3 because the "
        "frozen curation inputs document a cloning/activity claim, not purified-protein "
        "conversion or kinetics; the single `purified_activity` row is the only E3. No E1 row is "
        "declared: every activity-backed positive is a direct causal demonstration, so labelling "
        "one `indirect` would misstate the frozen evidence."
    )
    lines.append("")
    lines.append(
        "`independence_group` is the connected component of two documented relations: the same "
        "gene (identical Knoll 2009 Table 1 `seed_gi`, cross-checked between the frozen ledger "
        "`notes` and the amendment) and the same primary study (`primary_pmid`/`primary_doi`, or a "
        "`CITATION:` anchor for the documented non-indexed sources). The group id records which "
        "relation bound it: `IG_GI_<seed_gi>` when every member is one gene, otherwise "
        "`IG_<study>`. `functional_positive_count` is always the number of unique independence "
        "groups, never the row count."
    )
    lines.append("")
    lines.append(
        "`sequence_integrity` is evaluated from the bound FASTA with the project's existing rule "
        "(clean 20-residue alphabet + initiator Met = complete; clean but no initiator Met = "
        "partial; ambiguous residue or internal stop = unresolved). "
        "`discovery_training_eligible=true` requires a complete sequence and a consistent family "
        "call; `functional_calibration_eligible=true` requires a direct E2/E3 row with a "
        "documented independence group and a complete sequence."
    )
    lines.append("")
    lines.append("## 3. Grade distribution among the 30 positives")
    lines.append("")
    lines.append(_markdown_table(
        ["grade", "rows"],
        [[grade, grade_counts.get(grade, 0)] for grade in ("E3", "E2", "E1", "A")],
    ))
    lines.append("")
    lines.append(
        f"- direct functional-positive rows (E3/E2 + `assay_directness=direct`): "
        f"{counts['functional_positive_rows']}"
    )
    lines.append(
        f"- unique independence groups among them: {counts['functional_positive_groups']}"
    )
    lines.append(
        f"- rows eligible for functional calibration: {counts['functional_calibration_eligible_rows']}"
        f" (groups: {counts['functional_calibration_eligible_groups']})"
    )
    counting_not_eligible = [
        row for row in derivation_rows
        if row["counts_as_functional_positive"] == "true"
        and row["functional_calibration_eligible"] == "false"
    ]
    if counting_not_eligible:
        lines.append(
            "- direct positive rows that still may **not** be used for functional calibration "
            "(their sequence integrity is not `complete`), listed here so the gap is explicit:"
        )
        lines.append("")
        lines.append(_markdown_table(
            ["reference_id", "sequence_integrity", "sequence_integrity reason"],
            [
                [
                    row["reference_id"], row["sequence_integrity"],
                    "an ambiguous residue in the bound FASTA sequence means the record cannot be "
                    "verified as a complete CDS translation",
                ]
                for row in counting_not_eligible
            ],
        ))
    lines.append("")
    lines.append("## 4. Per-positive curation with the one-line basis")
    lines.append("")
    lines.append(_markdown_table(
        ["reference_id", "accession", "grade", "directness", "independence_group",
         "experiment_unit_id", "study_id", "counts", "one-line basis"],
        [
            [
                row["reference_id"], row["accession"], row["experimental_evidence_grade"],
                row["assay_directness"], row["independence_group"], row["experiment_unit_id"],
                row["study_id"], detail_by_reference[row["reference_id"]]["counts_as_functional_positive"],
                detail_by_reference[row["reference_id"]]["grade_basis"],
            ]
            for row in positives
        ],
    ))
    lines.append("")
    lines.append("### 4.1 Non-E3 downgrade basis (one line per row)")
    lines.append("")
    downgrades = [row for row in positives if row["experimental_evidence_grade"] != "E3"]
    lines.append(_markdown_table(
        ["reference_id", "grade", "why not E3"],
        [
            [
                row["reference_id"], row["experimental_evidence_grade"],
                detail_by_reference[row["reference_id"]]["grade_basis"],
            ]
            for row in downgrades
        ],
    ))
    lines.append("")
    lines.append("## 5. Independence-group derivation")
    lines.append("")
    grouped: dict[str, list[dict[str, str]]] = {}
    for row in derivation_rows:
        grouped.setdefault(row["independence_group"], []).append(row)
    lines.append(_markdown_table(
        ["independence_group", "binding", "members (reference_id / seed_gi / study)",
         "counting rows", "basis"],
        [
            [
                group_id,
                members[0]["independence_binding"],
                "; ".join(
                    f"{member['reference_id']} / gi{member['seed_gi']} / {member['study_key']}"
                    for member in members
                ),
                str(sum(1 for member in members if member["counts_as_functional_positive"] == "true")),
                members[0]["independence_basis"],
            ]
            for group_id, members in sorted(grouped.items())
        ],
    ))
    lines.append("")
    lines.append(
        f"Every group contributes exactly **1** to `functional_positive_count` (= "
        f"{counts['functional_positive_groups']} here), whatever its counting-row count; a group with "
        "0 counting rows is a duplicate-record group that adds no positive at all."
    )
    lines.append("")
    lines.append("### 5.1 Documented bindings used by the union")
    lines.append("")
    bindings: list[str] = []
    for member in derivation_rows:
        for entry in member["independence_documented_bindings"].split(" | "):
            if entry and entry not in bindings:
                bindings.append(entry)
    if bindings:
        for entry in bindings:
            lines.append(f"- {entry}")
    else:
        lines.append(
            "- no two positives share a gene or a primary study in this ledger, so every group is a "
            "single study"
        )
    lines.append("")
    lines.append("### 5.2 Study anchors")
    lines.append("")
    lines.append(_markdown_table(
        ["reference_id", "study_id", "anchor basis", "amendment study basis"],
        [
            [
                row["reference_id"], row["study_id"], row["study_id_basis"], row["study_basis"],
            ]
            for row in derivation_rows
        ],
    ))
    lines.append("")
    lines.append("## 6. Per-family independent-positive counts")
    lines.append("")
    lines.append(_markdown_table(
        ["phaded_family_id", "ledger rows", "positive rows", "unique independent positives",
         "distinct genera", "meets >= 3 bar"],
        [
            [
                family,
                family_rows.get(family, 0),
                positive_rows_per_family.get(family, 0),
                count,
                genus_counts.get(family, 0),
                "yes" if count >= payload["minimum_independent_positive_count"] else "no",
            ]
            for family, count in sorted(
                family_counts.items(),
                key=lambda item: (-item[1], item[0]),
            )
        ],
    ))
    lines.append("")
    lines.append(
        "Families meeting the >= "
        f"{payload['minimum_independent_positive_count']} independent-positive bar: "
        + (", ".join(payload["families_meeting_minimum_independent_positives"]) or "none")
    )
    lines.append("")
    for family in ("DED_hfam_4", "DED_hfam_53", "DED_hfam_52", "DED_hfam_70"):
        lines.append(
            f"- `{family}`: {family_rows.get(family, 0)} ledger rows, "
            f"{positive_rows_per_family.get(family, 0)} `experimental_positive` rows, "
            f"unique independent positives = **{family_counts.get(family, 0)}**, "
            f"distinct genera = {genus_counts.get(family, 0)}"
        )
    lines.append("")
    lines.append("## 7. Values left `pending`")
    lines.append("")
    lines.append(_markdown_table(
        ["column", "row class", "rows", "reason"],
        [
            [entry["column"], entry["row_class"], entry["count"], entry["reason"]]
            for entry in payload["pending_value_classes"]
        ],
    ))
    lines.append("")
    lines.append(
        f"- values left `pending` on the 30 `experimental_positive` rows: "
        f"{len(payload['positives_with_pending_values'])}"
    )
    if payload["positives_with_pending_values"]:
        lines.append(_markdown_table(
            ["reference_id", "column"],
            [
                [entry["reference_id"], entry["column"]]
                for entry in payload["positives_with_pending_values"]
            ],
        ))
    else:
        lines.append(
            "- every evidence column of every `experimental_positive` row carries a declared value; "
            "the 693 `annotation_only` rows carry the `pending` sentinel because no experimental "
            "unit exists for them"
        )
    lines.append("")
    lines.append("## 8. Acceptance gate (strict evidence-column mode)")
    lines.append("")
    lines.append(f"- validator: `{gate['validator']}`")
    lines.append(f"- mode: `{gate['evidence_column_mode']}`")
    lines.append(f"- status: `{gate['status']}`")
    lines.append(f"- rows: {gate['row_count']}")
    lines.append(f"- evidence columns declared: {', '.join(gate['evidence_columns_declared'])}")
    lines.append(f"- `pending_evidence_columns`: {gate['pending_evidence_columns'] or 'none'}")
    lines.append(f"- `functional_positive_count`: **{gate['functional_positive_count']}**")
    lines.append(
        f"- `family_functional_positive_counts`: "
        f"`{json.dumps(gate['family_functional_positive_counts'], sort_keys=True)}`"
    )
    lines.append(
        f"- `family_distinct_genus_counts`: "
        f"`{json.dumps(gate['family_distinct_genus_counts'], sort_keys=True)}`"
    )
    lines.append("")
    lines.append(
        "*This curation records evidence grades for a reference panel. It is candidate-only: no "
        "grade here asserts a verified PHB/PHA degradation phenotype for any GTDB candidate.*"
    )
    lines.append("")
    return "\n".join(lines)


# --- run ------------------------------------------------------------------------


def run(
    *,
    ledger: str | Path,
    provenance: str | Path,
    family_definitions: str | Path,
    reference_fasta: str | Path,
    output_dir: str | Path,
    evidence_columns: str | Path | None = None,
    derivation_table: str | Path | None = None,
    input_contract: str | Path | None = None,
    report: str | Path | None = None,
    validator_path: str | Path | None = None,
    expected_rows: int | None = None,
    expected_positives: int | None = None,
    expected_family_count: int = 38,
    run_id: str = "",
    generated_at: str | None = None,
    write_report: bool = True,
) -> dict:
    """Curate the frozen ledger and write the new run outputs.

    Returns the full payload (including the in-memory ``curated_rows`` and
    ``derivation_rows``); the JSON written to disk omits those two keys.
    """
    generated_at = generated_at or utc_now()
    output = Path(output_dir)
    if output.exists():
        if not output.is_dir():
            raise ValueError(f"output dir is not a directory: {output}")
        existing = sorted(entry.name for entry in output.iterdir())
        if existing:
            raise ValueError(
                f"output dir is not empty (refusing to overwrite): {output} -> {existing}"
            )

    ledger_path = Path(ledger)
    provenance_path = Path(provenance)
    definitions_path = Path(family_definitions)
    fasta_path = Path(reference_fasta)
    evidence_config = Path(evidence_columns) if evidence_columns is not None else EVIDENCE_COLUMNS_CONFIG
    validator = load_validator(validator_path)

    declared_columns = load_evidence_columns(evidence_config)
    if tuple(declared_columns) != tuple(getattr(validator, "EVIDENCE_COLUMNS", EVIDENCE_COLUMNS_FALLBACK)):
        raise ValueError(
            "evidence columns config disagrees with the validator's EVIDENCE_COLUMNS: "
            f"{tuple(declared_columns)}"
        )

    inputs = {
        "frozen_reference_ledger": ledger_path,
        "positive_primary_source_amendment": provenance_path,
        "family_definitions": definitions_path,
        "reference_fasta": fasta_path,
        "evidence_columns_config": evidence_config,
    }
    input_records = {
        name: {
            "path": str(path.resolve()),
            "size": path.stat().st_size if path.is_file() else None,
            "sha256": _sha256_file(path) if path.is_file() else None,
        }
        for name, path in inputs.items()
    }
    frozen_sha_before = input_records["frozen_reference_ledger"]["sha256"]

    ledger_fields, ledger_rows = read_tsv(ledger_path, "frozen reference ledger")
    provenance_fields, provenance_rows = read_tsv(
        provenance_path, "positive primary-source amendment"
    )
    read_tsv(definitions_path, "family definitions")
    sequences = read_fasta(fasta_path)

    if not ledger_rows:
        raise ValueError(f"frozen reference ledger has no data rows: {ledger_path}")
    missing_fasta = sorted(
        _text(row.get("reference_id")) for row in ledger_rows
        if _text(row.get("reference_id")) not in sequences
    )
    if missing_fasta:
        raise ValueError("reference FASTA has no record for ledger rows: " + ",".join(missing_fasta[:10]))

    positives = [row for row in ledger_rows if _text(row.get("evidence_status")) == EVIDENCE_STATUS_POSITIVE]
    if expected_rows is not None and len(ledger_rows) != expected_rows:
        raise ValueError(f"frozen ledger has {len(ledger_rows)} rows; expected {expected_rows}")
    if expected_positives is not None and len(positives) != expected_positives:
        raise ValueError(
            f"frozen ledger has {len(positives)} experimental positives; expected {expected_positives}"
        )

    curation = build_curation(
        ledger_fields, ledger_rows, provenance_fields, provenance_rows, sequences, declared_columns
    )
    curated_rows = curation["curated_rows"]
    derivation_rows = curation["derivation_rows"]

    output.mkdir(parents=True, exist_ok=True)
    curated_path = output / "phaded_reference_ledger_curated_v2.tsv"
    curated_fields = list(ledger_fields) + [
        column for column in declared_columns if column not in ledger_fields
    ]
    with curated_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=curated_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in curated_rows:
            writer.writerow({field: row.get(field, "") for field in curated_fields})

    # The validator is authoritative for the independence/gate numbers.
    gate_report = validator.validate_ledger(
        curated_path,
        definitions_path,
        fasta_path,
        evidence_column_mode="strict",
        expected_family_count=expected_family_count,
    )
    independence_summary = validator.summarize_independence(curated_rows)
    # Fail-closed self-check: the rows this script marked as counting must agree with
    # the validator's own independence accounting (unique independence group).
    counting = [row for row in derivation_rows if row["counts_as_functional_positive"] == "true"]
    counting_groups = {row["independence_group"] for row in counting}
    if (
        len(counting_groups) != independence_summary["functional_positive_count"]
        or len(counting) != independence_summary["functional_positive_row_count"]
    ):
        raise ValueError(
            "internal counting disagreement with "
            "validate_phaded_reference_ledger.summarize_independence: "
            f"rows={len(counting)} vs {independence_summary['functional_positive_row_count']}, "
            f"groups={len(counting_groups)} vs {independence_summary['functional_positive_count']}"
        )

    family_gate = gate_report["functional_calibration_gate"]
    family_row_counts = Counter(_text(row.get("phaded_family_id")) for row in curated_rows)
    family_positive_row_counts = Counter(
        _text(row.get("phaded_family_id")) for row in curated_rows
        if _text(row.get("evidence_status")) == EVIDENCE_STATUS_POSITIVE
    )
    eligible_groups = sorted({
        row["independence_group"] for row in curated_rows
        if _text(row.get("functional_calibration_eligible")) == "true"
    })
    counts = {
        "ledger_rows": len(curated_rows),
        "annotation_only_rows": sum(
            1 for row in curated_rows if _text(row.get("evidence_status")) == EVIDENCE_STATUS_ANNOTATION
        ),
        "experimental_positive_rows": len(positives),
        "positives_by_grade": dict(sorted(Counter(
            _text(row.get("experimental_evidence_grade")) for row in curated_rows
            if _text(row.get("evidence_status")) == EVIDENCE_STATUS_POSITIVE
        ).items())),
        "annotation_rows_by_grade": dict(sorted(Counter(
            _text(row.get("experimental_evidence_grade")) for row in curated_rows
            if _text(row.get("evidence_status")) == EVIDENCE_STATUS_ANNOTATION
        ).items())),
        "functional_positive_rows": independence_summary["functional_positive_row_count"],
        "functional_positive_groups": independence_summary["functional_positive_count"],
        "functional_calibration_eligible_rows": sum(
            1 for row in curated_rows if _text(row.get("functional_calibration_eligible")) == "true"
        ),
        "functional_calibration_eligible_groups": len(eligible_groups),
        "discovery_training_eligible_rows": sum(
            1 for row in curated_rows if _text(row.get("discovery_training_eligible")) == "true"
        ),
        "sequence_integrity_counts": dict(sorted(Counter(
            _text(row.get("sequence_integrity")) for row in curated_rows
        ).items())),
        "independence_group_values": len({
            row["independence_group"] for row in derivation_rows
        }),
        "experiment_unit_values": len({row["experiment_unit_id"] for row in derivation_rows}),
    }

    strict_gate = {
        "validator": str(Path(validator_path or VALIDATOR_PATH).resolve()),
        "evidence_column_mode": gate_report["evidence_column_mode"],
        "status": gate_report["status"],
        "row_count": gate_report["row_count"],
        "fasta_count": gate_report["fasta_count"],
        "evidence_status_counts": gate_report["evidence_status_counts"],
        "pending_evidence_columns": gate_report["pending_evidence_columns"],
        "evidence_columns_declared": gate_report["evidence_columns_declared"],
        "functional_positive_count": gate_report["functional_positive_count"],
        "functional_positive_row_count": gate_report["functional_positive_row_count"],
        "family_functional_positive_counts": family_gate["family_functional_positive_counts"],
        "family_distinct_genus_counts": family_gate["family_distinct_genus_counts"],
        "families_meeting_minimum": family_gate["families_meeting_minimum"],
        "families_below_minimum": family_gate["families_below_minimum"],
        "minimum_independent_positive_count": family_gate["minimum_independent_positive_count"],
        "positive_count_basis": family_gate["positive_count_basis"],
        "source_file_sha256": gate_report["source_file_sha256"],
        "evidence_columns_config_sha256": gate_report["evidence_columns_config_sha256"],
    }

    payload: dict[str, Any] = {
        "run_id": run_id,
        "generated_at": generated_at,
        "curation_schema_version": "F13-v1",
        "positive_count_basis": "unique_independence_group",
        "minimum_independent_positive_count": family_gate["minimum_independent_positive_count"],
        "independence_summary": independence_summary,
        "family_functional_positive_counts": family_gate["family_functional_positive_counts"],
        "family_distinct_genus_counts": family_gate["family_distinct_genus_counts"],
        "family_row_counts": dict(sorted(family_row_counts.items())),
        "family_positive_row_counts": dict(sorted(family_positive_row_counts.items())),
        "families_meeting_minimum_independent_positives": family_gate["families_meeting_minimum"],
        "families_below_minimum_independent_positives": family_gate["families_below_minimum"],
        "counts": counts,
        "pending_value_classes": curation["pending_value_classes"],
        "positives_with_pending_values": curation["positives_with_pending_values"],
        "strict_gate": strict_gate,
        "inputs": input_records,
        "outputs": {
            "curated_ledger": {
                "path": str(curated_path.resolve()),
                "sha256": _sha256_file(curated_path),
                "row_count": len(curated_rows),
                "column_count": len(curated_fields),
            }
        },
        "curated_rows": curated_rows,
        "derivation_rows": derivation_rows,
    }

    summary_path = output / "reference_evidence_independence_summary.json"
    on_disk = {key: value for key, value in payload.items() if key not in {"curated_rows", "derivation_rows"}}
    summary_path.write_text(
        json.dumps(on_disk, indent=2, ensure_ascii=False, sort_keys=False) + "\n", encoding="utf-8"
    )

    if derivation_table is not None:
        table_path = Path(derivation_table)
        table_path.parent.mkdir(parents=True, exist_ok=True)
        derivation_fields = [
            "reference_id", "accession", "accession_version", "accession_version_basis", "gene_id",
            "seed_gi", "seed_gi_source", "evidence_type", "confidence", "primary_pmid", "primary_doi",
            "study_key", "study_id", "study_id_basis", "study_basis", "experiment_unit_id",
            "independence_group", "independence_binding", "independence_members", "independence_basis",
            "independence_documented_bindings", "grade", "grade_basis", "assay_directness",
            "experimental_substrate_class", "experimental_substrate_class_basis",
            "experimental_system", "negative_control_description", "catalytic_mutant_evidence",
            "localization_evidence", "sequence_integrity", "discovery_training_eligible",
            "functional_calibration_eligible", "counts_as_functional_positive", "amendment_note",
        ]
        with table_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=derivation_fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in derivation_rows:
                writer.writerow({field: row.get(field, "") for field in derivation_fields})
        payload["outputs"]["derivation_table"] = {
            "path": str(table_path.resolve()),
            "sha256": _sha256_file(table_path),
            "row_count": len(derivation_rows),
        }

    if write_report and report is not None:
        # Written before the input contract so that every produced artifact is hashed there.
        report_path = Path(report)
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(
            build_report(
                run_id=run_id, generated_at=generated_at, payload=payload, inputs=input_records
            ),
            encoding="utf-8",
        )
        payload["outputs"]["curation_report"] = {
            "path": str(report_path.resolve()),
            "sha256": _sha256_file(report_path),
        }

    if input_contract is not None:
        contract_path = Path(input_contract)
        contract_path.parent.mkdir(parents=True, exist_ok=True)
        frozen_sha_after = _sha256_file(ledger_path)
        contract = {
            "schema_version": "1.0",
            "run_id": run_id,
            "generated_at": generated_at,
            "status": "candidate_only_reference_curation",
            "purpose": (
                "F13: grade the frozen 723-row PhaDED reference ledger with the 15 v2 evidence "
                "columns and count functional positives by unique independence group"
            ),
            "authorization": {
                "candidate_only_execution": True,
                "frozen_evidence_modified": False,
                "formal_scan_authorized": False,
                "formal_registry_modified": False,
                "candidate_deleted": False,
            },
            "inputs": input_records,
            "frozen_ledger_sha256_before": frozen_sha_before,
            "frozen_ledger_sha256_after": frozen_sha_after,
            "frozen_ledger_unchanged": frozen_sha_before == frozen_sha_after,
            "outputs": payload["outputs"],
            "environment": {"python_version": sys.version.split()[0]},
            "phenotype_boundary": (
                "Evidence grades recorded here are reference-panel curation. No grade asserts a "
                "verified PHB/PHA degradation phenotype for any GTDB candidate."
            ),
        }
        contract_path.write_text(
            json.dumps(contract, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
        )
        payload["input_contract"] = {
            "path": str(contract_path.resolve()),
            "sha256": _sha256_file(contract_path),
        }
        payload["frozen_ledger_sha256_before"] = frozen_sha_before
        payload["frozen_ledger_sha256_after"] = frozen_sha_after
        if frozen_sha_before != frozen_sha_after:
            raise ValueError(
                "the frozen reference ledger changed during the run (SHA-256 "
                f"{frozen_sha_before} -> {frozen_sha_after}); refusing to continue"
            )

    payload["console_summary"] = {
        "run_id": run_id,
        "status": "ok",
        "curated_ledger": payload["outputs"]["curated_ledger"],
        "counts": counts,
        "family_functional_positive_counts": payload["family_functional_positive_counts"],
        "family_distinct_genus_counts": payload["family_distinct_genus_counts"],
        "families_meeting_minimum_independent_positives": payload[
            "families_meeting_minimum_independent_positives"
        ],
        "strict_gate": {
            "status": strict_gate["status"],
            "evidence_column_mode": strict_gate["evidence_column_mode"],
            "pending_evidence_columns": strict_gate["pending_evidence_columns"],
            "functional_positive_count": strict_gate["functional_positive_count"],
        },
        "pending_value_classes": curation["pending_value_classes"],
        "positives_with_pending_values": curation["positives_with_pending_values"],
        "outputs": payload["outputs"],
    }
    return payload


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--ledger", required=True, help="frozen 723-row reference ledger TSV (read-only)")
    parser.add_argument("--provenance", required=True, help="30-row primary-source amendment TSV")
    parser.add_argument("--family-definitions", required=True, help="family definition TSV (validator input)")
    parser.add_argument("--reference-fasta", required=True, help="frozen reference FASTA (validator input)")
    parser.add_argument("--output-dir", required=True, help="empty/new directory for the curated outputs")
    parser.add_argument("--evidence-columns", default=None, help=f"evidence columns contract (default: {EVIDENCE_COLUMNS_CONFIG})")
    parser.add_argument("--derivation-table", default=None, help="write the auditable per-positive derivation TSV here")
    parser.add_argument("--input-contract", default=None, help="write the run input contract JSON here")
    parser.add_argument("--report", default=None, help="write the human-readable curation report here")
    parser.add_argument("--no-report", action="store_true", help="skip the markdown curation report")
    parser.add_argument("--validator", default=None, help=f"validator module (default: {VALIDATOR_PATH})")
    parser.add_argument("--expected-rows", type=int, default=None, help="refuse a ledger with a different row count")
    parser.add_argument("--expected-positives", type=int, default=None, help="refuse a ledger with a different positive count")
    parser.add_argument("--expected-family-count", type=int, default=38, help="family definitions count (validator input)")
    parser.add_argument("--run-id", default="", help="run identifier recorded in the outputs")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        payload = run(
            ledger=args.ledger,
            provenance=args.provenance,
            family_definitions=args.family_definitions,
            reference_fasta=args.reference_fasta,
            output_dir=args.output_dir,
            evidence_columns=args.evidence_columns,
            derivation_table=args.derivation_table,
            input_contract=args.input_contract,
            report=args.report,
            validator_path=args.validator,
            expected_rows=args.expected_rows,
            expected_positives=args.expected_positives,
            expected_family_count=args.expected_family_count,
            run_id=args.run_id,
            write_report=not args.no_report,
        )
    except ValueError as error:
        print(f"error: {error}", file=sys.stderr)
        return 2
    print(json.dumps(payload["console_summary"], indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
