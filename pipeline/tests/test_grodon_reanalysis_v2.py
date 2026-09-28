"""Tests for the gRodon re-analysis v2 correction module (plan Task 12).

Written before the implementation (test-first).  ``grodon_reanalysis_v2.py``
supersedes the two legacy scripts of
``deploy/20260920_phaded_grodon_growth_01/`` for the 2026-09-28 evidence-model
redesign.  The corrections under test are exactly the six points of the design
document's "下游统计" section:

1. group names are "candidate gene carrier" vs "candidate not detected under
   the defined search/quality conditions" (never a phenotype label), and the
   final positive set and the control set are asserted to be **disjoint**;
2. ``growth_rate_per_h`` is a verified conversion of the gRodon minimum
   doubling time (ln2 / d), applied exactly once, with an explicit resolver so a
   table that already carries a growth rate never gets a second conversion;
3. the **genus-level mean difference is the primary estimand**; bootstrap and
   permutation resample genera (clusters), never individual genomes; median,
   sign test, per-genome and compartment comparisons are secondary and
   labelled as such;
4. reused predictions enter the primary analysis only when tool, model
   settings, temperature handling, CDS preparation, ribosomal-marker procedure
   and output formula are all comparable; otherwise they are retained in a
   sensitivity table with a named reason;
5. selection coverage is quantified by named reason buckets whose accession
   counts are computed from sets, never by subtracting aggregates;
6. the 4,507 -> 4,441 = 66 difference is explained by reason buckets summing to
   exactly 66.

Every fixture here is synthetic.  Nothing in this module reads or writes
``runs/``, ``results/``, ``deploy/``, or executes gRodon/R.
"""

import importlib.util
import math
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "grodon_reanalysis_v2.py"
LEGACY_DEPLOY = (
    Path(__file__).resolve().parents[2] / "deploy" / "20260920_phaded_grodon_growth_01"
)
NEW_DEPLOY = (
    Path(__file__).resolve().parents[2] / "deploy" / "20260928_phaded_grodon_reanalysis_v2_01"
)

LN2 = math.log(2.0)


def load_module():
    spec = importlib.util.spec_from_file_location("grodon_reanalysis_v2", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def gid(i):
    """9-digit GTDB accession, e.g. GCA_000000001.1."""
    return f"GCA_000000{i:03d}.1"


def tax(
    genus="TestGenus",
    domain="Bacteria",
    phylum="Pseudomonadota",
    cls="Gammaproteobacteria",
    order="TestOrder",
    family="TestFamily",
    species="sp",
):
    return {
        "domain": domain,
        "phylum": phylum,
        "class": cls,
        "order": order,
        "family": family,
        "genus": genus,
        "species": species,
    }


def carrier_row(**extra):
    """Synthetic carrier (candidate-positive) entry for one genome.

    Deliberately carries no ``genome_id``: the key of the ``carriers`` mapping
    is the genome identifier, exactly as in the v2 carrier table.
    """
    row = {
        "candidate_family": "DED_hfam_70",
        "candidate_count": 4,
        "group": "intracellular",
    }
    row.update(extra)
    return row


def genome(genome_id, genus, status, growth, **extra):
    """One genus-balanced analysis row in the module's own column shape."""
    row = {
        "genome_id": genome_id,
        "genus": genus,
        "candidate_detection_status": status,
        "group": extra.pop("group", "intracellular"),
        "growth_rate_per_h": growth,
    }
    row.update(extra)
    return row


CARRIER = "candidate_gene_carrier"
CONTROL = "candidate_not_detected_under_defined_search"


def valid_manifest():
    """A minimal manifest with one carrier genus and one not-detected-only genus.

    ``exclusion`` is the pool-exclusion ledger and therefore contains every
    final positive, as the real v2 carrier table does.
    """
    module = load_module()
    carriers = {
        gid(1): carrier_row(),
        gid(2): carrier_row(group="extracellular"),
    }
    tax_by_gid = {
        gid(1): tax("GenA"),
        gid(2): tax("GenB"),
        gid(3): tax("GenA"),
        gid(4): tax("GenB"),
    }
    return module.build_manifest(
        carriers=carriers,
        exclusion=set(),
        tax_by_gid=tax_by_gid,
        seed=42,
        max_per_genus=1,
        file_exists=lambda _: True,
    )


def balanced_frame():
    """Two genera x two carriers/controls with a known, strong effect.

    G1: carriers 0.9/0.7, controls 0.2/0.4 -> genus delta +0.5
    G2: carriers 0.8/0.6, controls 0.1/0.3 -> genus delta +0.5
    """
    rows = [
        genome("G1_C1", "G1", CARRIER, 0.9),
        genome("G1_C2", "G1", CARRIER, 0.7),
        genome("G1_N1", "G1", CONTROL, 0.2),
        genome("G1_N2", "G1", CONTROL, 0.4),
        genome("G2_C1", "G2", CARRIER, 0.8),
        genome("G2_C2", "G2", CARRIER, 0.6),
        genome("G2_N1", "G2", CONTROL, 0.1),
        genome("G2_N2", "G2", CONTROL, 0.3),
    ]
    return pd.DataFrame(rows)


def null_frame():
    """Two genera whose deltas cancel exactly -> mean delta 0, CI must contain 0."""
    return pd.DataFrame(
        [
            genome("H1_C1", "H1", CARRIER, 0.80),
            genome("H1_N1", "H1", CONTROL, 0.30),
            genome("H2_C1", "H2", CARRIER, 0.10),
            genome("H2_N1", "H2", CONTROL, 0.60),
        ]
    )


def reference_profile(**extra):
    profile = {
        "tool": "gRodon",
        "tool_version": "2.4.0",
        "model_settings": "gRodon2_default_predictGrowth",
        "temperature_handling": "temperature_celsius=37;not_used_for_prediction",
        "cds_preparation": "pyrodigal_meta_v3.0.0_min_protein_30nt90",
        "ribosomal_marker_procedure": "pfam_ribosomal_hmm_cut_ga_min10",
        "output_formula": "growth_rate_per_h=ln2/doubling_time_h",
    }
    profile.update(extra)
    return profile


def reused_record(genome_id, **extra):
    row = {
        "genome_id": genome_id,
        "tool": "gRodon",
        "tool_version": "2.4.0",
        "model_settings": "gRodon2_default_predictGrowth",
        "temperature_handling": "temperature_celsius=37;not_used_for_prediction",
        "cds_preparation": "pyrodigal_meta_v3.0.0_min_protein_30nt90",
        "ribosomal_marker_procedure": "pfam_ribosomal_hmm_cut_ga_min10",
        "output_formula": "growth_rate_per_h=ln2/doubling_time_h",
        "doubling_time_h": 2.0,
        "growth_rate_per_h": LN2 / 2.0,
        "status": "ok",
    }
    row.update(extra)
    return row


class ManifestOverlapTests(unittest.TestCase):
    """Requirement 1: zero overlap between the final positive and control sets."""

    def test_final_positive_and_control_sets_are_disjoint(self):
        with self.assertRaisesRegex(ValueError, "positive/control overlap"):
            module = load_module()
            module.build_manifest(
                carriers={"G1": carrier_row()},
                exclusion=set(),
                tax_by_gid={"G1": tax(), "G2": tax()},
                seed=42,
                max_per_genus=0,
                file_exists=lambda _: True,
                control_pool={"G1"},
            )

    def test_overlap_guard_fires_when_a_carrier_is_declared_excluded(self):
        module = load_module()
        carriers = {gid(1): carrier_row(), gid(2): carrier_row()}
        tax_by_gid = {gid(i): tax() for i in range(1, 6)}
        with self.assertRaisesRegex(ValueError, "positive/control overlap"):
            module.build_manifest(
                carriers=carriers,
                exclusion={gid(2)},
                tax_by_gid=tax_by_gid,
                seed=42,
                max_per_genus=0,
                file_exists=lambda _: True,
            )

    def test_overlap_guard_fires_when_the_declared_control_pool_has_a_carrier(self):
        module = load_module()
        carriers = {gid(1): carrier_row(), gid(2): carrier_row()}
        tax_by_gid = {gid(i): tax() for i in range(1, 6)}
        with self.assertRaisesRegex(ValueError, "positive/control overlap"):
            module.build_manifest(
                carriers=carriers,
                exclusion=set(),
                tax_by_gid=tax_by_gid,
                seed=42,
                max_per_genus=2,
                file_exists=lambda _: True,
                control_pool={gid(1), gid(3), gid(4)},
            )

    def test_a_manifest_built_from_the_same_inputs_has_no_overlap(self):
        """The guard is narrow: the same inputs without contamination succeed."""
        module = load_module()
        carriers = {gid(1): carrier_row(), gid(2): carrier_row()}
        tax_by_gid = {gid(i): tax() for i in range(1, 6)}
        rows, stats = module.build_manifest(
            carriers=carriers,
            exclusion=set(),
            tax_by_gid=tax_by_gid,
            seed=42,
            max_per_genus=2,
            file_exists=lambda _: True,
            control_pool={gid(3), gid(4), gid(5)},
        )
        self.assertTrue(rows)
        self.assertEqual(module.carrier_genomes(rows) & module.control_genomes(rows), set())
        self.assertTrue(module.sets_are_disjoint(rows))
        self.assertEqual(stats["overlap_detected"], False)

    def test_max_per_genus_zero_means_no_pair_limit(self):
        module = load_module()
        carriers = {gid(i): carrier_row() for i in range(1, 4)}
        tax_by_gid = {gid(i): tax("GenA") for i in range(1, 10)}
        rows, stats = module.build_manifest(
            carriers=carriers,
            exclusion=set(),
            tax_by_gid=tax_by_gid,
            seed=42,
            max_per_genus=0,
            file_exists=lambda _: True,
        )
        # 3 carriers, 6 usable controls, no cap -> all 3 carriers are paired.
        self.assertEqual(stats["matched_pairs"], 3)
        self.assertEqual(stats["manifest_carrier"], 3)
        self.assertEqual(stats["manifest_control"], 3)
        self.assertEqual(stats["regime"], "matched_same_genus_1to1")
        self.assertTrue(module.sets_are_disjoint(rows))

    def test_unpaired_carriers_contribute_no_manifest_row(self):
        module = load_module()
        carriers = {gid(i): carrier_row() for i in range(1, 4)}
        tax_by_gid = {gid(i): tax("GenA") for i in range(1, 10)}
        rows, stats = module.build_manifest(
            carriers=carriers,
            exclusion=set(),
            tax_by_gid=tax_by_gid,
            seed=42,
            max_per_genus=1,
            file_exists=lambda _: True,
        )
        self.assertEqual(stats["manifest_rows"], 2)
        self.assertEqual(len(module.carrier_genomes(rows)), 1)
        self.assertEqual(stats["carrier_genomes_without_a_same_genus_control"], 2)
        self.assertTrue(module.sets_are_disjoint(rows))

    def test_overlap_guard_fires_when_a_carrier_is_declared_excluded(self):
        module = load_module()
        carriers = {gid(1): carrier_row(), gid(2): carrier_row()}
        tax_by_gid = {gid(i): tax() for i in range(1, 6)}
        with self.assertRaisesRegex(ValueError, "positive/control overlap"):
            module.build_manifest(
                carriers=carriers,
                exclusion={gid(2)},
                tax_by_gid=tax_by_gid,
                seed=42,
                max_per_genus=0,
                file_exists=lambda _: True,
            )

    def test_overlap_guard_fires_when_a_genome_is_both_carrier_and_excluded(self):
        module = load_module()
        carriers = {gid(1): carrier_row(), gid(2): carrier_row()}
        tax_by_gid = {gid(i): tax() for i in range(1, 6)}
        with self.assertRaisesRegex(ValueError, "positive/control overlap"):
            module.build_manifest(
                carriers=carriers,
                exclusion={gid(1)},
                tax_by_gid=tax_by_gid,
                seed=42,
                max_per_genus=2,
                file_exists=lambda _: True,
            )

    def test_control_pool_never_contains_a_carrier_or_an_excluded_genome(self):
        module = load_module()
        carriers = {gid(1): carrier_row()}
        tax_by_gid = {gid(i): tax("GenA") for i in range(1, 7)}
        rows, stats = module.build_manifest(
            carriers=carriers,
            exclusion={gid(3), gid(4)},
            tax_by_gid=tax_by_gid,
            seed=42,
            max_per_genus=0,
            file_exists=lambda _: True,
        )
        controls = module.control_genomes(rows)
        self.assertEqual(controls & set(carriers), set())
        self.assertEqual(controls & {gid(3), gid(4)}, set())
        self.assertTrue(controls <= {gid(2), gid(5), gid(6)})
        self.assertEqual(stats["manifest_carrier"], 1)
        self.assertEqual(stats["manifest_control"], 1)
        self.assertEqual(stats["matched_pairs"], 1)
        self.assertEqual(stats["exclusion_genomes"], 2)
        self.assertEqual(stats["control_pool_genomes"], 3)
        self.assertTrue(module.sets_are_disjoint(rows))

    def test_manifest_is_balanced_one_to_one_within_genus(self):
        module = load_module()
        carriers = {gid(i): carrier_row() for i in range(1, 4)}
        tax_by_gid = {gid(i): tax("GenA") for i in range(1, 9)}
        rows, stats = module.build_manifest(
            carriers=carriers,
            exclusion=set(),
            tax_by_gid=tax_by_gid,
            seed=42,
            max_per_genus=3,
            file_exists=lambda _: True,
        )
        self.assertEqual(stats["manifest_carrier"], 3)
        self.assertEqual(stats["manifest_control"], 3)
        self.assertEqual(stats["matched_pairs"], 3)
        self.assertEqual(stats["genera_with_a_control"], 1)
        self.assertEqual({row["genus"] for row in rows}, {"GenA"})
        pairs = {row["pair_id"] for row in rows if row["pair_id"]}
        self.assertEqual(len(pairs), 3)

    def test_carrier_without_a_same_genus_control_has_no_control_row(self):
        module = load_module()
        carriers = {gid(1): carrier_row(), gid(2): carrier_row()}
        tax_by_gid = {gid(1): tax("GenA"), gid(2): tax("Lonely"), gid(3): tax("GenA")}
        rows, stats = module.build_manifest(
            carriers=carriers,
            exclusion=set(),
            tax_by_gid=tax_by_gid,
            seed=42,
            max_per_genus=1,
            file_exists=lambda _: True,
        )
        self.assertEqual(module.carrier_genomes(rows), {gid(1)})
        self.assertEqual(module.control_genomes(rows), {gid(3)})
        self.assertTrue(module.sets_are_disjoint(rows))
        self.assertEqual(stats["carrier_genomes_without_a_same_genus_control"], 1)
        self.assertEqual(stats["genera_skipped_no_control"], 1)
        self.assertEqual(stats["manifest_carrier"], 1)
        self.assertEqual(stats["manifest_control"], 1)

    def test_manifest_seed_determinism(self):
        module = load_module()
        carriers = {gid(i): carrier_row() for i in range(1, 5)}
        tax_by_gid = {gid(i): tax("GenA") for i in range(1, 13)}
        first, _ = module.build_manifest(
            carriers=carriers, exclusion=set(), tax_by_gid=tax_by_gid, seed=42,
            max_per_genus=2, file_exists=lambda _: True,
        )
        second, _ = module.build_manifest(
            carriers=carriers, exclusion=set(), tax_by_gid=tax_by_gid, seed=42,
            max_per_genus=2, file_exists=lambda _: True,
        )
        self.assertEqual(
            [row["genome_id"] for row in first], [row["genome_id"] for row in second]
        )

    def test_genomes_without_a_fasta_are_not_selected(self):
        module = load_module()
        carriers = {gid(1): carrier_row()}
        tax_by_gid = {gid(i): tax("GenA") for i in range(1, 6)}
        missing = {gid(2)}
        rows, stats = module.build_manifest(
            carriers=carriers,
            exclusion=set(),
            tax_by_gid=tax_by_gid,
            seed=42,
            max_per_genus=1,
            file_exists=lambda path: Path(path).name.replace(
                "_genomic.fna.gz", ""
            ) not in missing,
        )
        self.assertEqual(module.control_genomes(rows), {gid(5)})
        self.assertEqual(stats["control_genomes_missing_fasta"], 1)

    def test_disjointness_helper_detects_an_overlap(self):
        module = load_module()
        rows = [
            {"genome_id": "A", "candidate_detection_status": CARRIER},
            {"genome_id": "A", "candidate_detection_status": CONTROL},
        ]
        self.assertFalse(module.sets_are_disjoint(rows))
        with self.assertRaisesRegex(ValueError, "positive/control overlap"):
            module.assert_sets_are_disjoint(rows)

    def test_unknown_max_per_genus_is_rejected(self):
        module = load_module()
        with self.assertRaises(ValueError):
            module.build_manifest(
                carriers={"G1": carrier_row()},
                exclusion={"G1"},
                tax_by_gid={"G1": tax()},
                seed=42,
                max_per_genus=-1,
                file_exists=lambda _: True,
            )


class ManifestLabelTests(unittest.TestCase):
    """Requirement 1: candidate-carrier wording, never a phenotype label."""

    def test_manifest_uses_candidate_carrier_language(self):
        module = load_module()
        rows, _ = valid_manifest()
        self.assertEqual(
            {row["candidate_detection_status"] for row in rows},
            {"candidate_gene_carrier", "candidate_not_detected_under_defined_search"},
        )

    def test_control_rows_are_explicitly_not_detected_not_negative(self):
        module = load_module()
        rows, _ = valid_manifest()
        controls = [row for row in rows if row["candidate_detection_status"] == CONTROL]
        self.assertTrue(controls)
        for row in controls:
            self.assertEqual(row["group"], "control")
            self.assertEqual(row["candidate_count"], 0)
            self.assertEqual(row["candidate_family"], "")
            self.assertEqual(row["control_basis"], "same_genus_no_candidate_detected")

    def test_legacy_phenotype_labels_are_not_emitted(self):
        module = load_module()
        rows, _ = valid_manifest()
        blob = " ".join(
            str(value) for row in rows for value in row.values()
        ).lower()
        for forbidden in ("degrader", "non-degrader", "non_degrader", "positive_control"):
            self.assertNotIn(forbidden, blob)

    def test_forbidden_legacy_labels_constant_is_exposed(self):
        module = load_module()
        self.assertIn("degrader", module.FORBIDDEN_LEGACY_LABELS)
        self.assertIn("non-degrader", module.FORBIDDEN_LEGACY_LABELS)
        self.assertEqual(module.CANDIDATE_CARRIER, "candidate_gene_carrier")
        self.assertEqual(
            module.CANDIDATE_NOT_DETECTED,
            "candidate_not_detected_under_defined_search",
        )

    def test_manifest_columns_are_frozen(self):
        module = load_module()
        rows, _ = valid_manifest()
        for row in rows:
            self.assertEqual(list(row.keys()), list(module.MANIFEST_COLUMNS))


class GrowthRateUnitTests(unittest.TestCase):
    """Requirement 2: one tested conversion, applied exactly once."""

    def test_growth_rate_matches_ln2_over_doubling_time(self):
        module = load_module()
        self.assertAlmostEqual(module.growth_rate_per_h(1.0), LN2)
        self.assertAlmostEqual(module.growth_rate_per_h(0.5), 2.0 * LN2)
        self.assertAlmostEqual(module.growth_rate_per_h(2.0), LN2 / 2.0)
        self.assertAlmostEqual(module.growth_rate_per_h(1.0), 0.693147, places=6)

    def test_non_positive_doubling_time_raises(self):
        module = load_module()
        for bad in (0.0, -1.0, -0.001):
            with self.assertRaisesRegex(ValueError, "doubling time must be positive"):
                module.growth_rate_per_h(bad)

    def test_resolver_converts_a_doubling_time_field(self):
        module = load_module()
        record = {"doubling_time_h": 1.0, "status": "ok"}
        self.assertAlmostEqual(module.resolve_growth_rate(record), LN2)
        self.assertEqual(module.resolve_growth_rate_source(record), "doubling_time_h")

    def test_resolver_accepts_the_gr_d_family_field_name(self):
        module = load_module()
        record = {"d": 0.5, "LowerCI": 0.4, "UpperCI": 0.6}
        self.assertAlmostEqual(module.resolve_growth_rate(record), 2.0 * LN2)

    def test_resolver_does_not_convert_an_already_rate_field(self):
        module = load_module()
        record = {"growth_rate_per_h": 0.25}
        self.assertAlmostEqual(module.resolve_growth_rate(record), 0.25)
        self.assertEqual(module.resolve_growth_rate_source(record), "growth_rate_per_h")
        # The value is returned unchanged: a second conversion would give ln2/0.25.
        self.assertNotAlmostEqual(module.resolve_growth_rate(record), LN2 / 0.25)

    def test_resolver_applies_the_conversion_exactly_once_on_the_doubling_path(self):
        module = load_module()
        record = {"doubling_time_h": 2.0}
        self.assertAlmostEqual(module.resolve_growth_rate(record), LN2 / 2.0)
        self.assertNotAlmostEqual(module.resolve_growth_rate(record), LN2 / (LN2 / 2.0))

    def test_resolver_checks_unit_when_present(self):
        module = load_module()
        for unit in ("1/h", "per_h", "h^-1", "h-1", "hour^-1"):
            self.assertAlmostEqual(
                module.resolve_growth_rate({"growth_rate_per_h": 0.25, "rate_unit": unit}),
                0.25,
            )
        with self.assertRaisesRegex(ValueError, "rate_unit"):
            module.resolve_growth_rate({"growth_rate_per_h": 0.25, "rate_unit": "per_day"})

    def test_resolver_raises_on_a_missing_field(self):
        module = load_module()
        with self.assertRaisesRegex(ValueError, "growth rate"):
            module.resolve_growth_rate({"status": "ok", "genome_id": "G1"})

    def test_resolver_raises_on_both_fields_carrying_inconsistent_units(self):
        module = load_module()
        # 0.18 is not ln2/2.0 = 0.3466 -- the stored rate belongs to a different
        # doubling-time definition, so the record is not comparable.
        with self.assertRaisesRegex(ValueError, "output_formula"):
            module.resolve_growth_rate(
                {"doubling_time_h": 2.0, "growth_rate_per_h": 0.18}
            )
        # Rounding inside the tolerance is not a unit error and must pass.
        rounded = module.resolve_growth_rate(
            {"doubling_time_h": 2.0, "growth_rate_per_h": 0.3465}
        )
        self.assertAlmostEqual(rounded, 0.3465)

    def test_resolver_rejects_a_rate_that_looks_like_a_second_conversion(self):
        module = load_module()
        # ln2/(ln2/2) = 2.0: the table converted an already-converted rate.
        with self.assertRaisesRegex(ValueError, "output_formula"):
            module.resolve_growth_rate(
                {"doubling_time_h": 2.0, "growth_rate_per_h": LN2 / (LN2 / 2.0)}
            )

    def test_resolver_raises_on_non_positive_doubling_time(self):
        module = load_module()
        with self.assertRaisesRegex(ValueError, "doubling time must be positive"):
            module.resolve_growth_rate({"doubling_time_h": 0.0})

    def test_resolver_raises_on_a_non_finite_doubling_time(self):
        module = load_module()
        for bad in (float("nan"), float("inf")):
            with self.assertRaisesRegex(ValueError, "doubling time"):
                module.resolve_growth_rate({"doubling_time_h": bad})

    def test_accession_inputs_accept_a_series_a_frame_and_a_mapping(self):
        module = load_module()
        base = module.explain_manifest_difference(
            eligible_positive_accessions={"A", "B", "C"},
            manifest_positive_accessions={"A"},
            reason_buckets={"missing_fasta": {"B"}, "no_same_genus_control": {"C"}},
        )
        self.assertEqual(base["status"], "ok")
        for eligible, manifest, buckets in (
            (
                pd.Series(["A", "B", "C"]),
                pd.Series(["A"]),
                {"missing_fasta": ["B"], "no_same_genus_control": ["C"]},
            ),
            (
                pd.DataFrame({"accession": ["A", "B", "C"]}),
                pd.DataFrame({"genome_id": ["A"]}),
                {"missing_fasta": {"B"}, "no_same_genus_control": {"C"}},
            ),
            (
                {"A": 1, "B": 1, "C": 1},
                {"A": 1},
                {"missing_fasta": {"B"}, "no_same_genus_control": {"C"}},
            ),
        ):
            payload = module.explain_manifest_difference(
                eligible_positive_accessions=eligible,
                manifest_positive_accessions=manifest,
                reason_buckets=buckets,
            )
            self.assertEqual(payload["eligible_positive_count"], 3)
            self.assertEqual(payload["difference"], 2)
            self.assertEqual(payload["reason_bucket_counts"], base["reason_bucket_counts"])
            self.assertEqual(payload["status"], "ok")

    def test_a_bare_accession_string_is_rejected(self):
        module = load_module()
        with self.assertRaisesRegex(ValueError, "collection of accessions"):
            module.explain_manifest_difference(
                eligible_positive_accessions="ABC",
                manifest_positive_accessions={"A"},
                reason_buckets={"missing_fasta": {"B"}},
            )

    def test_resolve_growth_rate_series_matches_the_scalar_resolver(self):
        module = load_module()
        table = pd.DataFrame(
            {
                "genome_id": ["G1", "G2", "G3"],
                "doubling_time_h": [1.0, 2.0, 4.0],
            }
        )
        resolved = module.resolve_growth_rate_series(table)
        self.assertAlmostEqual(float(resolved.loc["G1"]), LN2)
        self.assertAlmostEqual(float(resolved.loc["G2"]), LN2 / 2.0)
        self.assertAlmostEqual(float(resolved.loc["G3"]), LN2 / 4.0)

    def test_unit_constant_is_per_hour(self):
        module = load_module()
        self.assertEqual(module.GROWTH_RATE_UNIT, "1/h")


class GenusDependenceTests(unittest.TestCase):
    """Requirement 3: the genus is the dependence unit for every resampling."""

    def test_cluster_bootstrap_resamples_genera_not_genomes(self):
        module = load_module()
        table = pd.DataFrame(
            [
                genome("A1", "G1", CARRIER, 0.9),
                genome("A2", "G1", CARRIER, 0.8),
                genome("A3", "G1", CARRIER, 0.7),
                genome("A4", "G1", CONTROL, 0.2),
                genome("A5", "G1", CONTROL, 0.1),
                genome("A6", "G1", CONTROL, 0.3),
                genome("B1", "G2", CARRIER, 0.6),
                genome("B2", "G2", CONTROL, 0.2),
            ]
        )
        result = module.cluster_bootstrap_mean_difference(
            table, np.random.default_rng(42), n_boot=50
        )
        self.assertEqual(result["resampling_unit"], "genus")
        self.assertEqual(result["n_genera"], 2)
        self.assertEqual(result["n_resamples"], 50)
        self.assertEqual(result["resample_size"], 2)
        # Every bootstrap draw is a genus label drawn from the genus blocks.
        self.assertTrue(set(result["resampled_units"]) <= {"G1", "G2"})
        self.assertEqual(len(result["resampled_units"]), result["n_resamples"] * 2)
        self.assertAlmostEqual(result["estimate"], 0.5, places=6)
        self.assertLessEqual(result["ci_low"], result["estimate"])
        self.assertLessEqual(result["estimate"], result["ci_high"])

    def test_cluster_bootstrap_requires_a_genus_column(self):
        module = load_module()
        table = pd.DataFrame(
            [{"genome_id": "A", "candidate_detection_status": CARRIER, "growth_rate_per_h": 0.5}]
        )
        with self.assertRaisesRegex(ValueError, "genus"):
            module.cluster_bootstrap_mean_difference(table, np.random.default_rng(1))

    def test_genus_permutation_resamples_genera(self):
        module = load_module()
        table = pd.DataFrame(
            [
                genome("A1", "G1", CARRIER, 0.9),
                genome("A2", "G1", CONTROL, 0.2),
                genome("B1", "G2", CARRIER, 0.8),
                genome("B2", "G2", CONTROL, 0.1),
                genome("C1", "G3", CARRIER, 0.7),
                genome("C2", "G3", CONTROL, 0.3),
            ]
        )
        result = module.genus_level_permutation_test(
            table, np.random.default_rng(7), n_permutations=200
        )
        self.assertEqual(result["resampling_unit"], "genus")
        self.assertEqual(result["n_genera"], 3)
        self.assertEqual(result["n_permutations"], 200)
        self.assertAlmostEqual(result["observed"], 0.6, places=6)
        self.assertGreaterEqual(result["p_value"], 0.0)
        self.assertLessEqual(result["p_value"], 1.0)
        self.assertTrue(set(result["resampled_units"]) <= {"G1", "G2", "G3"})

    def test_genus_mean_difference_is_the_primary_estimand(self):
        module = load_module()
        estimates = module.genus_level_estimates(
            balanced_frame(), np.random.default_rng(20260928), n_boot=400
        )
        self.assertEqual(estimates["primary_estimand"], "genus_level_mean_difference")
        self.assertEqual(estimates["resampling_unit"], "genus")
        self.assertEqual(estimates["n_genera"], 2)
        self.assertAlmostEqual(estimates["mean_delta"], 0.5, places=6)
        self.assertEqual(estimates["effect_unit"], module.GROWTH_RATE_UNIT)
        self.assertTrue(estimates["difference_detected"])
        self.assertGreater(estimates["ci_low"], 0.0)

    def test_secondary_analyses_are_labelled_secondary(self):
        module = load_module()
        estimates = module.genus_level_estimates(
            balanced_frame(), np.random.default_rng(20260928), n_boot=200
        )
        secondary = set(estimates["secondary_analyses"])
        self.assertIn("genus_median_difference", secondary)
        self.assertIn("genus_sign_test", secondary)
        self.assertIn("genome_level_mann_whitney_u", secondary)
        self.assertIn("intracellular_vs_extracellular", secondary)
        for name in secondary:
            self.assertNotEqual(name, estimates["primary_estimand"])
        self.assertIn("median_delta", estimates)
        self.assertIn("sign_test_p", estimates)
        self.assertIn("genome_level_mann_whitney_p", estimates)

    def test_no_difference_wording_requires_absence_of_a_preregistered_margin(self):
        module = load_module()
        estimates = module.genus_level_estimates(
            null_frame(), np.random.default_rng(20260928), n_boot=200
        )
        self.assertFalse(estimates["difference_detected"])
        self.assertEqual(estimates["conclusion_basis"], "not_detected")
        self.assertIn("no difference detected under the current design", estimates["statement"])
        self.assertNotIn("equivalence", estimates["statement"].lower())

    def test_equivalence_margin_enables_a_bounded_statement(self):
        module = load_module()
        estimates = module.genus_level_estimates(
            null_frame(),
            np.random.default_rng(20260928),
            n_boot=400,
            equivalence_margin=0.5,
        )
        self.assertEqual(estimates["conclusion_basis"], "within_equivalence_margin")
        self.assertEqual(estimates["equivalence_margin"], 0.5)
        self.assertIn("equivalence", estimates["statement"].lower())
        for token in ("degrader", "non-degrader"):
            self.assertNotIn(token, estimates["statement"].lower())

    def test_absence_of_a_margin_does_not_produce_an_equivalence_claim(self):
        module = load_module()
        statement = module.equivalence_statement(
            (0.1, 0.2), mean_delta=0.15, equivalence_margin=None
        )
        self.assertEqual(statement["basis"], "not_detected")
        self.assertNotIn("equivalent", statement["statement"].lower())

    def test_margin_not_containing_the_interval_does_not_claim_equivalence(self):
        module = load_module()
        statement = module.equivalence_statement(
            (0.10, 0.90), mean_delta=0.5, equivalence_margin=0.2
        )
        self.assertEqual(statement["basis"], "outside_equivalence_margin")
        self.assertIn("margin", statement["statement"].lower())

    def test_negative_equivalence_margin_is_rejected(self):
        module = load_module()
        with self.assertRaises(ValueError):
            module.genus_level_estimates(
                balanced_frame(), np.random.default_rng(1), n_boot=50,
                equivalence_margin=-0.1,
            )

    def test_genus_deltas_balance_one_to_one_within_genus(self):
        module = load_module()
        rows = module.genus_level_deltas(balanced_frame())
        self.assertEqual([row["genus"] for row in rows], ["G1", "G2"])
        for row in rows:
            self.assertEqual(row["n_carrier"], 2)
            self.assertEqual(row["n_control"], 2)
            self.assertAlmostEqual(row["delta"], 0.5, places=6)

    def test_genus_deltas_accept_a_frame_filtered_from_a_larger_table(self):
        module = load_module()
        table = balanced_frame()
        filtered = table[table["group"].eq("intracellular")].copy()
        rows = module.genus_level_deltas(filtered)
        self.assertEqual([row["genus"] for row in rows], ["G1", "G2"])

    def test_group_comparison_is_genus_level_and_secondary(self):
        module = load_module()
        table = pd.DataFrame(
            [
                genome("I1", "G1", CARRIER, 0.60, group="intracellular"),
                genome("I2", "G1", CONTROL, 0.50, group="control"),
                genome("I3", "G2", CARRIER, 0.70, group="intracellular"),
                genome("I4", "G2", CONTROL, 0.60, group="control"),
                genome("E1", "G3", CARRIER, 1.00, group="extracellular"),
                genome("E2", "G3", CONTROL, 0.30, group="control"),
                genome("E3", "G4", CARRIER, 0.90, group="extracellular"),
                genome("E4", "G4", CONTROL, 0.30, group="control"),
            ]
        )
        result = module.intracellular_vs_extracellular(
            table, np.random.default_rng(3), n_permutations=200
        )
        self.assertEqual(result["analysis_role"], "secondary")
        self.assertEqual(result["resampling_unit"], "genus")
        self.assertEqual(result["n_intracellular_genera"], 2)
        self.assertEqual(result["n_extracellular_genera"], 2)
        self.assertAlmostEqual(result["mean_intracellular_delta"], 0.1, places=6)
        self.assertAlmostEqual(result["mean_extracellular_delta"], 0.65, places=6)
        self.assertAlmostEqual(result["observed"], -0.55, places=6)
        self.assertGreaterEqual(result["p_value"], 0.0)
        self.assertLessEqual(result["p_value"], 1.0)

    def test_mixed_genomes_are_excluded_from_the_compartment_comparison(self):
        module = load_module()
        table = pd.DataFrame(
            [
                genome("I1", "G1", CARRIER, 0.9, group="intracellular"),
                genome("I2", "G1", CONTROL, 0.2, group="control"),
                genome("X1", "G2", CARRIER, 0.9, group="mixed"),
                genome("X2", "G2", CONTROL, 0.2, group="control"),
            ]
        )
        result = module.intracellular_vs_extracellular(
            table, np.random.default_rng(3), n_permutations=50
        )
        self.assertEqual(result["n_extracellular_genera"], 0)


class ReusedPredictionAuditTests(unittest.TestCase):
    """Requirement 4: only comparable reused predictions enter the primary result."""

    def test_comparable_records_are_accepted(self):
        module = load_module()
        compatible, incompatible = module.audit_reused_predictions(
            [reused_record("G1"), reused_record("G2")], reference_profile()
        )
        self.assertEqual([row["genome_id"] for row in compatible], ["G1", "G2"])
        self.assertEqual(incompatible, [])

    def test_incompatible_tool_version_names_the_mismatched_field(self):
        module = load_module()
        compatible, incompatible = module.audit_reused_predictions(
            [reused_record("G1"), reused_record("G2", tool_version="1.0.0")],
            reference_profile(),
        )
        self.assertEqual([row["genome_id"] for row in compatible], ["G1"])
        self.assertEqual(len(incompatible), 1)
        self.assertEqual(incompatible[0]["genome_id"], "G2")
        self.assertEqual(incompatible[0]["mismatch_fields"], ["tool_version"])
        self.assertEqual(incompatible[0]["comparable"], False)
        self.assertIn("tool_version", incompatible[0]["reason"])

    def test_incompatible_output_formula_names_the_mismatched_field(self):
        module = load_module()
        record = reused_record(
            "G3",
            output_formula="growth_rate_per_h=1/doubling_time_h",
            growth_rate_per_h=1.0 / 2.0,
        )
        compatible, incompatible = module.audit_reused_predictions(
            [record], reference_profile()
        )
        self.assertEqual(compatible, [])
        self.assertEqual(incompatible[0]["mismatch_fields"], ["output_formula"])
        self.assertIn("output_formula", incompatible[0]["reason"])

    def test_each_comparability_dimension_is_individually_enforced(self):
        module = load_module()
        cases = {
            "tool": {"tool": "other-tool"},
            "tool_version": {"tool_version": "2.5.0"},
            "model_settings": {"model_settings": "legacy_default"},
            "temperature_handling": {"temperature_handling": "no_temperature_recorded"},
            "cds_preparation": {"cds_preparation": "prodigal_v2.6.3"},
            "ribosomal_marker_procedure": {"ribosomal_marker_procedure": "rfam_5S_blast"},
            "output_formula": {
                "output_formula": "growth_rate_per_h=ln2/doubling_time_h_minutes",
            },
        }
        for field, override in cases.items():
            with self.subTest(field=field):
                compatible, incompatible = module.audit_reused_predictions(
                    [reused_record("G1", **override)], reference_profile()
                )
                self.assertEqual(compatible, [])
                self.assertEqual(incompatible[0]["mismatch_fields"], [field])

    def test_a_record_can_mismatch_on_several_fields_at_once(self):
        module = load_module()
        record = reused_record(
            "G4", tool_version="1.0.0", cds_preparation="prodigal_v2.6.3"
        )
        _, incompatible = module.audit_reused_predictions([record], reference_profile())
        self.assertEqual(
            incompatible[0]["mismatch_fields"], ["cds_preparation", "tool_version"]
        )

    def test_incompatible_records_are_retained_in_the_sensitivity_table(self):
        module = load_module()
        records = [reused_record("G1"), reused_record("G2", tool_version="1.0.0")]
        audit = module.audit_reused_predictions_table(records, reference_profile())
        self.assertEqual(audit["n_records"], 2)
        self.assertEqual(audit["n_compatible"], 1)
        self.assertEqual(audit["n_incompatible"], 1)
        self.assertEqual(audit["primary_genome_ids"], ["G1"])
        self.assertEqual(
            [row["genome_id"] for row in audit["sensitivity_table"]], ["G2"]
        )
        for row in audit["sensitivity_table"]:
            self.assertTrue(row["reason"])
            self.assertEqual(row["used_in_primary_analysis"], False)
            self.assertEqual(row["used_in_sensitivity_analysis"], True)
        self.assertEqual(audit["status"], "ok")

    def test_sensitivity_rows_carry_the_prediction_payload(self):
        module = load_module()
        records = [reused_record("G2", tool_version="1.0.0", doubling_time_h=3.0,
                                 growth_rate_per_h=LN2 / 3.0)]
        audit = module.audit_reused_predictions_table(records, reference_profile())
        row = audit["sensitivity_table"][0]
        self.assertAlmostEqual(row["doubling_time_h"], 3.0)
        self.assertEqual(row["source_table"], "reused_predictions")

    def test_a_record_missing_a_comparability_field_is_incompatible(self):
        module = load_module()
        record = reused_record("G5")
        del record["temperature_handling"]
        compatible, incompatible = module.audit_reused_predictions(
            [record], reference_profile()
        )
        self.assertEqual(compatible, [])
        self.assertEqual(incompatible[0]["mismatch_fields"], ["temperature_handling"])

    def test_duplicate_reused_genomes_are_reported_not_silently_counted_twice(self):
        module = load_module()
        audit = module.audit_reused_predictions_table(
            [reused_record("G1"), reused_record("G1")], reference_profile()
        )
        self.assertEqual(audit["n_records"], 2)
        self.assertEqual(audit["n_unique_genomes"], 1)
        self.assertEqual(audit["duplicate_genome_ids"], ["G1"])


class SelectionCoverageTests(unittest.TestCase):
    """Requirement 5: coverage is quantified by named, non-aggregate reasons."""

    def coverage_table(self):
        return pd.DataFrame(
            [
                {
                    "genome_id": "A1", "taxonomy": "Pseudomonadota;GenA",
                    "candidate_family": "DED_hfam_70", "group": "intracellular",
                    "quality": "high", "selection_status": "included",
                },
                {
                    "genome_id": "A2", "taxonomy": "Pseudomonadota;GenA",
                    "candidate_family": "DED_hfam_70", "group": "intracellular",
                    "quality": "high", "selection_status": "included",
                },
                {
                    "genome_id": "B1", "taxonomy": "Bacillota;GenB",
                    "candidate_family": "DED_hfam_2", "group": "extracellular",
                    "quality": "medium", "selection_status": "excluded_no_same_genus_control",
                },
                {
                    "genome_id": "C1", "taxonomy": "Bacillota;GenB",
                    "candidate_family": "DED_hfam_2", "group": "extracellular",
                    "quality": "medium", "selection_status": "excluded_missing_fasta",
                },
                {
                    "genome_id": "D1", "taxonomy": "Actinomycetota;GenC",
                    "candidate_family": "DED_hfam_52", "group": "intracellular",
                    "quality": "low", "selection_status": "excluded_failed_prediction",
                },
            ]
        )

    def test_summary_counts_every_selection_status(self):
        module = load_module()
        summary = module.summarize_selection_coverage(self.coverage_table())
        self.assertEqual(summary["total_genomes"], 5)
        self.assertEqual(summary["included_count"], 2)
        self.assertEqual(summary["excluded_count"], 3)
        nonzero = {
            name: value for name, value in summary["status_counts"].items() if value
        }
        self.assertEqual(
            nonzero,
            {
                "included": 2,
                "excluded_no_same_genus_control": 1,
                "excluded_missing_fasta": 1,
                "excluded_failed_prediction": 1,
            },
        )
        self.assertEqual(
            sum(summary["status_counts"].values()), summary["total_genomes"]
        )

    def test_summary_stratifies_by_taxonomy_family_group_and_quality(self):
        module = load_module()
        summary = module.summarize_selection_coverage(self.coverage_table())
        self.assertEqual(summary["by_taxonomy"]["Pseudomonadota;GenA"]["included"], 2)
        self.assertEqual(
            summary["by_taxonomy"]["Bacillota;GenB"]["excluded_missing_fasta"], 1
        )
        self.assertEqual(summary["by_candidate_family"]["DED_hfam_70"]["total"], 2)
        self.assertEqual(
            summary["by_candidate_family"]["DED_hfam_2"][
                "excluded_no_same_genus_control"
            ],
            1,
        )
        self.assertEqual(summary["by_group"]["intracellular"]["included"], 2)
        self.assertEqual(summary["by_group"]["extracellular"]["excluded_count"], 2)
        self.assertEqual(summary["by_quality"]["low"]["excluded_failed_prediction"], 1)
        self.assertEqual(
            summary["by_quality"]["high"]["included_fraction"], 1.0
        )

    def test_unknown_selection_status_is_rejected(self):
        module = load_module()
        table = self.coverage_table()
        table.loc[0, "selection_status"] = "excluded_because_we_say_so"
        with self.assertRaisesRegex(ValueError, "selection_status"):
            module.summarize_selection_coverage(table)

    def test_missing_column_is_rejected(self):
        module = load_module()
        table = self.coverage_table().drop(columns=["quality"])
        with self.assertRaises(ValueError):
            module.summarize_selection_coverage(table)

    def test_require_no_control_status_is_a_named_bucket(self):
        module = load_module()
        summary = module.summarize_selection_coverage(self.coverage_table())
        self.assertEqual(
            summary["no_same_genus_control_count"],
            summary["status_counts"]["excluded_no_same_genus_control"],
        )
        self.assertEqual(summary["missing_fasta_count"], 1)
        self.assertEqual(summary["failed_prediction_count"], 1)


class ManifestDifferenceTests(unittest.TestCase):
    """Requirement 6: the 4,507 -> 4,441 = 66 gap from accession-level sets."""

    def split_66(self):
        """66 synthetic eligible positives that do not reach the manifest.

        Buckets: 30 genomes with no same-genus control, 20 missing FASTA,
        14 failed prediction (too few ribosomal hits), 2 detected as carriers
        only outside the balanced same-genus design.
        """
        reasons = {}
        reasons["no_same_genus_control"] = {f"NOCTL_{i:03d}" for i in range(30)}
        reasons["missing_fasta"] = {f"NOFASTA_{i:03d}" for i in range(20)}
        reasons["failed_prediction_too_few_ribosomal_hits"] = {
            f"RIBOHIT_{i:03d}" for i in range(14)
        }
        reasons["not_reaching_the_balanced_same_genus_design"] = {
            f"UNBALANCED_{i:03d}" for i in range(2)
        }
        return reasons

    def test_reason_buckets_sum_to_exactly_sixty_six(self):
        module = load_module()
        reasons = self.split_66()
        eligible = set().union(*reasons.values()) | {f"KEPT_{i:03d}" for i in range(10)}
        manifest = eligible - set().union(*reasons.values())
        payload = module.explain_manifest_difference(
            eligible_positive_accessions=eligible,
            manifest_positive_accessions=manifest,
            reason_buckets=reasons,
        )
        self.assertEqual(payload["status"], "ok")
        self.assertEqual(payload["eligible_positive_count"], 76)
        self.assertEqual(payload["manifest_positive_count"], 10)
        self.assertEqual(payload["difference"], 66)
        self.assertEqual(
            sum(payload["reason_bucket_counts"].values()), payload["difference"]
        )
        self.assertEqual(
            payload["reason_bucket_counts"],
            {
                "failed_prediction_too_few_ribosomal_hits": 14,
                "missing_fasta": 20,
                "no_same_genus_control": 30,
                "not_reaching_the_balanced_same_genus_design": 2,
            },
        )
        self.assertEqual(payload["unexplained_count"], 0)
        self.assertEqual(payload["unexplained_accessions"], [])
        self.assertEqual(payload["overlap_count"], 0)
        self.assertEqual(payload["source"], "accession_level_sets")

    def test_counts_come_from_accession_sets_not_from_subtracting_aggregates(self):
        module = load_module()
        eligible = {f"E{i:03d}" for i in range(40)}
        # The manifest has the SAME positive count as the eligible set, but a
        # different membership: 28 accessions are missing and 28 unrelated ones
        # were counted instead.  Aggregate subtraction reports "0 difference";
        # accession-level sets expose the 28 unexplained and the 28 ineligible.
        manifest = (eligible - {f"E{i:03d}" for i in range(12, 40)}) | {
            f"X{i:03d}" for i in range(28)
        }
        self.assertEqual(len(manifest), len(eligible))
        payload = module.explain_manifest_difference(
            eligible_positive_accessions=eligible,
            manifest_positive_accessions=manifest,
            reason_buckets={"no_same_genus_control": {f"X{i:03d}" for i in range(28)}},
        )
        self.assertEqual(payload["eligible_positive_count"], 40)
        self.assertEqual(payload["manifest_positive_count"], 40)
        self.assertEqual(payload["difference"], 28)
        self.assertEqual(payload["status"], "mismatch")
        self.assertEqual(payload["manifest_not_eligible_count"], 28)
        self.assertEqual(
            payload["manifest_not_eligible_accessions"], [f"X{i:03d}" for i in range(28)]
        )
        self.assertEqual(payload["unexplained_count"], 28)
        self.assertEqual(
            payload["unexplained_accessions"], [f"E{i:03d}" for i in range(12, 40)]
        )
        self.assertEqual(payload["bucket_members_outside_eligible_count"], 28)
        self.assertEqual(payload["source"], "accession_level_sets")

    def test_bucket_members_outside_the_eligible_set_are_reported(self):
        module = load_module()
        payload = module.explain_manifest_difference(
            eligible_positive_accessions={"A", "B", "C"},
            manifest_positive_accessions={"A", "B", "C"},
            reason_buckets={"no_same_genus_control": {"Z"}},
        )
        self.assertEqual(payload["difference"], 0)
        self.assertEqual(payload["bucket_members_outside_eligible_count"], 1)
        self.assertEqual(payload["bucket_members_outside_eligible_accessions"], ["Z"])

    def test_membership_is_decided_by_accession_identity(self):
        module = load_module()
        eligible = {"A", "B", "C", "D", "E", "F"}
        manifest = {"A", "B", "C", "D", "E"}
        payload = module.explain_manifest_difference(
            eligible_positive_accessions=eligible,
            manifest_positive_accessions=manifest,
            reason_buckets={"missing_fasta": {"F"}},
        )
        self.assertEqual(payload["difference"], 1)
        self.assertEqual(payload["reason_bucket_counts"], {"missing_fasta": 1})
        self.assertEqual(payload["reason_bucket_accessions"]["missing_fasta"], ["F"])

    def test_an_unexplained_genome_is_reported(self):
        module = load_module()
        payload = module.explain_manifest_difference(
            eligible_positive_accessions={"A", "B", "C"},
            manifest_positive_accessions={"A"},
            reason_buckets={"missing_fasta": {"B"}},
        )
        self.assertEqual(payload["status"], "mismatch")
        self.assertEqual(payload["difference"], 2)
        self.assertEqual(payload["unexplained_count"], 1)
        self.assertEqual(payload["unexplained_accessions"], ["C"])

    def test_a_genome_explained_by_two_buckets_is_a_mismatch(self):
        module = load_module()
        payload = module.explain_manifest_difference(
            eligible_positive_accessions={"A", "B"},
            manifest_positive_accessions={"A"},
            reason_buckets={"missing_fasta": {"B"}, "failed_prediction": {"B"}},
        )
        self.assertEqual(payload["status"], "mismatch")
        self.assertEqual(payload["overlap_count"], 1)
        self.assertEqual(payload["overlap_accessions"], ["B"])

    def test_expected_counts_pin_the_plan_numbers(self):
        module = load_module()
        reasons = self.split_66()
        eligible = set().union(*reasons.values())
        payload = module.explain_manifest_difference(
            eligible_positive_accessions=eligible,
            manifest_positive_accessions=set(),
            reason_buckets=reasons,
            expected_difference=66,
        )
        self.assertEqual(payload["difference"], 66)
        self.assertEqual(payload["expected_difference"], 66)
        self.assertEqual(payload["status"], "ok")

    def test_expected_counts_that_do_not_match_are_a_mismatch(self):
        module = load_module()
        payload = module.explain_manifest_difference(
            eligible_positive_accessions={"A", "B"},
            manifest_positive_accessions={"A"},
            reason_buckets={"missing_fasta": {"B"}},
            expected_difference=66,
        )
        self.assertEqual(payload["status"], "mismatch")
        self.assertEqual(payload["expected_difference"], 66)
        self.assertEqual(payload["difference"], 1)

    def test_manifest_positives_outside_the_eligible_set_are_a_mismatch(self):
        module = load_module()
        payload = module.explain_manifest_difference(
            eligible_positive_accessions={"A"},
            manifest_positive_accessions={"A", "Z"},
            reason_buckets={},
        )
        self.assertEqual(payload["status"], "mismatch")
        self.assertEqual(payload["manifest_not_eligible_accessions"], ["Z"])

    def test_empty_eligible_set_is_rejected(self):
        module = load_module()
        with self.assertRaises(ValueError):
            module.explain_manifest_difference(
                eligible_positive_accessions=set(),
                manifest_positive_accessions=set(),
                reason_buckets={"missing_fasta": {"A"}},
            )


class ModuleBoundaryTests(unittest.TestCase):
    def test_module_is_pure_and_lives_outside_the_untouched_deploy(self):
        module = load_module()
        self.assertFalse(hasattr(module, "RUN_GRODON"))
        self.assertIn("deploy/20260928_phaded_grodon_reanalysis_v2_01", module.__doc__)
        self.assertEqual(Path(module.__file__).resolve(), SCRIPT.resolve())

    def test_legacy_deploy_scripts_still_exist_and_stay_untouched(self):
        self.assertTrue((LEGACY_DEPLOY / "scripts" / "build_grodon_manifest.py").is_file())
        self.assertTrue((LEGACY_DEPLOY / "scripts" / "grodon_group_stats.py").is_file())

    def test_the_new_dated_deploy_directory_is_not_created_by_this_task(self):
        self.assertFalse(NEW_DEPLOY.exists())

    def test_importing_the_module_does_not_import_grodon_or_r_interfaces(self):
        module = load_module()
        for name in ("subprocess", "shutil", "pyrodigal"):
            self.assertFalse(hasattr(module, name), name)


if __name__ == "__main__":
    unittest.main()
