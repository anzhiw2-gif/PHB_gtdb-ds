"""Tests for the gRodon 66-record manifest-delta reconciliation (plan Task F2).

Written **before** the implementation (failing-test-first).  ``reconcile_phaded_
grodon_manifest_delta.py`` closes, or refuses to fake, the gap between the
documented 4,507 "eligible positives" and the 4,441 manifest positives of the
frozen run ``20260920_phaded_grodon_growth_01``.

What these tests pin down
-------------------------

1. **Buckets are accession sets.**  Every bucket count is ``len(set)`` and the
   counts must sum to exactly the difference.  Two buckets may never claim the
   same accession.
2. **Aggregate subtraction is not evidence.**  The negative control builds a
   manifest whose *totals* match the eligible set exactly but whose *members*
   differ; aggregate subtraction reports a difference of 0 and hides everything,
   while set logic must expose both the disqualified and the unexplained
   accessions.  This is the whole point of the module.
3. **Unexplained accessions are never silently dropped.**  They are listed by
   accession in the summary and force a non-``ok`` status.
4. **An unresolvable difference is reported as pending, never invented.**  A
   bucket that cannot be reconstructed from the frozen evidence contributes
   ``0`` accessions and names the missing input; the identity check then
   compares the *verified* difference, not the aggregate hypothesis.
5. **Missing or malformed inputs raise.**  A missing file, a missing column or a
   mismatched declared count is an error, not a zero.

Every fixture is synthetic.  Nothing here reads or writes ``runs/``,
``results/`` or ``deploy/``, executes gRodon/R, or opens the network.
"""

import contextlib
import importlib.util
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "reconcile_phaded_grodon_manifest_delta.py"
)
REUSE_GRODON = Path(__file__).resolve().parents[1] / "scripts" / "grodon_reanalysis_v2.py"
REUSE_FLOW = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "reconcile_phaded_candidate_flow.py"
)


def load_module(name, path):
    """Load a module by file path, reusing an already-loaded copy if present.

    Two modules load each other's helpers by path, so a plain second load would
    register a second module object under the same name and can break
    ``@dataclass``-decorated types.  Caching by name keeps the identity stable.
    """
    cached = sys.modules.get(name)
    if cached is not None and getattr(cached, "__file__", None) == str(path):
        return cached
    scripts_dir = str(path.parent)
    if scripts_dir not in sys.path:
        sys.path.insert(0, scripts_dir)
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


module = load_module("reconcile_phaded_grodon_manifest_delta", SCRIPT)


# --------------------------------------------------------------------- helpers


def gid(i):
    """9-digit GTDB accession, e.g. GCA_000000001.1."""
    return "GCA_%09d.1" % i


def write_text(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8", newline="\n")


def write_table(path, header, rows):
    lines = ["\t".join(header)]
    for row in rows:
        lines.append("\t".join("" if c is None else str(c) for c in row))
    write_text(path, "\n".join(lines) + "\n")


MANIFEST_HEADER = [
    "genome_id",
    "phaZ_status",
    "phaZ_validated_count",
    "major_subtype",
    "group",
    "domain",
    "genus",
    "genome_path",
]
FAILED_HEADER = ["genome_id", "phaZ_status", "status", "error", "genome_path"]
DEGRADER_HEADER = [
    "genome_id",
    "major_superfamily",
    "n_superfamilies",
    "n_candidates",
    "group",
]


def manifest_rows(positive_ids, negative_ids):
    rows = []
    for i in positive_ids:
        rows.append(
            [i, "phaZ_positive", 1, "intracellular nPHASCL without lipase box",
             "intracellular", "Bacteria", "TestGenus", "/g/%s.fna.gz" % i]
        )
    for i in negative_ids:
        rows.append(
            [i, "phaZ_negative", 0, "none", "control", "Bacteria", "TestGenus",
             "/g/%s.fna.gz" % i]
        )
    return rows


def fixture(tmp, eligible, manifest_pos, manifest_neg, failed_pos=(),
            failed_neg=(), degrader=(), expected=None):
    """Materialise a synthetic frozen-run shaped directory tree."""
    root = Path(tmp)
    write_table(
        root / "results" / "grodon_growth_manifest_newdeg40k.tsv",
        MANIFEST_HEADER,
        manifest_rows(manifest_pos, manifest_neg),
    )
    write_table(
        root / "results" / "grodon_failed_genomes_newdeg40k.tsv",
        FAILED_HEADER,
        [
            [i, "phaZ_positive", "failed", "too_few_ribosomal_hits:7", "/g/%s.fna.gz" % i]
            for i in failed_pos
        ]
        + [
            [i, "phaZ_negative", "failed", "too_few_ribosomal_hits:9", "/g/%s.fna.gz" % i]
            for i in failed_neg
        ],
    )
    write_table(
        root / "inputs" / "degrader_genomes.tsv",
        DEGRADER_HEADER,
        [[i, "intracellular nPHASCL without lipase box", 1, 1, "intracellular"] for i in degrader],
    )
    write_text(
        root / "inputs" / "pool_genomes_exclusion.txt",
        "".join("%s\n" % i for i in (expected or {}).get("exclusion", [])),
    )
    stats = {
        "degrader_input_genomes": len(degrader),
        "exclusion_genomes": len((expected or {}).get("exclusion", [])),
        "manifest_rows": len(manifest_pos) + len(manifest_neg),
        "manifest_positive": len(manifest_pos),
        "manifest_negative": len(manifest_neg),
        "degrader_genomes_excluded_no_control": (expected or {}).get("excluded_no_control"),
        "degrader_genomes_missing_fasta": (expected or {}).get("missing_fasta", 0),
    }
    if stats["degrader_genomes_excluded_no_control"] is None:
        stats["degrader_genomes_excluded_no_control"] = (
            len(degrader) - len(manifest_pos)
        )
    write_text(root / "results" / "manifest_stats.json", json.dumps(stats, indent=2) + "\n")
    return root


def tree(tmp):
    """Read the synthetic run tree into the module's own input bundle."""
    return module.load_run_inputs(Path(tmp))


# ------------------------------------------------------------------- the tests


class ReconcileManifestDeltaTest(unittest.TestCase):
    """The pure set-algebra function: buckets must partition the difference."""

    def test_bucket_counts_sum_exactly_to_the_difference(self):
        eligible = {gid(i) for i in range(1, 11)}
        manifest = {gid(i) for i in range(1, 6)}
        buckets = {
            "balance_dropped": {gid(6), gid(7)},
            "prediction_failed": {gid(8)},
            "dedupe": {gid(9), gid(10)},
        }
        rows, payload = module.reconcile_manifest_delta(
            eligible, manifest, buckets
        )
        self.assertEqual(payload["eligible_count"], 10)
        self.assertEqual(payload["manifest_positive_count"], 5)
        self.assertEqual(payload["difference"], 5)
        self.assertEqual(payload["sum_of_bucket_counts"], 5)
        self.assertEqual(
            payload["bucket_counts"],
            {"balance_dropped": 2, "prediction_failed": 1, "dedupe": 2},
        )
        self.assertTrue(payload["identity_check"]["passes"])
        self.assertEqual(payload["identity_check"]["difference"], 5)
        self.assertEqual(payload["unexplained_count"], 0)
        self.assertEqual(payload["overlap_count"], 0)
        self.assertEqual(payload["status"], "ok")
        self.assertEqual(len(rows), 5)
        self.assertEqual(
            sorted({row["accession"] for row in rows}),
            sorted(eligible - manifest),
        )
        for row in rows:
            self.assertIn(row["bucket"], buckets)
            self.assertEqual(row["direction"], "eligible_not_in_manifest")
            self.assertTrue(row["evidence"])

    def test_no_accession_in_two_buckets(self):
        eligible = {gid(i) for i in range(1, 6)}
        manifest = {gid(1)}
        buckets = {
            "a": {gid(2), gid(3)},
            "b": {gid(3), gid(4)},
        }
        with self.assertRaises(module.BucketOverlapError) as caught:
            module.reconcile_manifest_delta(eligible, manifest, buckets)
        self.assertIn(gid(3), str(caught.exception))

    def test_negative_control_counts_match_but_membership_differs(self):
        """Aggregate subtraction sees nothing; set logic must see everything.

        ``eligible`` and ``manifest`` have the same size, so any aggregate
        subtraction reports a difference of 0.  The members differ: 3 eligible
        genomes are absent from the manifest, and 3 manifest genomes were never
        eligible.  Set algebra must expose both sets, and the result must not be
        a silently clean ``ok``.
        """
        eligible = {gid(i) for i in range(1, 8)}
        manifest = {gid(i) for i in range(5, 12)}
        self.assertEqual(len(eligible), len(manifest))

        rows, payload = module.reconcile_manifest_delta(
            eligible, manifest, {"no_bucket_claims_these": set()}
        )
        # Aggregate subtraction is blind here; set algebra is not.
        self.assertEqual(payload["aggregate_subtraction_difference"], 0)
        self.assertEqual(payload["difference"], 4)
        self.assertEqual(payload["manifest_positive_count"], 7)
        self.assertEqual(payload["eligible_count"], 7)
        self.assertEqual(payload["unexplained_count"], 4)
        self.assertEqual(payload["manifest_not_eligible_count"], 4)
        self.assertEqual(
            payload["unexplained_accessions"], sorted(gid(i) for i in range(1, 5))
        )
        self.assertEqual(
            payload["manifest_not_eligible_accessions"],
            sorted(gid(i) for i in range(8, 12)),
        )
        self.assertNotEqual(payload["status"], "ok")
        # The identity is over named accessions only: no bucket names any of the
        # three, so it must fail rather than pass on the coincidence that the
        # totals match.
        self.assertFalse(payload["identity_check"]["passes"])
        self.assertEqual(payload["identity_check"]["sum_of_bucket_counts"], 0)
        self.assertEqual(
            sorted(row["accession"] for row in rows),
            sorted(eligible - manifest) + sorted(manifest - eligible),
        )

    def test_unexplained_accessions_are_listed_never_dropped(self):
        eligible = {gid(i) for i in range(1, 6)}
        manifest = {gid(1)}
        buckets = {"partial": {gid(2)}}
        _, payload = module.reconcile_manifest_delta(eligible, manifest, buckets)
        self.assertEqual(payload["unexplained_count"], 3)
        self.assertEqual(
            payload["unexplained_accessions"], sorted([gid(3), gid(4), gid(5)])
        )
        self.assertEqual(payload["sum_of_bucket_counts"], 1)
        self.assertEqual(payload["difference"], 4)
        self.assertFalse(payload["identity_check"]["passes"])
        self.assertEqual(payload["status"], "mismatch")

    def test_pending_bucket_contributes_zero_and_names_the_missing_input(self):
        """An unresolvable bucket is reported beside the identity, never inside it."""
        eligible = {gid(i) for i in range(1, 6)}
        manifest = {gid(1)}
        _, payload = module.reconcile_manifest_delta(
            eligible,
            manifest,
            {"exactly_one_reason": {gid(2)}},
            pending_buckets={
                "shuffled_out_of_the_1to1_balance": {
                    "count": 63,
                    "missing_input": "GTDB bac120_taxonomy_r232.tsv + genome FASTA presence",
                    "reason": "selection order comes from random.shuffle(seed=42)",
                }
            },
        )
        self.assertEqual(payload["bucket_counts"]["exactly_one_reason"], 1)
        self.assertEqual(
            payload["pending_buckets"]["shuffled_out_of_the_1to1_balance"]["count"],
            63,
        )
        self.assertEqual(
            payload["pending_buckets"]["shuffled_out_of_the_1to1_balance"]["accessions"],
            0,
        )
        self.assertEqual(payload["pending_count"], 63)
        self.assertEqual(payload["sum_of_bucket_counts"], 1)
        # The identity is over named accessions only and it genuinely fails: the
        # pending count must not be able to make it pass.
        self.assertEqual(payload["difference"], 4)
        self.assertEqual(payload["named_difference_count"], 1)
        self.assertFalse(payload["identity_check"]["passes"])
        self.assertEqual(payload["identity_check"]["pending_excluded_count"], 63)
        self.assertEqual(payload["status"], "mismatch")

    def test_pending_bucket_is_absent_from_a_fully_named_difference(self):
        eligible = {gid(i) for i in range(1, 6)}
        manifest = {gid(1)}
        _, payload = module.reconcile_manifest_delta(
            eligible,
            manifest,
            {"named": {gid(2), gid(3), gid(4), gid(5)}},
            pending_buckets={
                "may_or_may_not_be_eligible": {
                    "count": 66,
                    "missing_input": "GTDB taxonomy + FASTA presence map",
                    "reason": "the scalar sums two reasons and names no accession",
                }
            },
        )
        self.assertTrue(payload["identity_check"]["passes"])
        self.assertEqual(payload["identity_check"]["sum_of_bucket_counts"], 4)
        self.assertEqual(payload["identity_check"]["pending_excluded_count"], 66)
        self.assertEqual(payload["status"], "ok")
        self.assertEqual(payload["unexplained_count"], 0)

    def test_pending_bucket_may_carry_its_own_accessions(self):
        eligible = {gid(i) for i in range(1, 6)}
        manifest = {gid(1)}
        _, payload = module.reconcile_manifest_delta(
            eligible,
            manifest,
            {"named": {gid(2), gid(3)}},
            pending_buckets={
                "partially_resolved": {
                    "count": 2,
                    "accessions": [gid(4), gid(5)],
                    "missing_input": "",
                    "reason": "",
                }
            },
        )
        spec = payload["pending_buckets"]["partially_resolved"]
        self.assertEqual(spec["count"], 2)
        self.assertEqual(spec["accessions"], 2)
        self.assertTrue(spec["resolvable"])

    def test_pending_bucket_without_a_missing_input_is_rejected(self):
        with self.assertRaises(ValueError) as caught:
            module.reconcile_manifest_delta(
                {gid(1), gid(2)},
                {gid(1)},
                {},
                pending_buckets={"guessed": {"count": 4, "missing_input": "", "reason": ""}},
            )
        self.assertIn("missing input", str(caught.exception))

    def test_empty_eligible_set_is_an_error(self):
        with self.assertRaises(ValueError):
            module.reconcile_manifest_delta(set(), {gid(1)}, {})

    def test_duplicate_accessions_within_one_bucket_collapse(self):
        rows, payload = module.reconcile_manifest_delta(
            {gid(1), gid(2)}, {gid(1)}, {"b": [gid(2), gid(2), gid(2)]}
        )
        self.assertEqual(payload["bucket_counts"]["b"], 1)
        self.assertEqual(len(rows), 1)


class LoadRunInputsTest(unittest.TestCase):
    """Reading the frozen-run shaped tree by set operations, never by subtraction."""

    def test_reads_counts_and_accession_sets(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture(
                tmp,
                eligible=None,
                manifest_pos=[gid(i) for i in range(1, 6)],
                manifest_neg=[gid(i) for i in range(101, 106)],
                failed_pos=[gid(6), gid(7)],
                failed_neg=[gid(8)],
                degrader=[gid(i) for i in range(1, 11)],
                expected={"exclusion": [gid(9)], "excluded_no_control": 3},
            )
            bundle = tree(tmp)
        self.assertEqual(len(bundle.degrader_input), 10)
        self.assertEqual(len(bundle.manifest_positive), 5)
        self.assertEqual(len(bundle.manifest_negative), 5)
        self.assertEqual(len(bundle.failed_positive), 2)
        self.assertEqual(len(bundle.failed_negative), 1)
        self.assertEqual(bundle.exclusion, {gid(9)})
        self.assertEqual(bundle.declared["degrader_input_genomes"], 10)
        self.assertEqual(bundle.declared["manifest_positive"], 5)

    def test_manifest_positive_and_negative_are_disjoint(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture(
                tmp,
                eligible=None,
                manifest_pos=[gid(i) for i in range(1, 4)],
                manifest_neg=[gid(i) for i in range(101, 104)],
                degrader=[gid(i) for i in range(1, 7)],
            )
            bundle = tree(tmp)
        self.assertFalse(bundle.manifest_positive & bundle.manifest_negative)

    def test_manifest_with_overlapping_arms_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture(
                tmp,
                eligible=None,
                manifest_pos=[gid(1), gid(2)],
                manifest_neg=[gid(3)],
                degrader=[gid(1), gid(2)],
            )
            write_table(
                Path(tmp) / "results" / "grodon_growth_manifest_newdeg40k.tsv",
                MANIFEST_HEADER,
                manifest_rows([gid(1), gid(2)], [gid(1), gid(3)]),
            )
            with self.assertRaises(ValueError) as caught:
                tree(tmp)
        self.assertIn("overlap", str(caught.exception).lower())

    def test_failed_table_without_phaz_status_raises(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture(
                tmp,
                eligible=None,
                manifest_pos=[gid(1)],
                manifest_neg=[gid(101)],
                degrader=[gid(1)],
            )
            write_table(
                Path(tmp) / "results" / "grodon_failed_genomes_newdeg40k.tsv",
                ["genome_id", "status", "error"],
                [[gid(2), "failed", "too_few_ribosomal_hits:3"]],
            )
            with self.assertRaises(ValueError) as caught:
                tree(tmp)
        self.assertIn("phaZ_status", str(caught.exception))

    def test_manifest_without_phaz_status_raises(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture(
                tmp,
                eligible=None,
                manifest_pos=[gid(1)],
                manifest_neg=[gid(101)],
                degrader=[gid(1)],
            )
            write_table(
                Path(tmp) / "results" / "grodon_growth_manifest_newdeg40k.tsv",
                ["genome_id", "genus"],
                [[gid(1), "TestGenus"]],
            )
            with self.assertRaises(ValueError) as caught:
                tree(tmp)
        self.assertIn("phaZ_status", str(caught.exception))

    def test_declared_count_mismatch_raises(self):
        """A manifest_stats.json that contradicts the table is an error."""
        with tempfile.TemporaryDirectory() as tmp:
            fixture(
                tmp,
                eligible=None,
                manifest_pos=[gid(i) for i in range(1, 4)],
                manifest_neg=[gid(i) for i in range(101, 104)],
                degrader=[gid(i) for i in range(1, 7)],
            )
            stats_path = Path(tmp) / "results" / "manifest_stats.json"
            stats = json.loads(stats_path.read_text(encoding="utf-8"))
            stats["manifest_positive"] = 99
            write_text(stats_path, json.dumps(stats, indent=2) + "\n")
            with self.assertRaises(ValueError) as caught:
                tree(tmp)
        self.assertIn("manifest_positive", str(caught.exception))

    def test_missing_input_raises(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "results").mkdir(parents=True)
            with self.assertRaises(FileNotFoundError) as caught:
                module.load_run_inputs(root)
        self.assertIn("grodon_growth_manifest", str(caught.exception))

    def test_missing_manifest_stats_raises(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture(
                tmp,
                eligible=None,
                manifest_pos=[gid(1)],
                manifest_neg=[gid(101)],
                degrader=[gid(1)],
            )
            (Path(tmp) / "results" / "manifest_stats.json").unlink()
            with self.assertRaises(FileNotFoundError):
                module.load_run_inputs(Path(tmp))


class NoControlDefinitionTest(unittest.TestCase):
    """The exclusion-set hypothesis and its bound are computed, never assumed."""

    def test_exclusion_set_is_not_the_no_control_set(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture(
                tmp,
                eligible=None,
                manifest_pos=[gid(i) for i in range(1, 6)],
                manifest_neg=[gid(i) for i in range(101, 106)],
                degrader=[gid(i) for i in range(1, 11)],
                expected={"exclusion": [gid(i) for i in range(9, 11)], "excluded_no_control": 5},
            )
            bundle = tree(tmp)
        evidence = module.hypothesis_evidence(bundle)
        # degrader - exclusion is 8, but the run excluded 5 -> the exclusion file
        # is not the "no control" ledger and the residual hypothesis differs.
        self.assertEqual(evidence["degrader_minus_exclusion_count"], 8)
        self.assertEqual(evidence["declared_excluded_no_control"], 5)
        self.assertEqual(evidence["declared_residual_eligible"], 5)
        self.assertFalse(evidence["exclusion_file_is_the_no_control_set"])
        self.assertEqual(evidence["material_absent_from_manifest_count"], 5)


class AccountingClosureTest(unittest.TestCase):
    """``run_analysis`` closes the documented ledger without inventing members.

    The fixture reproduces the frozen run's *shape*: a declared excluded scalar
    that is smaller than the material remainder, a failure table split by
    ``phaZ_status``, and measured buckets covering what is actually in the
    manifest.  The identities under test are the ones the real run must satisfy.
    """

    def _analyse(self, tmp, declared_excluded):
        fixture(
            tmp,
            eligible=None,
            manifest_pos=[gid(i) for i in range(1, 11)],
            manifest_neg=[gid(i) for i in range(101, 111)],
            failed_pos=[gid(1), gid(2), gid(3)],
            failed_neg=[gid(101)],
            degrader=[gid(i) for i in range(1, 21)],
            expected={"exclusion": [gid(11)], "excluded_no_control": declared_excluded},
        )
        return module.run_analysis(Path(tmp))

    def test_declared_eligible_ledger_closes_with_one_named_residue(self):
        """eligible = manifest + failed + residue, with the residue left pending."""
        with tempfile.TemporaryDirectory() as tmp:
            analysis = self._analyse(tmp, declared_excluded=4)
        payload = analysis["payload"]
        ledger = payload["declared_eligible_ledger"]
        # degrader 20, declared excluded 4 -> eligible 16; manifest 10, failed 3
        # -> residue 3, which closes 16 = 10 + 3 + 3.
        self.assertEqual(analysis["evidence"]["declared_residual_eligible"], 16)
        self.assertEqual(ledger["components"]["manifest_positive_measured"], 10)
        self.assertEqual(ledger["components"]["failed_positive_measured"], 3)
        self.assertEqual(ledger["components"]["residue_count_without_accessions"], 3)
        self.assertEqual(ledger["sum_of_components"], 16)
        self.assertTrue(ledger["closes"])
        # No accession may be emitted for the residue.
        residue_spec = payload["pending_buckets"][
            "residue_unresolvable_from_the_frozen_evidence"
        ]
        self.assertEqual(residue_spec["count"], 3)
        self.assertEqual(residue_spec["accessions"], 0)
        self.assertFalse(residue_spec["resolvable"])
        self.assertTrue(residue_spec["missing_input"])

    def test_material_ledger_residue_is_reported_separately(self):
        """The input-side reading of the same scalar is reported, not merged."""
        with tempfile.TemporaryDirectory() as tmp:
            analysis = self._analyse(tmp, declared_excluded=4)
        material = analysis["payload"]["material_ledger_residue"]
        # absent_from_manifest = 20 - 10 = 10; declared excluded covered by the
        # manifest = 3 failed positives -> declared eligible and absent = 1.
        self.assertEqual(material["material_absent_from_manifest_count"], 10)
        self.assertEqual(material["declared_excluded_covered_by_the_manifest"], 3)
        self.assertEqual(material["declared_eligible_and_absent_count"], 1)
        self.assertEqual(material["residue_count"], 9)
        self.assertIn("note", material)

    def test_residue_is_never_folded_into_a_bucket(self):
        with tempfile.TemporaryDirectory() as tmp:
            analysis = self._analyse(tmp, declared_excluded=4)
        payload = analysis["payload"]
        # The residue contributes zero: the reason-bucket partition stays empty
        # because no bucket can name these accessions.
        self.assertEqual(payload["sum_of_bucket_counts"], 0)
        self.assertEqual(payload["pending_count"], 3 + 3)
        self.assertEqual(payload["pending_accessions_total"], 3)

    def test_accounted_total_covers_the_whole_input(self):
        with tempfile.TemporaryDirectory() as tmp:
            analysis = self._analyse(tmp, declared_excluded=1)
        payload = analysis["payload"]
        bundle = analysis["bundle"]
        self.assertEqual(
            payload["measured_total"] + payload["difference"], len(bundle.degrader_input)
        )
        self.assertEqual(payload["measured_total"], len(bundle.manifest_positive))

    def test_failed_positive_bucket_carries_its_accessions(self):
        with tempfile.TemporaryDirectory() as tmp:
            analysis = self._analyse(tmp, declared_excluded=1)
        spec = analysis["payload"]["pending_buckets"]["prediction_failed_positive_arm"]
        self.assertEqual(spec["count"], 3)
        self.assertEqual(spec["accessions"], 3)
        self.assertTrue(spec["resolvable"])

    def test_identities_hold_and_the_naive_sum_is_reported_as_false(self):
        with tempfile.TemporaryDirectory() as tmp:
            analysis = self._analyse(tmp, declared_excluded=4)
        payload = analysis["payload"]
        identities = payload["identities"]
        for name in (
            "input_split",
            "positive_arm_split",
            "control_arm_split",
            "documented_eligible_count",
        ):
            self.assertTrue(identities[name]["holds"], name)
            self.assertEqual(identities[name]["lhs"], identities[name]["rhs"], name)
        # The sum recorded in the handoff document must be shown as NOT holding.
        naive = identities["naive_sum_that_does_not_hold"]
        self.assertFalse(naive["holds"])
        self.assertNotEqual(naive["lhs"], naive["rhs"])
        self.assertTrue(naive["why"])

    def test_failed_positive_is_a_subset_of_the_manifest_positive_arm(self):
        with tempfile.TemporaryDirectory() as tmp:
            analysis = self._analyse(tmp, declared_excluded=4)
        bundle = analysis["bundle"]
        self.assertTrue(bundle.failed_positive <= bundle.manifest_positive)
        self.assertTrue(bundle.failed_negative <= bundle.manifest_negative)
        # ...and therefore NOT part of the material remainder.
        self.assertFalse(bundle.failed_positive & bundle.absent_from_manifest)
        self.assertFalse(bundle.failed_negative & bundle.absent_from_manifest)

        with self.assertRaises(module.BucketOverlapError):
            module.reconcile_manifest_delta(
                {gid(1), gid(2), gid(3)},
                {gid(1)},
                {"reason": {gid(2)}},
                measured_buckets={"measured": {gid(2)}},
            )

    def test_measured_bucket_overlapping_a_reason_bucket_is_rejected(self):
        with self.assertRaises(module.BucketOverlapError):
            module.reconcile_manifest_delta(
                {gid(1), gid(2), gid(3)},
                {gid(1)},
                {"reason": {gid(2)}},
                measured_buckets={"measured": {gid(2)}},
            )

    def test_measured_buckets_may_not_overlap_each_other(self):
        with self.assertRaises(module.BucketOverlapError):
            module.reconcile_manifest_delta(
                {gid(1), gid(2)},
                {gid(1)},
                {},
                measured_buckets={"a": {gid(1)}, "b": {gid(1)}},
            )

    def test_measured_bucket_rows_are_marked_in_manifest(self):
        rows, _ = module.reconcile_manifest_delta(
            {gid(1), gid(2)},
            {gid(1)},
            {},
            measured_buckets={"measured": {gid(1)}},
        )
        measured_rows = [row for row in rows if row["bucket"] == "measured"]
        self.assertEqual(len(measured_rows), 1)
        self.assertEqual(measured_rows[0]["direction"], "in_manifest")
        self.assertTrue(measured_rows[0]["evidence"])


class BucketTableCoverageTest(unittest.TestCase):
    """The bucket table accounts for the whole input, exactly once per accession."""

    def test_rows_are_unique_and_cover_the_symmetric_difference(self):
        eligible = {gid(i) for i in range(1, 8)}
        manifest = {gid(i) for i in range(5, 12)}
        rows, _ = module.reconcile_manifest_delta(eligible, manifest, {})
        accessions = [row["accession"] for row in rows]
        self.assertEqual(len(accessions), len(set(accessions)))
        # A plain call emits the two difference directions; the accessions the
        # two sets share carry no disposition and get no row.
        self.assertEqual(set(accessions), eligible ^ manifest)

    def test_measured_buckets_extend_coverage_to_the_whole_eligible_set(self):
        eligible = {gid(i) for i in range(1, 8)}
        manifest = {gid(i) for i in range(5, 12)}
        rows, payload = module.reconcile_manifest_delta(
            eligible, manifest, {}, measured_buckets={"in_manifest": manifest}
        )
        accessions = [row["accession"] for row in rows]
        self.assertEqual(len(accessions), len(set(accessions)))
        self.assertEqual(set(accessions), eligible | manifest)
        self.assertEqual(payload["measured_total"], len(manifest))
        # Every accession in the union gets exactly one row, and the ones in both
        # directions record the second bucket on the row instead of a second row.
        self.assertEqual(len(rows), len(eligible | manifest))
        by_accession = {row["accession"]: row for row in rows}
        # Eligible but not in the manifest: explained by no bucket.
        for accession in sorted(eligible - manifest):
            self.assertEqual(by_accession[accession]["bucket"], "unexplained")
            self.assertEqual(by_accession[accession]["direction"], "eligible_not_in_manifest")
        # In both: the measured bucket names the row; no accessor collision occurs
        # because a shared accession is in neither difference direction.
        for accession in sorted(eligible & manifest):
            row = by_accession[accession]
            self.assertEqual(row["bucket"], "in_manifest")
            self.assertEqual(row["direction"], "in_manifest")
            self.assertNotIn("also_in_buckets", row)
        # Manifest only: named by the reverse direction, flagged as also measured.
        for accession in sorted(manifest - eligible):
            row = by_accession[accession]
            self.assertEqual(row["bucket"], "manifest_not_eligible")
            self.assertEqual(row["also_in_buckets"], ["in_manifest"])

    def test_no_row_is_emitted_for_a_pending_bucket(self):
        rows, _ = module.reconcile_manifest_delta(
            {gid(1), gid(2)},
            {gid(1)},
            {},
            pending_buckets={
                "unknown": {"count": 5, "missing_input": "x.tsv", "reason": "not on disk"}
            },
        )
        self.assertEqual([row["accession"] for row in rows], [gid(2)])
        self.assertEqual(
            [row["bucket"] for row in rows], ["unexplained"]
        )


class RunDirectoryTest(unittest.TestCase):
    """Writing the run refuses to overwrite evidence and records its contract."""

    def _payload(self):
        _, payload = module.reconcile_manifest_delta(
            {gid(i) for i in range(1, 11)},
            {gid(i) for i in range(1, 6)},
            {"b1": {gid(6), gid(7)}, "b2": {gid(8)}, "b3": {gid(9), gid(10)}},
        )
        return payload

    def test_non_empty_output_dir_is_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "run"
            out.mkdir()
            (out / "leftover.txt").write_text("x", encoding="utf-8")
            with self.assertRaises(FileExistsError) as caught:
                module.write_run(
                    out,
                    self._payload(),
                    [{"accession": gid(6), "bucket": "b1", "evidence": "e"}],
                    input_rows=[],
                    provenance={},
                )
            self.assertIn("not empty", str(caught.exception))

    def test_creates_expected_layout_and_round_trips_json(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "20260928_run"
            module.write_run(
                out,
                self._payload(),
                [
                    {"accession": gid(6), "bucket": "b1", "evidence": "e1"},
                    {"accession": gid(7), "bucket": "b1", "evidence": "e1"},
                    {"accession": gid(8), "bucket": "b2", "evidence": "e2"},
                    {"accession": gid(9), "bucket": "b3", "evidence": "e3"},
                    {"accession": gid(10), "bucket": "b3", "evidence": "e3"},
                ],
                input_rows=[
                    {"path": "a.tsv", "size": 1, "sha256": "0" * 64, "role": "r"}
                ],
                provenance={"run_id": "20260928_run"},
            )
            for sub in ("logs", "inputs", "results"):
                self.assertTrue((out / sub).is_dir(), sub)
            payload = json.loads(
                (out / "results" / "grodon_66_reconciliation.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertTrue(payload["identity_check"]["passes"])
            self.assertEqual(payload["sum_of_bucket_counts"], payload["difference"])
            buckets = (out / "results" / "grodon_66_buckets.tsv").read_text(
                encoding="utf-8"
            )
            self.assertIn("accession\tbucket\tevidence", buckets)
            self.assertEqual(len(buckets.strip().splitlines()), 6)
            contract = json.loads(
                (out / "input_contract.json").read_text(encoding="utf-8")
            )
            self.assertEqual(contract["schema_version"], "1.0")
            self.assertEqual(contract["run_id"], "20260928_run")
            self.assertEqual(len(contract["inputs"]), 1)

    def test_empty_bucket_table_still_writes_a_header(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "run2"
            module.write_run(
                out, self._payload(), [], input_rows=[], provenance={}
            )
            text = (out / "results" / "grodon_66_buckets.tsv").read_text(
                encoding="utf-8"
            )
            self.assertEqual(text.strip().splitlines(), ["accession\tbucket\tevidence"])


class MainCliTest(unittest.TestCase):
    """``main`` is a checked entry point: int return, JSON summary, no writes in dry-run."""

    def _fixture(self, tmp):
        fixture(
            tmp,
            eligible=None,
            manifest_pos=[gid(i) for i in range(1, 11)],
            manifest_neg=[gid(i) for i in range(101, 111)],
            failed_pos=[gid(1), gid(2), gid(3)],
            failed_neg=[gid(101)],
            degrader=[gid(i) for i in range(1, 21)],
            expected={"exclusion": [gid(11)], "excluded_no_control": 4},
        )
        return Path(tmp)

    def _run_main(self, argv):
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            code = module.main(argv)
        return code, buffer.getvalue()

    def test_dry_run_returns_zero_and_emits_the_documented_closures(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self._fixture(tmp)
            code, out = self._run_main(["--run-root", str(root), "--dry-run"])
        self.assertEqual(code, 0)
        payload = json.loads(out)
        reading = payload["difference_reading"]
        self.assertEqual(reading["material_difference"]["count"], 10)
        self.assertEqual(reading["hypothesised_difference"]["count"], 3)
        self.assertTrue(reading["bucket_sum_closes_the_documented_figures"]["closes"])
        self.assertEqual(
            reading["bucket_sum_closes_the_documented_figures"]["sum"],
            reading["bucket_sum_closes_the_documented_figures"]["documented_eligible_count"],
        )
        # The two readings must be distinguishable, never conflated.
        self.assertNotEqual(
            reading["material_difference"]["basis"],
            reading["hypothesised_difference"]["basis"],
        )

    def test_dry_run_writes_nothing(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self._fixture(tmp)
            before = sorted(
                (str(p.relative_to(root)), p.stat().st_size)
                for p in root.rglob("*")
                if p.is_file()
            )
            self._run_main(["--run-root", str(root), "--dry-run"])
            after = sorted(
                (str(p.relative_to(root)), p.stat().st_size)
                for p in root.rglob("*")
                if p.is_file()
            )
        self.assertEqual(before, after)

    def test_missing_input_returns_one_without_traceback(self):
        with tempfile.TemporaryDirectory() as tmp:
            code, _ = self._run_main(["--run-root", str(tmp), "--dry-run"])
        self.assertEqual(code, 1)

    def test_run_id_with_path_traversal_is_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self._fixture(tmp)
            code, _ = self._run_main(
                ["--run-root", str(root), "--run-id", "../evil", "--output-dir", str(root / "out")]
            )
        self.assertEqual(code, 1)


class ReuseTest(unittest.TestCase):
    """The module must reuse the two existing scripts, not re-implement them."""

    def test_grodon_checker_is_imported_by_identity(self):
        flow = load_module("reconcile_phaded_candidate_flow", REUSE_FLOW)
        module = load_module("reconcile_phaded_grodon_manifest_delta", SCRIPT)
        self.assertIs(module.explain_manifest_difference, load_module(
            "grodon_reanalysis_v2", REUSE_GRODON
        ).explain_manifest_difference)
        self.assertIs(module.check_grodon_manifest_delta, flow.check_grodon_manifest_delta)

    def test_check_grodon_manifest_delta_reaches_ok_on_the_reconciled_sets(self):
        flow = load_module("reconcile_phaded_candidate_flow", REUSE_FLOW)
        eligible = {gid(i) for i in range(1, 11)}
        manifest = {gid(i) for i in range(1, 6)}
        result = flow.check_grodon_manifest_delta(
            eligible_positives=eligible,
            manifest_positives=manifest,
            exclusion_reasons={
                gid(6): "shuffled_out_of_the_1to1_balance",
                gid(7): "shuffled_out_of_the_1to1_balance",
                gid(8): "shuffled_out_of_the_1to1_balance",
                gid(9): "shuffled_out_of_the_1to1_balance",
                gid(10): "shuffled_out_of_the_1to1_balance",
            },
            prediction_failures={},
            expected_eligible_count=10,
            expected_manifest_count=5,
            expected_difference=5,
        )
        self.assertEqual(result["status"], "ok", result.get("problems"))
        self.assertEqual(result["reason_bucket_counts"], {"shuffled_out_of_the_1to1_balance": 5})

    def test_check_grodon_manifest_delta_requires_a_real_eligible_set(self):
        """Without a per-accession eligible set the check stays a mismatch.

        This is the frozen-run situation: only aggregates exist, so the checker
        must not be handed a fabricated eligible list.
        """
        flow = load_module("reconcile_phaded_candidate_flow", REUSE_FLOW)
        manifest = {gid(i) for i in range(1, 6)}
        result = flow.check_grodon_manifest_delta(
            eligible_positives=manifest,
            manifest_positives=manifest,
            exclusion_reasons={},
            prediction_failures={},
            expected_eligible_count=10,
            expected_manifest_count=5,
            expected_difference=5,
        )
        self.assertEqual(result["status"], "mismatch")
        self.assertIn("eligible_positive_count", result["count_expectation_failures"])


class StdlibOnlyTest(unittest.TestCase):
    """The new analysis script is standard library only."""

    def test_source_has_no_third_party_imports(self):
        source = SCRIPT.read_text(encoding="utf-8")
        for banned in ("import numpy", "import pandas", "import scipy", "from numpy", "from pandas"):
            self.assertNotIn(banned, source, banned)

    def test_main_returns_int_and_never_touches_frozen_paths(self):
        source = SCRIPT.read_text(encoding="utf-8")
        self.assertIn("def main(", source)
        self.assertNotIn("deploy/2026", source)


if __name__ == "__main__":
    unittest.main()
