"""The v2 calibration readiness must be append-only over the frozen readiness.

``build_phaded_calibration_readiness_v2.py`` (F15 of
``docs/superpowers/plans/2026-09-28-phaded-evidence-model-redesign-followup.md``)
starts from the frozen amendment readiness TSV and adds the two qualification
columns the 2026-09-28 redesign (plan Task 5) requires:

* ``qualified_e2_e3_positive_count`` -- unique ``independence_group`` count of
  ledger rows graded E2/E3 with ``functional_calibration_eligible=true``;
* ``independent_genus_count`` -- distinct genera among those groups.

These tests pin the properties the run depends on, before the builder exists:

* **append-only**: every frozen column of every frozen row survives verbatim, in
  the frozen row order, as a tab prefix of the emitted line, so the frozen
  byte-level regression on the untouched input can never move;
* the two new columns are derived *only* from the curated ledger, and the counts
  reproduce ``reference_evidence_independence_summary.json`` exactly when that
  summary is supplied;
* the profile frame is closed in both directions: a ledger family without a
  readiness profile, or a readiness profile without a ledger family, fails
  closed instead of being dropped;
* rows that the frozen readiness does not cover are materialised from the frozen
  profile manifest plus cited trained-profile evidence -- never invented -- and
  an unresolvable or ambiguous profile fails closed;
* the promotion scope is explicit: a deferred profile is named on the command
  line with a reason and is marked in the readiness itself.
"""

import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "pipeline" / "scripts" / "build_phaded_calibration_readiness_v2.py"

FROZEN_FIELDS = (
    "profile_id", "profile_kind", "phaded_superfamily", "phaded_family_id",
    "existing_experimental_positive", "external_bound_positive",
    "external_bound_negative", "external_bound_challenge", "heldout_positive",
    "family_resolved_negative_status", "challenge_status", "calibration_decision",
    "new_family_call_created",
)
LEDGER_FIELDS = (
    "reference_id", "accession", "organism", "phaded_superfamily", "phaded_family_id",
    "independence_group", "experimental_evidence_grade", "functional_calibration_eligible",
)
MANIFEST_FIELDS = (
    "profile_id", "profile_kind", "phaded_superfamily", "phaded_family_id",
    "training_accessions", "training_count", "model_status",
)
TRAINED_FIELDS = (
    "profile_id", "external_bound_negative", "external_bound_negative_accessions",
    "family_resolved_negative_status", "calibration_decision", "evidence_citation",
)
HELDOUT_FIELDS = ("family", "heldout", "evalue", "pass", "note")
CHALLENGE_FIELDS = ("family", "challenge", "evalue", "rejected")


def load_module():
    spec = importlib.util.spec_from_file_location("build_phaded_calibration_readiness_v2", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_tsv(path: Path, fields, rows) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})
    return path


def read_tsv(path: Path):
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def frozen_row(profile_id, kind, superfamily, family_id, **overrides):
    row = {
        "profile_id": profile_id, "profile_kind": kind, "phaded_superfamily": superfamily,
        "phaded_family_id": family_id, "existing_experimental_positive": "1",
        "external_bound_positive": "0", "external_bound_negative": "0",
        "external_bound_challenge": "4", "heldout_positive": "0",
        "family_resolved_negative_status": "missing", "challenge_status": "sufficient",
        "calibration_decision": "reference_only_insufficient_panel",
        "new_family_call_created": "false",
    }
    row.update(overrides)
    return row


def ledger_row(family, group, grade, eligible, organism):
    return {
        "reference_id": f"{family}_{group[:6]}", "accession": "X.1", "organism": organism,
        "phaded_superfamily": "sf", "phaded_family_id": family,
        "independence_group": group, "experimental_evidence_grade": grade,
        "functional_calibration_eligible": eligible,
    }


def manifest_row(profile_id, kind, superfamily, family_id, count, status="trained"):
    return {
        "profile_id": profile_id, "profile_kind": kind, "phaded_superfamily": superfamily,
        "phaded_family_id": family_id, "training_accessions": "X.1", "training_count": str(count),
        "model_status": status,
    }


class Fixture:
    """Synthetic project frame: 3 families (A/B/C) + 2 superfamilies (SF1/SF2)."""

    def __init__(self, root: Path):
        self.root = root
        self.frozen = write_tsv(root / "frozen_readiness.tsv", FROZEN_FIELDS, [
            frozen_row("family_DED_hfam_901_aaaa", "family", "sf1", "DED_hfam_901",
                       existing_experimental_positive="9"),
            frozen_row("superfamily_sf1_bbbb", "superfamily", "sf1", "",
                       existing_experimental_positive="0", external_bound_challenge="0",
                       challenge_status="missing",
                       calibration_decision="planned_not_run_superfamily_requires_family_resolution"),
        ])
        self.ledger = write_tsv(root / "ledger.tsv", LEDGER_FIELDS, [
            ledger_row("DED_hfam_901", "IG_9011", "E2", "true", "Bacillus subtilis"),
            ledger_row("DED_hfam_901", "IG_9011", "E2", "true", "Bacillus subtilis"),
            ledger_row("DED_hfam_902", "IG_9021", "E3", "true", "Comamonas testosteroni"),
            ledger_row("DED_hfam_902", "IG_9022", "E2", "true", "Delftia acidovorans"),
            ledger_row("DED_hfam_902", "IG_9023", "E2", "true", "Ralstonia pickettii"),
            ledger_row("DED_hfam_902", "IG_9024", "E1", "true", "Pseudomonas putida"),
            ledger_row("DED_hfam_902", "IG_9025", "E2", "false", "Cupriavidus necator"),
            ledger_row("DED_hfam_903", "IG_9031", "E2", "true", "Alcaligenes faecalis"),
            ledger_row("DED_hfam_903", "IG_9032", "A", "false", "Bacillus cereus"),
        ])
        self.manifest = write_tsv(root / "profile_manifest.tsv", MANIFEST_FIELDS, [
            manifest_row("family_DED_hfam_901_aaaa", "family", "sf1", "DED_hfam_901", 1, "reference_only"),
            manifest_row("family_DED_hfam_902_bbbb", "family", "sf2", "DED_hfam_902", 7),
            manifest_row("family_DED_hfam_903_cccc", "family", "sf2", "DED_hfam_903", 3),
            manifest_row("superfamily_sf1_bbbb", "superfamily", "sf1", "", 1, "reference_only"),
            manifest_row("superfamily_sf2_dddd", "superfamily", "sf2", "", 7),
        ])
        self.trained = write_tsv(root / "trained_evidence.tsv", TRAINED_FIELDS, [
            {"profile_id": "family_DED_hfam_902_bbbb", "external_bound_negative": "2",
             "external_bound_negative_accessions": "Q84C08;Q51718",
             "family_resolved_negative_status": "sufficient",
             "calibration_decision": "candidate_gate_passed_not_promoted",
             "evidence_citation": "runs/x/results/gate_determination.json::gate_legs.negative_count"},
            {"profile_id": "family_DED_hfam_903_cccc", "external_bound_negative": "1",
             "external_bound_negative_accessions": "Q88N36",
             "family_resolved_negative_status": "sufficient",
             "calibration_decision": "reference_only_insufficient_panel",
             "evidence_citation": "runs/x/results/negatives.json::family_resolved_negatives[0]"},
        ])
        self.heldout = write_tsv(root / "heldout_validation.tsv", HELDOUT_FIELDS, [
            {"family": "hfam_902", "heldout": "B1", "evalue": "1e-50", "pass": "pass", "note": ""},
            {"family": "hfam_903", "heldout": "C1", "evalue": "", "pass": "leakage_blocked", "note": "max_id=0.9"},
        ])
        self.challenge = write_tsv(root / "challenge_validation.tsv", CHALLENGE_FIELDS, [
            {"family": "hfam_902", "challenge": "P1", "evalue": "999", "rejected": "yes"},
            {"family": "hfam_902", "challenge": "P2", "evalue": "999", "rejected": "yes"},
            {"family": "hfam_902", "challenge": "P3", "evalue": "999", "rejected": "yes"},
            {"family": "hfam_903", "challenge": "P1", "evalue": "1e-9", "rejected": "NO"},
        ])
        self.summary = root / "independence_summary.json"
        self.summary.write_text(json.dumps({
            "family_functional_positive_counts": {
                "DED_hfam_901": 1, "DED_hfam_902": 3, "DED_hfam_903": 1,
            },
            "family_distinct_genus_counts": {
                "DED_hfam_901": 1, "DED_hfam_902": 3, "DED_hfam_903": 1,
            },
        }), encoding="utf-8")
        self.out = root / "out"

    def build(self, module, **kwargs):
        options = dict(
            frozen_readiness=self.frozen, curated_ledger=self.ledger,
            profile_manifest=self.manifest, trained_evidence=self.trained,
            heldout_validation=self.heldout, challenge_validation=self.challenge,
            independence_summary=self.summary,
        )
        options.update(kwargs)
        output_dir = options.pop("output_dir", self.out)
        run_id = options.pop("run_id", "test")
        defer_profiles = options.pop("defer_profiles", ())
        defer_reason = options.pop("defer_reason", "")
        return module.build(output_dir=output_dir, run_id=run_id,
                            defer_profiles=defer_profiles, defer_reason=defer_reason, **options)


class AppendOnlyTests(unittest.TestCase):
    def test_frozen_columns_survive_verbatim_in_frozen_order(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture = Fixture(Path(tmp))
            module = load_module()
            fixture.build(module)
            frozen_lines = fixture.frozen.read_text(encoding="utf-8").splitlines()
            emitted = (fixture.out / "calibration_readiness_v2.tsv").read_text(encoding="utf-8").splitlines()
            self.assertEqual(emitted[0].split("\t")[:len(FROZEN_FIELDS)], list(FROZEN_FIELDS))
            for index, line in enumerate(frozen_lines):
                self.assertEqual(emitted[index].split("\t")[:len(FROZEN_FIELDS)], line.split("\t"),
                                 f"frozen line {index} moved: append-only violated")

    def test_frozen_values_are_not_recomputed_from_the_manifest(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture = Fixture(Path(tmp))
            module = load_module()
            fixture.build(module)
            fields, rows = read_tsv(fixture.out / "calibration_readiness_v2.tsv")
            by_id = {row["profile_id"]: row for row in rows}
            row = by_id["family_DED_hfam_901_aaaa"]
            # frozen 9 beats the manifest's training_count of 1: the frozen value is authoritative
            self.assertEqual(row["existing_experimental_positive"], "9")
            self.assertEqual(row["calibration_decision"], "reference_only_insufficient_panel")

    def test_materialised_rows_follow_the_frozen_rows(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture = Fixture(Path(tmp))
            module = load_module()
            fixture.build(module)
            fields, rows = read_tsv(fixture.out / "calibration_readiness_v2.tsv")
            self.assertEqual(len(rows), 5)
            self.assertEqual([row["profile_id"] for row in rows[:2]],
                             ["family_DED_hfam_901_aaaa", "superfamily_sf1_bbbb"])
            self.assertEqual([row["profile_id"] for row in rows[2:]],
                             ["family_DED_hfam_902_bbbb", "family_DED_hfam_903_cccc", "superfamily_sf2_dddd"])


class QualificationColumnTests(unittest.TestCase):
    def test_counts_come_from_the_ledger_only(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture = Fixture(Path(tmp))
            module = load_module()
            fixture.build(module)
            fields, rows = read_tsv(fixture.out / "calibration_readiness_v2.tsv")
            counts = {row["phaded_family_id"]: (row["qualified_e2_e3_positive_count"],
                                                row["independent_genus_count"]) for row in rows}
            # A: 1 unique E2 group; B: IG_9024 is E1 and IG_9025 is not calibration eligible
            self.assertEqual(counts["DED_hfam_901"], ("1", "1"))
            self.assertEqual(counts["DED_hfam_902"], ("3", "3"))
            self.assertEqual(counts["DED_hfam_903"], ("1", "1"))
            self.assertEqual(counts[""], ("0", "0"))  # superfamily rows carry no family resolution

    def test_qualification_source_names_file_and_columns(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture = Fixture(Path(tmp))
            module = load_module()
            fixture.build(module)
            fields, rows = read_tsv(fixture.out / "calibration_readiness_v2.tsv")
            source = rows[0]["qualification_source"]
            self.assertIn("ledger.tsv::", source)
            for column in ("experimental_evidence_grade", "functional_calibration_eligible",
                           "independence_group", "organism"):
                self.assertIn(column, source)

    def test_mapping_table_names_groups_genera_and_source(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture = Fixture(Path(tmp))
            module = load_module()
            fixture.build(module)
            fields, rows = read_tsv(fixture.out / "family_qualification_mapping.tsv")
            by_family = {row["phaded_family_id"]: row for row in rows}
            self.assertEqual(set(by_family), {"DED_hfam_901", "DED_hfam_902", "DED_hfam_903"})
            row = by_family["DED_hfam_902"]
            self.assertEqual(row["qualified_e2_e3_positive_count"], "3")
            self.assertEqual(row["independent_genus_count"], "3")
            self.assertEqual(row["contributing_independence_groups"], "IG_9021;IG_9022;IG_9023")
            self.assertEqual(row["contributing_genera"], "Comamonas;Delftia;Ralstonia")
            self.assertIn("independence_group", row["source_column"])
            self.assertIn("ledger.tsv", row["source_file"])

    def test_independence_summary_is_cross_checked(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture = Fixture(Path(tmp))
            fixture.summary.write_text(json.dumps({
                "family_functional_positive_counts": {"DED_hfam_901": 1, "DED_hfam_902": 9, "DED_hfam_903": 1},
                "family_distinct_genus_counts": {"DED_hfam_901": 1, "DED_hfam_902": 9, "DED_hfam_903": 1},
            }), encoding="utf-8")
            module = load_module()
            with self.assertRaises(module.ReadinessReconciliationError):
                fixture.build(module)


class ReconciliationFailClosedTests(unittest.TestCase):
    def test_ledger_family_without_profile_fails_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture = Fixture(Path(tmp))
            fields, rows = read_tsv(fixture.ledger)
            rows.append(ledger_row("DED_hfam_904", "IG_9041", "E2", "true", "Bacillus subtilis"))
            write_tsv(fixture.ledger, fields, rows)
            module = load_module()
            with self.assertRaises(module.ReadinessReconciliationError) as error:
                fixture.build(module)
            self.assertIn("DED_hfam_904", str(error.exception))

    def test_profile_without_ledger_family_fails_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture = Fixture(Path(tmp))
            fields, rows = read_tsv(fixture.manifest)
            rows.append(manifest_row("family_DED_hfam_E_eeee", "family", "sf2", "DED_hfam_905", 1))
            write_tsv(fixture.manifest, fields, rows)
            module = load_module()
            with self.assertRaises(module.ReadinessReconciliationError) as error:
                fixture.build(module)
            self.assertIn("DED_hfam_905", str(error.exception))

    def test_unresolvable_materialised_profile_fails_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture = Fixture(Path(tmp))
            fields, rows = read_tsv(fixture.trained)
            write_tsv(fixture.trained, fields, [row for row in rows
                                                if row["profile_id"] != "family_DED_hfam_903_cccc"])
            module = load_module()
            with self.assertRaises(module.ReadinessReconciliationError) as error:
                fixture.build(module)
            self.assertIn("family_DED_hfam_903_cccc", str(error.exception))

    def test_profile_present_twice_fails_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture = Fixture(Path(tmp))
            fields, rows = read_tsv(fixture.trained)
            rows.append(dict(rows[0]))
            write_tsv(fixture.trained, fields, rows)
            module = load_module()
            with self.assertRaises(module.ReadinessReconciliationError):
                fixture.build(module)

    def test_profile_in_both_frozen_and_materialised_fails_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture = Fixture(Path(tmp))
            fields, rows = read_tsv(fixture.trained)
            rows.append({"profile_id": "family_DED_hfam_901_aaaa", "external_bound_negative": "1",
                         "external_bound_negative_accessions": "X",
                         "family_resolved_negative_status": "sufficient",
                         "calibration_decision": "candidate_gate_passed_not_promoted",
                         "evidence_citation": "runs/x::y"})
            write_tsv(fixture.trained, fields, rows)
            module = load_module()
            with self.assertRaises(module.ReadinessReconciliationError):
                fixture.build(module)

    def test_duplicate_profile_id_in_manifest_fails_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture = Fixture(Path(tmp))
            fields, rows = read_tsv(fixture.manifest)
            rows.append(dict(rows[0]))
            write_tsv(fixture.manifest, fields, rows)
            module = load_module()
            with self.assertRaises(module.ReadinessReconciliationError):
                fixture.build(module)

    def test_unknown_grade_or_eligibility_token_fails_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture = Fixture(Path(tmp))
            fields, rows = read_tsv(fixture.ledger)
            rows.append(ledger_row("DED_hfam_901", "IG_9019", "E2", "maybe", "Bacillus subtilis"))
            write_tsv(fixture.ledger, fields, rows)
            module = load_module()
            with self.assertRaises(module.ReadinessReconciliationError):
                fixture.build(module)

    def test_missing_frozen_column_fails_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture = Fixture(Path(tmp))
            fields, rows = read_tsv(fixture.frozen)
            write_tsv(fixture.frozen, [f for f in fields if f != "challenge_status"], rows)
            module = load_module()
            with self.assertRaises(module.ReadinessReconciliationError):
                fixture.build(module)

    def test_missing_heldout_row_fails_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture = Fixture(Path(tmp))
            fields, rows = read_tsv(fixture.heldout)
            write_tsv(fixture.heldout, fields, [row for row in rows if row["family"] != "hfam_903"])
            module = load_module()
            with self.assertRaises(module.ReadinessReconciliationError) as error:
                fixture.build(module)
            self.assertIn("DED_hfam_903", str(error.exception))

    def test_defer_without_reason_fails_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture = Fixture(Path(tmp))
            module = load_module()
            with self.assertRaises(module.ReadinessReconciliationError):
                fixture.build(module, defer_profiles=("family_DED_hfam_902_bbbb",))

    def test_defer_unknown_profile_fails_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture = Fixture(Path(tmp))
            module = load_module()
            with self.assertRaises(module.ReadinessReconciliationError):
                fixture.build(module, defer_profiles=("family_DED_hfam_999_zzzz",),
                              defer_reason="because")

    def test_existing_output_files_fail_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture = Fixture(Path(tmp))
            module = load_module()
            fixture.out.mkdir(parents=True, exist_ok=True)
            (fixture.out / "calibration_readiness_v2.tsv").write_text("x", encoding="utf-8")
            with self.assertRaises(FileExistsError):
                fixture.build(module)


class PromotionScopeTests(unittest.TestCase):
    def test_default_scope_is_the_authorized_frame(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture = Fixture(Path(tmp))
            module = load_module()
            fixture.build(module)
            fields, rows = read_tsv(fixture.out / "calibration_readiness_v2.tsv")
            self.assertEqual({row["promotion_scope"] for row in rows},
                             {module.SCOPE_AUTHORIZED})

    def test_deferred_profile_is_marked_with_its_reason(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture = Fixture(Path(tmp))
            module = load_module()
            result = fixture.build(module, defer_profiles=("family_DED_hfam_902_bbbb",),
                                   defer_reason="design spec special family: weak records graded separately")
            fields, rows = read_tsv(fixture.out / "calibration_readiness_v2.tsv")
            scope = {row["profile_id"]: row["promotion_scope"] for row in rows}
            self.assertEqual(scope["family_DED_hfam_902_bbbb"], module.SCOPE_DEFERRED)
            self.assertEqual(scope["family_DED_hfam_903_cccc"], module.SCOPE_AUTHORIZED)
            reconciliation = json.loads(
                (fixture.out / "readiness_reconciliation.json").read_text(encoding="utf-8"))
            self.assertEqual([row["profile_id"] for row in reconciliation["promotion_scope"]["deferred"]],
                             ["family_DED_hfam_902_bbbb"])
            self.assertEqual(reconciliation["promotion_scope"]["deferred"][0]["reason"],
                             "design spec special family: weak records graded separately")
            self.assertEqual(result["promotion_scope"]["authorized_count"], 4)


class ReconciliationReportTests(unittest.TestCase):
    def test_report_accounts_for_every_profile_and_ledger_family(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture = Fixture(Path(tmp))
            module = load_module()
            fixture.build(module)
            report = json.loads(
                (fixture.out / "readiness_reconciliation.json").read_text(encoding="utf-8"))
            self.assertEqual(report["counts"]["frame_profiles"], 5)
            self.assertEqual(report["counts"]["frozen_readiness_rows"], 2)
            self.assertEqual(report["counts"]["materialised_rows"], 3)
            self.assertEqual(report["counts"]["ledger_families"], 3)
            self.assertEqual(report["coverage"]["ledger_families_without_profile"], [])
            self.assertEqual(report["coverage"]["profiles_without_ledger_family"], ["superfamily_sf1_bbbb",
                                                                                   "superfamily_sf2_dddd"])
            self.assertTrue(report["inputs"]["frozen_readiness"]["sha256"])

    def test_materialised_family_rows_carry_their_leg_provenance(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture = Fixture(Path(tmp))
            module = load_module()
            fixture.build(module)
            fields, rows = read_tsv(fixture.out / "calibration_readiness_v2.tsv")
            row = {r["profile_id"]: r for r in rows}["family_DED_hfam_902_bbbb"]
            self.assertEqual(row["readiness_row_origin"], "materialised_frozen_manifest_and_trained_evidence")
            self.assertEqual(row["existing_experimental_positive"], "7")  # manifest training_count
            self.assertEqual(row["heldout_positive"], "1")               # heldout_validation pass
            self.assertEqual(row["external_bound_challenge"], "3")       # challenge_validation rows
            self.assertEqual(row["challenge_status"], "sufficient")
            self.assertEqual(row["external_bound_negative"], "2")
            self.assertEqual(row["new_family_call_created"], "false")
            self.assertIn("B1", row["readiness_notes"])

    def test_materialised_superfamily_row_uses_the_frozen_superfamily_convention(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture = Fixture(Path(tmp))
            module = load_module()
            fixture.build(module)
            fields, rows = read_tsv(fixture.out / "calibration_readiness_v2.tsv")
            row = {r["profile_id"]: r for r in rows}["superfamily_sf2_dddd"]
            self.assertEqual(row["external_bound_negative"], "0")
            self.assertEqual(row["family_resolved_negative_status"], "missing")
            self.assertEqual(row["challenge_status"], "missing")
            self.assertEqual(row["calibration_decision"],
                             "planned_not_run_superfamily_requires_family_resolution")
            self.assertEqual(row["readiness_row_origin"],
                             "materialised_frozen_manifest_superfamily_resolution_not_run")

    def test_challenge_failure_is_not_laundered_into_sufficient(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture = Fixture(Path(tmp))
            module = load_module()
            fixture.build(module)
            fields, rows = read_tsv(fixture.out / "calibration_readiness_v2.tsv")
            row = {r["profile_id"]: r for r in rows}["family_DED_hfam_903_cccc"]
            self.assertEqual(row["challenge_status"], "not_sufficient")
            self.assertEqual(row["heldout_positive"], "0")


class DeterminismTests(unittest.TestCase):
    def test_two_builds_are_byte_identical(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture = Fixture(Path(tmp))
            module = load_module()
            fixture.build(module, output_dir=Path(tmp) / "out_a")
            fixture.build(module, output_dir=Path(tmp) / "out_b")
            for name in ("calibration_readiness_v2.tsv", "family_qualification_mapping.tsv"):
                self.assertEqual((Path(tmp) / "out_a" / name).read_bytes(),
                                 (Path(tmp) / "out_b" / name).read_bytes(), name)


if __name__ == "__main__":
    unittest.main()
