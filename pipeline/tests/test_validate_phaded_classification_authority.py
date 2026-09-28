"""Task 1 invariants: PhaDED superfamily/family layers must not be conflated.

The core scientific invariant is that the 8 PhaDED superfamilies are a
*functional prior* hard framework (Knoll 2009: 28 experimentally validated seed
sequences assigned by function/localisation), while the 38 homologous families
are a *secondary* 2009 sequence-similarity + phylogenetic clustering layer taken
from a database frozen at v1.1 since 2009.  Only the functional-prior layer may
be registry eligible.  Conflating the two layers is the defect this test locks
down.
"""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "pipeline" / "scripts" / "validate_phaded_classification_authority.py"
RUN_ID = "20260917_phaded_classification_authority_01"
RUN_DIR = REPO_ROOT / "runs" / RUN_ID
AUTHORITY_TABLE = RUN_DIR / "inputs" / "phaded_classification_authority.tsv"
VALIDATION_JSON = RUN_DIR / "results" / "authority_validation.json"

FAMILY_DEFINITIONS = (
    REPO_ROOT
    / "runs"
    / "20260909_phaded_reference_evidence_01"
    / "inputs"
    / "phaded_family_definitions.tsv"
)
PROFILE_MANIFEST = (
    REPO_ROOT
    / "runs"
    / "20260915_phaded_profile_gap_audit_01"
    / "inputs"
    / "profile_manifest.tsv"
)

# Independently measured in this run (see status doc); re-measured by the test.
FAMILY_DEFINITIONS_SHA256 = (
    "9c9a21a778390873b6377903f90bb29a400747351081aac570f60b4fcbeea62a"
)
PROFILE_MANIFEST_SHA256 = (
    "50f584eab55516c7874dbedaed459b366f6fbed9b8be63a1844aabbc1a4b42c6"
)

REQUIRED_COLUMNS = [
    "layer_id",
    "layer_kind",
    "source_evidence_type",
    "source_version",
    "registry_eligible",
    "gate_profile_id",
    "notes",
]

EXPECTED_SUPERFAMILY_COUNT = 8
EXPECTED_FAMILY_COUNT = 38
EXPECTED_CANDIDATE_UNIVERSE = 109087


def _load_module():
    spec = importlib.util.spec_from_file_location(
        "validate_phaded_classification_authority", SCRIPT
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class SourceEvidenceTest(unittest.TestCase):
    """Ground the authority table in the two upstream artifacts it binds."""

    def test_family_definitions_hash_matches_independently_measured_value(self):
        self.assertTrue(FAMILY_DEFINITIONS.is_file())
        self.assertEqual(
            sha256_file(FAMILY_DEFINITIONS),
            FAMILY_DEFINITIONS_SHA256,
            "phaded_family_definitions.tsv changed; the authority table is bound "
            "to this exact snapshot hash",
        )

    def test_profile_manifest_hash_is_stable_across_all_copies(self):
        self.assertTrue(PROFILE_MANIFEST.is_file())
        self.assertEqual(sha256_file(PROFILE_MANIFEST), PROFILE_MANIFEST_SHA256)
        copies = sorted(REPO_ROOT.rglob("profile_manifest.tsv"))
        self.assertGreaterEqual(len(copies), 1)
        for copy in copies:
            self.assertEqual(
                sha256_file(copy),
                PROFILE_MANIFEST_SHA256,
                f"a profile_manifest.tsv copy diverges: {copy}",
            )

    def test_family_definitions_declare_8_superfamilies_and_38_families(self):
        rows = read_rows(FAMILY_DEFINITIONS)
        superfamilies = {row["phaded_superfamily"] for row in rows}
        families = {row["phaded_family_id"] for row in rows}
        self.assertEqual(len(superfamilies), EXPECTED_SUPERFAMILY_COUNT)
        self.assertEqual(len(families), EXPECTED_FAMILY_COUNT)
        self.assertEqual(len(rows), EXPECTED_FAMILY_COUNT)


class AuthorityTableTest(unittest.TestCase):
    def setUp(self):
        self.assertTrue(
            AUTHORITY_TABLE.is_file(), f"missing authority table: {AUTHORITY_TABLE}"
        )
        with AUTHORITY_TABLE.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            self.fieldnames = list(reader.fieldnames or [])
            self.rows = list(reader)

    def _by_kind(self, kind: str) -> list[dict[str, str]]:
        return [row for row in self.rows if row["layer_kind"] == kind]

    def test_required_columns_present_in_declared_order(self):
        self.assertEqual(self.fieldnames, REQUIRED_COLUMNS)

    def test_row_counts_are_8_superfamilies_and_38_families(self):
        self.assertEqual(len(self._by_kind("superfamily")), EXPECTED_SUPERFAMILY_COUNT)
        self.assertEqual(len(self._by_kind("family")), EXPECTED_FAMILY_COUNT)
        self.assertEqual(len(self.rows), 46)

    def test_layer_ids_are_unique_and_non_empty(self):
        ids = [row["layer_id"] for row in self.rows]
        self.assertTrue(all(ids), "layer_id must never be empty")
        self.assertEqual(len(set(ids)), len(ids), "layer_id must be unique")

    # --- core invariant #1: the hard functional-prior framework -----------
    def test_all_superfamilies_are_functional_prior_and_registry_eligible(self):
        rows = self._by_kind("superfamily")
        self.assertEqual(len(rows), EXPECTED_SUPERFAMILY_COUNT)
        for row in rows:
            self.assertEqual(
                row["source_evidence_type"],
                "functional_prior",
                f"{row['layer_id']}: superfamilies come from the function/"
                "localisation prior on the 28 experimentally validated seeds",
            )
            self.assertEqual(
                row["registry_eligible"],
                "true",
                f"{row['layer_id']}: the functional-prior layer is registry eligible",
            )

    # --- core invariant #2: the frozen 2009 clustering refinement layer ---
    def test_all_families_are_sequence_clustering_2009_frozen_and_ineligible(self):
        rows = self._by_kind("family")
        self.assertEqual(len(rows), EXPECTED_FAMILY_COUNT)
        for row in rows:
            self.assertEqual(
                row["source_evidence_type"],
                "sequence_clustering_2009",
                f"{row['layer_id']}: families are a 2009 sequence-similarity + "
                "phylogenetic clustering layer, not a functional prior",
            )
            self.assertEqual(
                row["source_version"],
                "Knoll2009_v1.1_frozen",
                f"{row['layer_id']}: family layer must be bound to the frozen "
                "Knoll 2009 v1.1 snapshot",
            )
            self.assertEqual(
                row["registry_eligible"],
                "false",
                f"{row['layer_id']}: THIS IS THE CORE INVARIANT — a 2009 "
                "clustering refinement layer must never be registry eligible",
            )

    def test_evidence_type_and_registry_eligibility_are_never_crossed(self):
        for row in self.rows:
            pair = (row["source_evidence_type"], row["registry_eligible"])
            self.assertIn(
                pair,
                {("functional_prior", "true"), ("sequence_clustering_2009", "false")},
                f"{row['layer_id']}: layer mismatch — {pair} crosses the "
                "superfamily/family hierarchy",
            )

    def test_every_row_is_bound_to_the_frozen_ded_snapshot(self):
        for row in self.rows:
            self.assertEqual(row["source_version"], "Knoll2009_v1.1_frozen")
            self.assertTrue(row["notes"].strip(), f"{row['layer_id']}: empty notes")

    # --- core invariant #3: exact correspondence with the DED import ------
    def test_layer_ids_correspond_one_to_one_with_family_definitions(self):
        definitions = read_rows(FAMILY_DEFINITIONS)
        expected_superfamilies = {row["phaded_superfamily"] for row in definitions}
        expected_families = {row["phaded_family_id"] for row in definitions}
        observed_superfamilies = {
            row["layer_id"] for row in self._by_kind("superfamily")
        }
        observed_families = {row["layer_id"] for row in self._by_kind("family")}
        self.assertEqual(observed_superfamilies - expected_superfamilies, set())
        self.assertEqual(expected_superfamilies - observed_superfamilies, set())
        self.assertEqual(observed_families - expected_families, set())
        self.assertEqual(expected_families - observed_families, set())

    def test_gate_profile_ids_resolve_to_the_matching_profile_kind(self):
        manifest = {
            row["profile_id"]: row for row in read_rows(PROFILE_MANIFEST)
        }
        for row in self.rows:
            profile_id = row["gate_profile_id"]
            self.assertIn(profile_id, manifest, f"{row['layer_id']}: unresolvable")
            self.assertEqual(manifest[profile_id]["profile_kind"], row["layer_kind"])

    def test_family_gate_profiles_confirm_the_superfamily_parentage(self):
        definitions = {
            row["phaded_family_id"]: row["phaded_superfamily"]
            for row in read_rows(FAMILY_DEFINITIONS)
        }
        manifest = {row["profile_id"]: row for row in read_rows(PROFILE_MANIFEST)}
        superfamily_ids = {row["layer_id"] for row in self._by_kind("superfamily")}
        for row in self._by_kind("family"):
            parent = manifest[row["gate_profile_id"]]["phaded_superfamily"]
            self.assertEqual(parent, definitions[row["layer_id"]])
            self.assertIn(parent, superfamily_ids)


class FailClosedTest(unittest.TestCase):
    """The validator must reject a conflated table, not merely describe it."""

    def setUp(self):
        self.module = _load_module()
        self.rows = read_rows(AUTHORITY_TABLE)

    def _validate(self, rows):
        return self.module.validate_authority_rows(
            rows,
            family_definitions_path=FAMILY_DEFINITIONS,
            profile_manifest_path=PROFILE_MANIFEST,
        )

    def test_pristine_table_passes(self):
        report = self._validate(self.rows)
        self.assertEqual(report["status"], "verified")

    def test_family_promoted_to_registry_eligible_is_rejected(self):
        rows = [dict(row) for row in self.rows]
        victim = next(row for row in rows if row["layer_kind"] == "family")
        victim["registry_eligible"] = "true"
        with self.assertRaises(self.module.AuthorityValidationError) as ctx:
            self._validate(rows)
        self.assertIn(victim["layer_id"], str(ctx.exception))

    def test_superfamily_demoted_to_clustering_layer_is_rejected(self):
        rows = [dict(row) for row in self.rows]
        victim = next(row for row in rows if row["layer_kind"] == "superfamily")
        victim["source_evidence_type"] = "sequence_clustering_2009"
        victim["registry_eligible"] = "false"
        with self.assertRaises(self.module.AuthorityValidationError):
            self._validate(rows)

    def test_wrong_source_version_is_rejected(self):
        rows = [dict(row) for row in self.rows]
        rows[0]["source_version"] = "Knoll2009_live"
        with self.assertRaises(self.module.AuthorityValidationError):
            self._validate(rows)

    def test_extra_layer_id_is_rejected(self):
        rows = [dict(row) for row in self.rows]
        extra = dict(rows[0])
        extra["layer_id"] = "DED_hfam_999"
        extra["layer_kind"] = "family"
        extra["source_evidence_type"] = "sequence_clustering_2009"
        extra["source_version"] = "Knoll2009_v1.1_frozen"
        extra["registry_eligible"] = "false"
        rows.append(extra)
        with self.assertRaises(self.module.AuthorityValidationError):
            self._validate(rows)

    def test_missing_layer_id_is_rejected(self):
        rows = [dict(row) for row in self.rows]
        victim = next(row for row in rows if row["layer_kind"] == "family")
        rows.remove(victim)
        with self.assertRaises(self.module.AuthorityValidationError):
            self._validate(rows)

    def test_missing_required_column_is_rejected(self):
        rows = [{k: v for k, v in row.items() if k != "notes"} for row in self.rows]
        with self.assertRaises(self.module.AuthorityValidationError):
            self._validate(rows)

    def test_empty_registry_eligible_value_is_rejected(self):
        rows = [dict(row) for row in self.rows]
        rows[0]["registry_eligible"] = ""
        with self.assertRaises(self.module.AuthorityValidationError):
            self._validate(rows)


class GeneratorTest(unittest.TestCase):
    """Generation must be deterministic and reproduce the committed table."""

    def setUp(self):
        self.module = _load_module()

    def test_generated_rows_reproduce_the_committed_table(self):
        generated = self.module.build_authority_rows(
            family_definitions_path=FAMILY_DEFINITIONS,
            profile_manifest_path=PROFILE_MANIFEST,
        )
        self.assertEqual(generated, read_rows(AUTHORITY_TABLE))

    def test_generation_is_deterministic(self):
        first = self.module.build_authority_rows(
            family_definitions_path=FAMILY_DEFINITIONS,
            profile_manifest_path=PROFILE_MANIFEST,
        )
        second = self.module.build_authority_rows(
            family_definitions_path=FAMILY_DEFINITIONS,
            profile_manifest_path=PROFILE_MANIFEST,
        )
        self.assertEqual(first, second)

    def test_generator_refuses_an_unparseable_provenance_manifest(self):
        with tempfile.TemporaryDirectory() as tmp:
            broken = Path(tmp) / "profile_manifest.tsv"
            broken.write_text("profile_id\tprofile_kind\n", encoding="utf-8")
            with self.assertRaises(self.module.AuthorityValidationError):
                self.module.build_authority_rows(
                    family_definitions_path=FAMILY_DEFINITIONS,
                    profile_manifest_path=broken,
                )


class RunArtifactTest(unittest.TestCase):
    def test_run_layout_and_contract_exist(self):
        for name in ("logs", "inputs", "results"):
            self.assertTrue((RUN_DIR / name).is_dir(), f"missing {name}/")
        contract_path = RUN_DIR / "input_contract.json"
        self.assertTrue(contract_path.is_file())
        contract = json.loads(contract_path.read_text(encoding="utf-8"))
        self.assertEqual(contract["run_id"], RUN_ID)
        recorded = contract["inputs"]["phaded_classification_authority"]
        self.assertEqual(recorded["status"], "verified")
        self.assertEqual(recorded["sha256"], sha256_file(AUTHORITY_TABLE))
        self.assertEqual(recorded["size"], AUTHORITY_TABLE.stat().st_size)

    def test_validation_json_declares_verified_status_and_counts(self):
        report = json.loads(VALIDATION_JSON.read_text(encoding="utf-8"))
        self.assertEqual(report["status"], "verified")
        self.assertEqual(
            report["layer_counts"],
            {"superfamily": 8, "family": 38, "total": 46},
        )
        self.assertEqual(
            report["registry_eligible_counts"], {"true": 8, "false": 38}
        )
        self.assertTrue(
            all(item["passed"] for item in report["invariants"]),
            "every recorded invariant must have passed",
        )

    def test_validation_json_declares_zero_candidate_impact(self):
        report = json.loads(VALIDATION_JSON.read_text(encoding="utf-8"))
        impact = report["candidate_universe_impact"]
        self.assertEqual(impact["candidates_total"], EXPECTED_CANDIDATE_UNIVERSE)
        self.assertEqual(impact["rows_changed"], 0)
        self.assertIn("candidate-only", report["boundary_statement"])

    def test_validation_json_records_every_source_file_hash(self):
        report = json.loads(VALIDATION_JSON.read_text(encoding="utf-8"))
        recorded = {item["role"]: item for item in report["source_files"]}
        for role in (
            "phaded_family_definitions",
            "profile_manifest",
            "phaded_classification_authority",
        ):
            self.assertIn(role, recorded)
            self.assertEqual(recorded[role]["status"], "verified")
            self.assertRegex(recorded[role]["sha256"], r"^[0-9a-f]{64}$")
        self.assertEqual(
            recorded["phaded_family_definitions"]["sha256"],
            FAMILY_DEFINITIONS_SHA256,
        )


if __name__ == "__main__":
    unittest.main()
