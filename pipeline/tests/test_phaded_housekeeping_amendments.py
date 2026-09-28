"""Tests for the Task 0 housekeeping amendment runner.

These tests pin three invariants:

* the provenance amendment records the real 09-16 server execution facts while
  keeping unmeasured server resources explicitly ``pending`` (never ``0``);
* the ambiguity-count discrepancy is explained by layer precedence in
  ``build_phaded_subtype_matrix.build_matrix`` and fully accounted from data;
* every historical input is read-only: no historical manifest is rewritten and
  the recorded SHA-256 of each read file is stable across the read.
"""

from __future__ import annotations

import csv
import hashlib
import json
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import run_phaded_housekeeping_amendments as hk  # noqa: E402


RUN_ID = "20260917_phaded_housekeeping_amendments_01"

MATRIX_HEADER = [
    "accession",
    "profile_evidence_status",
    "subtype_call",
    "subtype_confidence",
    "domain_evidence_status",
    "manual_review_status",
    "assignment_status",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_matrix(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.DictWriter(handle, fieldnames=MATRIX_HEADER, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def matrix_row(accession: str, profile_status: str, subtype_call: str, *, domain: str = "domain_supported_partial", manual: str = "not_in_19_candidate_review", confidence: str = "unresolved", assignment: str = "ambiguous_family") -> dict[str, str]:
    return {
        "accession": accession,
        "profile_evidence_status": profile_status,
        "subtype_call": subtype_call,
        "subtype_confidence": confidence,
        "domain_evidence_status": domain,
        "manual_review_status": manual,
        "assignment_status": assignment,
    }


class HousekeepingFixture(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.root = Path(self._tmp.name)
        self.run_dir = self.root / "runs" / RUN_ID
        for name in ("", "logs", "inputs", "results"):
            (self.run_dir / name).mkdir(parents=True, exist_ok=True)

        self.signalp = self.root / "signalp_run_summary.json"
        self.signalp.write_text(
            json.dumps({"status": "completed", "returncode": 0, "command": ["signalp6", "--mode", "fast"]}),
            encoding="utf-8",
        )
        self.foldseek = self.root / "foldseek_merge_manifest.json"
        self.foldseek.write_text(json.dumps({"status": "completed_candidate_only"}), encoding="utf-8")
        self.matrix = self.root / "foldseek_subtype_matrix.tsv"
        write_matrix(self.matrix, [matrix_row("acc1", "profile_trained_hit", "dPHASCL1_like_candidate")])

        self.historical_contract = self.root / "runs" / "20260916_phaded_full_library_signalp_01"
        self.historical_contract.mkdir(parents=True)
        self.historical_contract_json = self.historical_contract / "input_contract.json"
        self.historical_contract_json.write_text(
            json.dumps(
                {
                    "run_id": "20260916_phaded_full_library_signalp_01",
                    "generated_at": "2026-09-16",
                    "status": "completed_candidate_only",
                    "authorization": {
                        "candidate_only_execution": True,
                        "formal_scan_authorized": False,
                        "formal_registry_modified": False,
                        "server_execution_started": False,
                    },
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )

    def provenance(self, **overrides):
        kwargs = {
            "run_id": RUN_ID,
            "run_dir": self.run_dir,
            "signalp_summary": self.signalp,
            "foldseek_manifest": self.foldseek,
            "final_matrix": self.matrix,
            "server": "<SERVER_HOST>",
            "account": "<SERVER_USER>",
            "signalp_threads": 32,
        }
        kwargs.update(overrides)
        return hk.run_provenance_amendment(**kwargs)


class ProvenanceAmendmentTests(HousekeepingFixture):
    def test_records_server_dated_deploy_without_inventing_server_resources(self) -> None:
        manifest = self.provenance()

        self.assertEqual(manifest["execution"]["mode"], "server_dated_deploy")
        self.assertEqual(manifest["execution"]["signalp_returncode"], 0)
        self.assertIs(manifest["correction"]["no_historical_run_overwritten"], True)
        self.assertEqual(manifest["execution"]["signalp_threads"], 32)
        # The server was never contacted from this task: unmeasured resources stay pending.
        self.assertEqual(manifest["execution"]["logical_cpu"], "pending")
        self.assertEqual(manifest["execution"]["available_memory"], "pending")

        amendment = self.run_dir / "results" / "provenance_amendment.json"
        readme = self.run_dir / "results" / "provenance_readme.md"
        self.assertTrue(amendment.is_file())
        self.assertTrue(readme.is_file())

        contract = json.loads((self.run_dir / "input_contract.json").read_text(encoding="utf-8"))
        self.assertEqual(contract["run_id"], RUN_ID)
        self.assertEqual(contract["execution"]["mode"], "server_dated_deploy")

    def test_records_authoritative_final_output_hashes(self) -> None:
        self.provenance()
        trail = json.loads(
            (self.run_dir / "results" / "provenance_outputs.json").read_text(encoding="utf-8")
        )
        amendment = self.run_dir / "results" / "provenance_amendment.json"
        readme = self.run_dir / "results" / "provenance_readme.md"

        self.assertEqual(trail["outputs"]["provenance_amendment.json"]["sha256"], sha256(amendment))
        self.assertEqual(trail["outputs"]["provenance_readme.md"]["sha256"], sha256(readme))
        self.assertEqual(trail["outputs"]["provenance_amendment.json"]["size"], amendment.stat().st_size)
        # The readme is written once, so its embedded self-hash must match the file.
        self.assertIs(trail["outputs"]["provenance_readme.md"]["embedded_self_hash"]["matches_final_file"], True)
        # The amendment json is rewritten after it hashes itself; the trail records the
        # truth about the embedded value rather than silently trusting it.
        self.assertIsInstance(
            trail["outputs"]["provenance_amendment.json"]["embedded_self_hash"]["matches_final_file"], bool
        )

    def test_binds_every_provenance_input_with_size_and_sha256(self) -> None:
        manifest = self.provenance()
        for key, path in (
            ("signalp_run_summary", self.signalp),
            ("foldseek_merge_manifest", self.foldseek),
            ("final_foldseek_matrix", self.matrix),
        ):
            record = manifest["inputs"][key]
            self.assertEqual(record["size"], path.stat().st_size)
            self.assertEqual(record["sha256"], sha256(path))

    def test_rejects_unsuccessful_signalp_summary(self) -> None:
        self.signalp.write_text(json.dumps({"status": "completed", "returncode": 1}), encoding="utf-8")
        with self.assertRaises(ValueError):
            self.provenance()

    def test_requires_candidate_only_completed_foldseek_manifest(self) -> None:
        self.foldseek.write_text(json.dumps({"status": "running"}), encoding="utf-8")
        with self.assertRaises(ValueError):
            self.provenance()


class HistoricalContractSnapshotTests(HousekeepingFixture):
    def test_snapshot_is_read_only_and_records_authorization(self) -> None:
        before_hash = sha256(self.historical_contract_json)
        before_mtime = self.historical_contract_json.stat().st_mtime_ns

        snapshot = hk.write_historical_contract_snapshot(
            run_dir=self.run_dir, run_id=RUN_ID, contract_paths=[self.historical_contract_json]
        )

        record = snapshot["contracts"]["20260916_phaded_full_library_signalp_01"]
        self.assertEqual(record["sha256"], before_hash)
        self.assertEqual(record["size"], self.historical_contract_json.stat().st_size)
        self.assertIs(record["authorization"]["server_execution_started"], False)
        self.assertIs(record["read_only_verified"], True)
        self.assertTrue(snapshot["access"] == "read_only")

        self.assertEqual(sha256(self.historical_contract_json), before_hash)
        self.assertEqual(self.historical_contract_json.stat().st_mtime_ns, before_mtime)
        self.assertEqual(
            (self.run_dir / "results" / "historical_contract_authorization_snapshot.json").is_file(), True
        )
        self.assertEqual((self.run_dir / "inputs" / "input_manifest.json").is_file(), True)


class CountDiscrepancyTests(HousekeepingFixture):
    def _fixture(self, rows: list[dict[str, str]], profile_counts: dict[str, int], subtype_counts: dict[str, int]) -> tuple[Path, Path]:
        matrix = self.root / "count_matrix.tsv"
        write_matrix(matrix, rows)
        manifest = self.root / "reconciliation_manifest.json"
        manifest.write_text(
            json.dumps(
                {
                    "run_id": "20260916_phaded_full_library_reconciliation_amendment_02",
                    "summary": {
                        "profile_evidence_status_counts": profile_counts,
                        "subtype_call_counts": subtype_counts,
                    },
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        return manifest, matrix

    def test_explains_delta_by_hold_branch_precedence(self) -> None:
        rows = [
            matrix_row("a1", "profile_ambiguous_family", "ambiguous_family_within_superfamily"),
            matrix_row("a2", "profile_ambiguous_family", "hold_architecture_conflict", domain="domain_conflict", confidence="hold"),
            matrix_row("a3", "profile_ambiguous_superfamily", "hold_architecture_conflict", domain="domain_conflict", confidence="hold", assignment="ambiguous_superfamily"),
            matrix_row("a4", "profile_ambiguous_superfamily", "ambiguous_superfamily", assignment="ambiguous_superfamily"),
        ]
        manifest, matrix = self._fixture(
            rows,
            {"profile_ambiguous_family": 2, "profile_ambiguous_superfamily": 2},
            {"ambiguous_family_within_superfamily": 1, "ambiguous_superfamily": 1, "hold_architecture_conflict": 2},
        )

        result = hk.analyze_count_discrepancy(manifest_path=manifest, matrix_path=matrix)

        self.assertEqual(result["status"], "resolved_layer_precedence")
        by_layer = {item["profile_evidence_status"]: item for item in result["resolutions"]}
        family = by_layer["profile_ambiguous_family"]
        self.assertEqual((family["profile_layer_count"], family["subtype_layer_count"], family["delta"]), (2, 1, 1))
        self.assertEqual(family["reclassified_to"], {"hold_architecture_conflict": 1})
        self.assertIs(family["delta_fully_accounted"], True)
        superfamily = by_layer["profile_ambiguous_superfamily"]
        self.assertEqual((superfamily["profile_layer_count"], superfamily["subtype_layer_count"], superfamily["delta"]), (2, 1, 1))
        self.assertIs(superfamily["delta_fully_accounted"], True)
        self.assertEqual({r["accession"] for r in result["reclassified_rows"]}, {"a2", "a3"})
        self.assertEqual(result["unexplained_reclassified_count"], 0)
        self.assertEqual(result["manifest_matrix_count_mismatches"], {})

    def test_flags_unexplained_reclassification(self) -> None:
        rows = [
            matrix_row("b1", "profile_ambiguous_family", "dPHASCL1_like_candidate"),
        ]
        manifest, matrix = self._fixture(
            rows,
            {"profile_ambiguous_family": 1},
            {"ambiguous_family_within_superfamily": 0},
        )

        result = hk.analyze_count_discrepancy(manifest_path=manifest, matrix_path=matrix)

        self.assertEqual(result["status"], "unresolved")
        self.assertEqual(result["unexplained_reclassified_count"], 1)
        self.assertEqual(result["reclassified_rows"][0]["resolution"], "unexplained")

    def test_flags_manifest_matrix_mismatch(self) -> None:
        rows = [matrix_row("c1", "profile_ambiguous_family", "ambiguous_family_within_superfamily")]
        manifest, matrix = self._fixture(
            rows,
            {"profile_ambiguous_family": 5},
            {"ambiguous_family_within_superfamily": 1},
        )

        result = hk.analyze_count_discrepancy(manifest_path=manifest, matrix_path=matrix)

        self.assertEqual(result["status"], "unresolved")
        self.assertIn("profile_evidence_status_counts.profile_ambiguous_family", result["manifest_matrix_count_mismatches"])

    def test_writes_amendment_and_leaves_historical_manifest_untouched(self) -> None:
        rows = [
            matrix_row("d1", "profile_ambiguous_family", "ambiguous_family_within_superfamily"),
            matrix_row("d2", "profile_ambiguous_family", "hold_architecture_conflict", domain="domain_conflict", confidence="hold"),
        ]
        manifest, matrix = self._fixture(
            rows,
            {"profile_ambiguous_family": 2},
            {"ambiguous_family_within_superfamily": 1, "hold_architecture_conflict": 1},
        )
        before_hash = sha256(manifest)
        before_mtime = manifest.stat().st_mtime_ns

        result = hk.write_count_discrepancy_amendment(run_dir=self.run_dir, run_id=RUN_ID, manifest_path=manifest, matrix_path=matrix)

        self.assertEqual(sha256(manifest), before_hash)
        self.assertEqual(manifest.stat().st_mtime_ns, before_mtime)
        self.assertIs(result["historical_manifest_unmodified"], True)
        for name in ("count_discrepancy_amendment.json", "count_discrepancy_report.md", "reclassified_ambiguous_hold_rows.tsv"):
            self.assertTrue((self.run_dir / "results" / name).is_file(), name)
        written = json.loads((self.run_dir / "results" / "count_discrepancy_amendment.json").read_text(encoding="utf-8"))
        self.assertEqual(written["status"], "resolved_layer_precedence")
        report = self.run_dir / "results" / "count_discrepancy_report.md"
        self.assertEqual(written["outputs"]["count_discrepancy_report.md"]["sha256"], sha256(report))
        trail = json.loads(
            (self.run_dir / "results" / "count_discrepancy_outputs.json").read_text(encoding="utf-8")
        )
        self.assertEqual(
            trail["outputs"]["count_discrepancy_amendment.json"]["sha256"],
            sha256(self.run_dir / "results" / "count_discrepancy_amendment.json"),
        )
        self.assertEqual(
            trail["outputs"]["reclassified_ambiguous_hold_rows.tsv"]["sha256"],
            sha256(self.run_dir / "results" / "reclassified_ambiguous_hold_rows.tsv"),
        )


if __name__ == "__main__":
    unittest.main()
