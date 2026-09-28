import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

from pipeline.scripts.amend_phaded_provenance import build_provenance_amendment


REPO = Path(__file__).resolve().parents[2]
RUN_CONTEXT_SCRIPT = REPO / "pipeline" / "scripts" / "run_context.py"
#: the core keys ``run_context.build_input_contract`` always emits
REQUIRED_CONTRACT_KEYS = (
    "schema_version",
    "run_id",
    "run_dir",
    "generated_at",
    "status",
    "gtdb",
    "inputs",
)
#: the top-level status domain of the shared contract implementation
CONTRACT_STATUS_DOMAIN = ("verified", "pending")


def load_run_context():
    """Load the shared contract implementation by path, as the script does."""
    spec = importlib.util.spec_from_file_location("run_context_conformance", RUN_CONTEXT_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class ProvenanceAmendmentTests(unittest.TestCase):
    def test_records_server_execution_without_relabeling_phenotype(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            signalp = root / "signalp.json"
            foldseek = root / "foldseek.json"
            matrix = root / "matrix.tsv"
            signalp.write_text(json.dumps({"status": "completed", "returncode": 0}), encoding="utf-8")
            foldseek.write_text(json.dumps({"status": "completed_candidate_only"}), encoding="utf-8")
            matrix.write_text("accession\nA\n", encoding="utf-8")
            result = build_provenance_amendment(
                run_id="20260916_test_provenance",
                output_dir=root / "out",
                signalp_summary=signalp,
                foldseek_manifest=foldseek,
                final_matrix=matrix,
                server_observation={
                    "server": "<SERVER_HOST>",
                    "account": "<SERVER_USER>",
                    "logical_cpu": 80,
                    "signalp_threads": 40,
                    "gpu_execution": "not_requested",
                },
            )
            self.assertEqual(result["execution"]["mode"], "server_dated_deploy")
            self.assertEqual(result["execution"]["server"], "<SERVER_HOST>")
            self.assertEqual(result["execution"]["signalp_returncode"], 0)
            self.assertTrue((root / "out" / "provenance_amendment.json").exists())
            self.assertIn("candidate-only", result["phenotype_boundary"])


class AmendmentInputContractConformanceTests(unittest.TestCase):
    """V2 closure: the amendment must write a ``run_context``-consistent contract.

    The Phase 1 P1 defect existed because this script wrote ``input_contract.json``
    from a bespoke dict (``schema_version`` / ``run_id`` / ``status`` / ``inputs`` /
    ``execution`` / ``phenotype_boundary``): no ``run_dir``, no ``generated_at``, no
    ``gtdb`` block, input entries without a ``status``, and a top-level ``status``
    holding the *completion* word ``completed_candidate_only`` instead of the
    contract domain ``{verified, pending}``.  The tests below are the conformance
    assertion that was missing, so the defect cannot come back on the next run.
    """

    RISK_LABEL = "completed_candidate_only"

    def setUp(self):
        self.run_context = load_run_context()
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.signalp = self.root / "signalp.json"
        self.foldseek = self.root / "foldseek.json"
        self.matrix = self.root / "matrix.tsv"
        self.signalp.write_text(
            json.dumps({"status": "completed", "returncode": 0}), encoding="utf-8"
        )
        self.foldseek.write_text(
            json.dumps({"status": "completed_candidate_only"}), encoding="utf-8"
        )
        self.matrix.write_text("accession\nA\n", encoding="utf-8")
        self.run_id = "20260916_test_provenance"
        self.output_dir = self.root / "out"
        self.manifest = build_provenance_amendment(
            run_id=self.run_id,
            output_dir=self.output_dir,
            signalp_summary=self.signalp,
            foldseek_manifest=self.foldseek,
            final_matrix=self.matrix,
            server_observation={
                "server": "<SERVER_HOST>",
                "account": "<SERVER_USER>",
                "logical_cpu": 80,
                "signalp_threads": 40,
                "gpu_execution": "not_requested",
            },
        )
        self.contract_path = self.root / "input_contract.json"
        self.contract = json.loads(self.contract_path.read_text(encoding="utf-8"))
        self.canonical = self.run_context.build_input_contract(
            self.root,
            run_id=self.run_id,
            gtdb_inputs={name: None for name in self.run_context.GTDB_INPUT_NAMES},
            inputs={
                "signalp_run_summary": self.signalp,
                "foldseek_merge_manifest": self.foldseek,
                "final_foldseek_matrix": self.matrix,
            },
        )

    def test_the_contract_is_written_next_to_the_amendment(self):
        self.assertTrue(self.contract_path.is_file())
        self.assertTrue((self.output_dir / "provenance_amendment.json").is_file())

    def test_contract_carries_every_run_context_core_key(self):
        for key in REQUIRED_CONTRACT_KEYS:
            self.assertIn(key, self.contract, key)
        self.assertEqual(self.contract["run_id"], self.run_id)
        self.assertEqual(self.contract["schema_version"], self.canonical["schema_version"])
        self.assertEqual(self.contract["run_dir"], str(self.root.resolve()))
        self.assertTrue(self.contract["generated_at"])

    def test_contract_matches_the_shared_implementation_for_the_same_inputs(self):
        """The written contract must be the shared one, not a bespoke look-alike."""
        self.assertEqual(self.contract["gtdb"], self.canonical["gtdb"])
        self.assertEqual(sorted(self.contract["inputs"]), sorted(self.canonical["inputs"]))
        for name, entry in self.contract["inputs"].items():
            self.assertEqual(sorted(entry), ["path", "sha256", "size", "status"], name)
            self.assertEqual(entry, self.canonical["inputs"][name], name)
        self.assertEqual(self.contract["status"], self.canonical["status"])

    def test_gtdb_slots_are_pending_without_invented_hashes(self):
        self.assertEqual(
            sorted(self.contract["gtdb"]), sorted(self.run_context.GTDB_INPUT_NAMES)
        )
        for name in self.run_context.GTDB_INPUT_NAMES:
            slot = self.contract["gtdb"][name]
            self.assertEqual(slot["status"], "pending", name)
            self.assertIsNone(slot["path"], name)
            self.assertIsNone(slot["sha256"], name)
            self.assertIsNone(slot["size"], name)

    def test_every_input_entry_declares_a_verified_status_and_hash(self):
        self.assertTrue(self.contract["inputs"])
        for name, entry in self.contract["inputs"].items():
            self.assertIn("status", entry, name)
            self.assertEqual(entry["status"], "verified", name)
            self.assertEqual(len(entry["sha256"]), 64, name)
            self.assertGreater(entry["size"], 0, name)

    def test_top_level_status_stays_inside_the_contract_domain(self):
        self.assertIn(self.contract["status"], CONTRACT_STATUS_DOMAIN)
        self.assertNotEqual(self.contract["status"], self.RISK_LABEL)

    def test_completion_state_is_carried_in_its_own_field(self):
        self.assertEqual(self.contract["run_completion_status"], self.RISK_LABEL)
        self.assertEqual(self.manifest["status"], self.RISK_LABEL)
        # the completion word must not leak into the contract's status fields
        self.assertNotEqual(self.contract["status"], self.RISK_LABEL)
        for name, entry in self.contract["inputs"].items():
            self.assertNotEqual(entry["status"], self.RISK_LABEL, name)

    def test_provenance_payload_is_preserved_beside_the_contract_core(self):
        self.assertEqual(self.contract["execution"], self.manifest["execution"])
        self.assertEqual(self.contract["phenotype_boundary"], self.manifest["phenotype_boundary"])
        self.assertEqual(
            self.contract["run_inputs_status"],
            "verified"
            if all(item["status"] == "verified" for item in self.contract["inputs"].values())
            else "pending",
        )


if __name__ == "__main__":
    unittest.main()
