"""Tests for the bound PhaDED candidate-profile hmmsearch runner (F3, test-first).

Task 5 of the evidence-model redesign split model state into two independent
fields: ``model_layer`` (what the model *is*) and ``functional_calibration_status``
(what has been *functionally calibrated*).  The legacy ``model_status`` column was
deliberately left unchanged for backward compatibility, so it still reads
``trained`` for HMMs that are only DISCOVERY-layer models (recall-only,
uncalibrated, measured at 78.69% sensitivity on the 563 x1 confounder hits — no
discriminating power).  Any decision keyed on ``model_status == "trained"``
therefore promotes a recall-only model into a discriminating family model, which
AGENTS.md forbids: the discovery layer only recalls and must never filter, drop
or demote a candidate, and must never produce a family call.

Load-bearing invariants pinned here:

* the decision column is ``model_layer``; ``model_status`` is never consulted;
* ``sequence_family_hmm_validated`` / ``calibrated_candidate_model`` are the only
  discriminating layers and the only ones scored into ``phaded_profile_scores.tsv``;
* a ``discovery_hmm_uncalibrated`` profile is never scored into the discriminating
  score table, and even when explicitly requested is written to a separate
  recall-only table that carries its layer;
* ``reference_query_only`` means there is no HMM at all: it is rejected as
  unscoreable instead of being scored (fail-closed), and it is reported;
* a manifest that predates ``model_layer`` (column absent) is treated as
  ``reference_query_only`` / non-discriminating and says so — "trained in the
  legacy column" is never read as "validated";
* ``model_layer`` and ``functional_calibration_status`` are independent: neither
  is inferred from the other.

Fixtures are synthetic and written to a temp directory; no project data is read
and no HMMER binary is executed (the external command runner is injected).
"""

from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import run_phaded_candidate_hmmsearch as module  # noqa: E402


REFERENCE_ONLY = "reference_query_only"
DISCOVERY = "discovery_hmm_uncalibrated"
VALIDATED = "sequence_family_hmm_validated"
CALIBRATED = "calibrated_candidate_model"

PROFILE_COLUMNS = [
    "profile_id", "profile_kind", "phaded_superfamily", "phaded_family_id",
    "model_layer", "model_status", "functional_calibration_status",
]

CANDIDATE_FASTA = (
    ">GCA_1|CTX_1_1 candidate\nMKWV\n"
    ">GCA_2|CTX_2_1 candidate\nMKWV\n"
)

#: One row per discriminator class.  ``legacy_trained`` is the trap: it carries
#: the legacy ``model_status=trained`` word while declaring the discovery layer.
PROFILES = [
    {"profile_id": "family_DED_hfam_2_6cf580a7126e", "profile_kind": "family",
     "phaded_superfamily": "intracellular nPHASCL with lipase box",
     "phaded_family_id": "DED_hfam_2", "model_layer": DISCOVERY,
     "model_status": "trained", "functional_calibration_status": "not_function_calibrated"},
    {"profile_id": "family_DED_hfam_70_9818be7f78e3", "profile_kind": "family",
     "phaded_superfamily": "extracellular native-SCL/PhaZ7-like",
     "phaded_family_id": "DED_hfam_70", "model_layer": VALIDATED,
     "model_status": "trained", "functional_calibration_status": "candidate_gate_passed_not_promoted"},
    {"profile_id": "superfamily_intracellular_nPHAMCL_b13be41db69c", "profile_kind": "superfamily",
     "phaded_superfamily": "intracellular nPHAMCL", "phaded_family_id": "",
     "model_layer": CALIBRATED, "model_status": "trained",
     "functional_calibration_status": "calibrated_candidate_model"},
    {"profile_id": "family_DED_hfam_3_43646f83e298", "profile_kind": "family",
     "phaded_superfamily": "periplasmic PHA depolymerases",
     "phaded_family_id": "DED_hfam_3", "model_layer": REFERENCE_ONLY,
     "model_status": "reference_only", "functional_calibration_status": "not_function_calibrated"},
]

TRAINED_ID = "family_DED_hfam_70_9818be7f78e3"
CALIBRATED_ID = "superfamily_intracellular_nPHAMCL_b13be41db69c"
DISCOVERY_ID = "family_DED_hfam_2_6cf580a7126e"
REFERENCE_ONLY_ID = "family_DED_hfam_3_43646f83e298"

TBL_TEMPLATE = (
    "# target name        accession  query name  accession  E-value  score  bias\n"
    "GCA_1|CTX_1_1        -          {profile}   -          {evalue}  120.0  0.1\n"
)


class FakeProcess:
    """Minimal ``subprocess.CompletedProcess`` lookalike."""

    def __init__(self, stdout: str = "", stderr: str = "") -> None:
        self.stdout = stdout
        self.stderr = stderr
        self.returncode = 0


class FakeHmmsearch:
    """Injected command runner: records argv and fakes tblout/version output."""

    def __init__(self, evalue: str = "1e-30") -> None:
        self.evalue = evalue
        self.calls: list[list[str]] = []

    def __call__(self, command, capture_output=True, text=True, check=True):
        argv = [str(item) for item in command]
        self.calls.append(argv)
        if "-h" in argv:
            return FakeProcess(stdout="\n# HMMER 3.4 (Aug 2023); http://hmmer.org/\n")
        tblout = Path(argv[argv.index("--tblout") + 1])
        profile = Path(argv[-2]).stem
        tblout.parent.mkdir(parents=True, exist_ok=True)
        tblout.write_text(TBL_TEMPLATE.format(profile=profile, evalue=self.evalue),
                          encoding="utf-8", newline="\n")
        return FakeProcess(stdout="[ok]\n")

    def scored_profiles(self) -> list[str]:
        out = []
        for argv in self.calls:
            if "-h" in argv:
                continue
            out.append(Path(argv[-2]).stem)
        return out


def write_manifest(path: Path, profiles, columns) -> Path:
    """Write a synthetic profile manifest; missing column values become empty."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t",
                                lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for profile in profiles:
            writer.writerow({column: profile.get(column, "") for column in columns})
    return path


class TempFixture(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="dsh_hmmsearch_f3_"))

    def make_deploy(self, profiles=None, columns=None, with_layer_column=True):
        deploy = self.tmp / "deploy"
        inputs = deploy / "inputs"
        inputs.mkdir(parents=True, exist_ok=True)
        (inputs / "candidate_union.faa").write_text(CANDIDATE_FASTA, encoding="utf-8",
                                                    newline="\n")
        if columns is None:
            columns = PROFILE_COLUMNS if with_layer_column else [
                column for column in PROFILE_COLUMNS if column != "model_layer"
            ]
        write_manifest(inputs / "profile_manifest.tsv",
                       PROFILES if profiles is None else profiles, columns)
        for profile in (PROFILES if profiles is None else profiles):
            if profile.get("model_layer") in {VALIDATED, CALIBRATED, DISCOVERY}:
                (inputs / f"{profile['profile_id']}.hmm").write_text("HMMER3/f\n",
                                                                     encoding="utf-8")
        return deploy


class LayerVocabularyTests(TempFixture):
    """The four layers and the discriminating/non-discriminating split."""

    def test_layer_vocabulary_is_exactly_the_four_governed_layers(self):
        self.assertEqual(
            {REFERENCE_ONLY, DISCOVERY, VALIDATED, CALIBRATED},
            module.MODEL_LAYERS,
        )

    def test_discriminating_and_non_discriminating_layers_partition_the_vocabulary(self):
        discriminating = set(module.DISCRIMINATING_MODEL_LAYERS)
        non_discriminating = set(module.NON_DISCRIMINATING_MODEL_LAYERS)
        self.assertEqual(module.MODEL_LAYERS, discriminating | non_discriminating)
        self.assertEqual(set(), discriminating & non_discriminating)
        # The discovery layer is recall-only by definition; the calibrated layer
        # exists only after an explicitly authorized promotion.
        self.assertIn(DISCOVERY, non_discriminating)
        self.assertNotIn(DISCOVERY, discriminating)
        self.assertIn(VALIDATED, discriminating)
        self.assertIn(CALIBRATED, discriminating)

    def test_unrecognized_declared_layer_fails_closed(self):
        with self.assertRaisesRegex(ValueError, "model_layer"):
            module.resolve_model_layer({"profile_id": "p1", "model_layer": "discovery_hmm"})

    def test_declared_layer_is_never_derived_from_model_status(self):
        # A row that still says "trained" in the legacy column, with no layer
        # column at all, stays non-discriminating.
        self.assertEqual(module.resolve_model_layer({"profile_id": "p1", "model_status": "trained"}),
                         REFERENCE_ONLY)
        # And a reference_only legacy word never demotes an explicit validated layer.
        self.assertEqual(
            module.resolve_model_layer({"profile_id": "p1", "model_layer": VALIDATED,
                                        "model_status": "reference_only"}),
            VALIDATED,
        )


class ManifestLayerTests(TempFixture):
    def test_legacy_manifest_without_layer_column_is_reference_query_only(self):
        manifest = self.make_deploy(with_layer_column=False) / "inputs" / "profile_manifest.tsv"
        rows, report = module.read_profile_manifest(manifest)
        self.assertFalse(report["model_layer_column_present"])
        self.assertEqual(report["layer_source"], "absent_model_layer_column")
        self.assertEqual({row["model_layer"] for row in rows}, {REFERENCE_ONLY})
        self.assertEqual(report["fallback_layer"], REFERENCE_ONLY)
        # The legacy word is still carried verbatim but is not the decision.
        self.assertEqual({row["model_status"] for row in rows}, {"trained", "reference_only"})
        self.assertEqual(report["model_layer_counts"], {REFERENCE_ONLY: 4})
        self.assertIn(REFERENCE_ONLY, report["fallback_note"])

    def test_declared_layer_column_is_read_verbatim(self):
        manifest = self.make_deploy() / "inputs" / "profile_manifest.tsv"
        rows, report = module.read_profile_manifest(manifest)
        by_id = {row["profile_id"]: row for row in rows}
        self.assertTrue(report["model_layer_column_present"])
        self.assertEqual(report["layer_source"], "model_layer_column")
        self.assertEqual(by_id[TRAINED_ID]["model_layer"], VALIDATED)
        self.assertEqual(by_id[DISCOVERY_ID]["model_layer"], DISCOVERY)
        self.assertEqual(report["model_layer_counts"][DISCOVERY], 1)

    def test_absent_layer_column_note_is_machine_readable(self):
        manifest = self.make_deploy(with_layer_column=False) / "inputs" / "profile_manifest.tsv"
        _, report = module.read_profile_manifest(manifest)
        # The fallback must be stated in the output, not silently assumed.
        self.assertIn("non-discriminating", report["fallback_note"])
        self.assertIn("never assumed validated", report["fallback_note"])


class SelectionTests(TempFixture):
    def selection(self, **kwargs):
        manifest = self.make_deploy() / "inputs" / "profile_manifest.tsv"
        return module.select_scoreable_profiles(module.read_profile_manifest(manifest)[0], **kwargs)

    def test_only_validated_layers_are_scored_by_default(self):
        selection = self.selection()
        self.assertEqual([row["profile_id"] for row in selection["validated"]],
                         [TRAINED_ID, CALIBRATED_ID])
        # discovery is not scored unless explicitly requested, and the
        # HMM-less reference profile is rejected rather than scored
        self.assertEqual(selection["discovery"], [])
        self.assertEqual([row["profile_id"] for row in selection["rejected"]],
                         [REFERENCE_ONLY_ID])

    def test_discovery_layer_is_never_discriminating_even_when_requested(self):
        selection = self.selection(include_discovery_layer=True)
        self.assertEqual([row["profile_id"] for row in selection["validated"]],
                         [TRAINED_ID, CALIBRATED_ID])
        self.assertEqual([row["profile_id"] for row in selection["discovery"]], [DISCOVERY_ID])
        # A discovery-layer profile can never appear among the discriminating ones.
        self.assertNotIn(DISCOVERY_ID, [row["profile_id"] for row in selection["validated"]])

    def test_reference_query_only_profile_is_rejected_as_unscoreable(self):
        selection = self.selection()
        self.assertEqual([row["profile_id"] for row in selection["rejected"]], [REFERENCE_ONLY_ID])
        self.assertEqual(selection["rejected"][0]["unscoreable_reason"],
                         module.UNSCOREABLE_REASON_REFERENCE_ONLY)
        self.assertEqual(selection["rejected"][0]["model_layer"], REFERENCE_ONLY)

    def test_legacy_manifest_trained_rows_are_rejected_not_scored(self):
        manifest = self.make_deploy(with_layer_column=False) / "inputs" / "profile_manifest.tsv"
        rows, _ = module.read_profile_manifest(manifest)
        selection = module.select_scoreable_profiles(rows, include_discovery_layer=True)
        self.assertEqual(selection["validated"], [])
        self.assertEqual(selection["discovery"], [])
        self.assertEqual(len(selection["rejected"]), len(rows))
        self.assertTrue(all(row["model_layer"] == REFERENCE_ONLY
                            for row in selection["rejected"]))


class RunTests(TempFixture):
    def run_script(self, *, include_discovery=False, with_layer_column=True, cpu=2):
        deploy = self.make_deploy(with_layer_column=with_layer_column)
        runner = FakeHmmsearch()
        run_dir = self.tmp / "run"
        result = module.run(deploy, run_dir, cpu=cpu, hmmsearch="hmmsearch",
                            runner=runner, include_discovery_layer=include_discovery)
        return run_dir, runner, result

    def test_validated_layers_are_the_only_profiles_executed(self):
        _, runner, result = self.run_script()
        self.assertEqual(runner.scored_profiles(), [TRAINED_ID, CALIBRATED_ID])
        self.assertEqual(result["validated_profiles"], 2)
        self.assertNotIn(DISCOVERY_ID, runner.scored_profiles())
        self.assertNotIn(REFERENCE_ONLY_ID, runner.scored_profiles())

    def test_discovery_layer_never_reaches_the_discriminating_score_table(self):
        run_dir, runner, result = self.run_script(include_discovery=True)
        self.assertIn(DISCOVERY_ID, runner.scored_profiles())
        self.assertEqual(result["discovery_profiles"], 1)
        with (run_dir / "results" / "phaded_profile_scores.tsv").open(encoding="utf-8") as handle:
            scored = {row["profile_id"] for row in csv.DictReader(handle, delimiter="\t")}
        self.assertEqual(scored, {TRAINED_ID, CALIBRATED_ID})
        self.assertNotIn(DISCOVERY_ID, scored)
        # no score row may exist for an unscoreable profile either
        self.assertNotIn(REFERENCE_ONLY_ID, scored)

    def test_discovery_scores_are_written_to_a_labelled_recall_only_table(self):
        run_dir, _, _ = self.run_script(include_discovery=True)
        recall_path = run_dir / "results" / "phaded_profile_scores_recall_only.tsv"
        self.assertTrue(recall_path.exists())
        with recall_path.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual({row["profile_id"] for row in rows}, {DISCOVERY_ID})
        self.assertEqual({row["model_layer"] for row in rows}, {DISCOVERY})
        self.assertEqual({row["recall_only"] for row in rows}, {"true"})
        self.assertEqual({row["family_call"] for row in rows}, {""})

    def test_provenance_reports_layers_without_reading_model_status(self):
        run_dir, _, _ = self.run_script(include_discovery=True)
        provenance = json.loads(
            (run_dir / "results" / "hmmsearch_provenance.json").read_text(encoding="utf-8")
        )
        self.assertEqual(provenance["validated_profile_count"], 2)
        self.assertEqual(provenance["discovery_profile_count"], 1)
        self.assertEqual(provenance["unscoreable_profiles"],
                         [{"profile_id": REFERENCE_ONLY_ID, "model_layer": REFERENCE_ONLY,
                           "reason": module.UNSCOREABLE_REASON_REFERENCE_ONLY}])
        self.assertEqual(provenance["model_layer_counts"][DISCOVERY], 1)
        self.assertTrue(provenance["model_layer_column_present"])
        self.assertEqual(provenance["layer_source"], "model_layer_column")
        self.assertEqual(provenance["fallback_note"], "")
        # The legacy column is not a decision input anywhere in the provenance.
        self.assertNotIn("trained_profile_count", provenance)

    def test_audit_table_lists_every_profile_with_its_layer(self):
        run_dir, _, _ = self.run_script(include_discovery=True)
        audit = run_dir / "results" / "model_layer_audit.tsv"
        self.assertTrue(audit.exists())
        with audit.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        by_id = {row["profile_id"]: row for row in rows}
        self.assertEqual(by_id[TRAINED_ID]["usage"], "validated_sequence_family")
        self.assertEqual(by_id[CALIBRATED_ID]["usage"], "validated_sequence_family")
        self.assertEqual(by_id[DISCOVERY_ID]["usage"], "recall_only")
        self.assertEqual(by_id[REFERENCE_ONLY_ID]["usage"], "rejected_unscoreable")
        self.assertEqual(by_id[DISCOVERY_ID]["model_status"], "trained")
        self.assertEqual(by_id[DISCOVERY_ID]["discriminating"], "false")

    def test_legacy_manifest_scores_nothing_and_says_why(self):
        run_dir, runner, result = self.run_script(with_layer_column=False)
        self.assertEqual(runner.scored_profiles(), [])
        self.assertEqual(result["validated_profiles"], 0)
        self.assertEqual(result["discovery_profiles"], 0)
        self.assertEqual(result["unscoreable_profiles"], 4)
        self.assertFalse(result["model_layer_column_present"])
        self.assertEqual(result["layer_source"], "absent_model_layer_column")
        provenance = json.loads(
            (run_dir / "results" / "hmmsearch_provenance.json").read_text(encoding="utf-8")
        )
        self.assertIn("non-discriminating", provenance["fallback_note"])
        self.assertIn("never assumed validated", provenance["fallback_note"])
        self.assertEqual(provenance["unscoreable_profiles"],
                         [{"profile_id": profile["profile_id"], "model_layer": REFERENCE_ONLY,
                           "reason": module.UNSCOREABLE_REASON_REFERENCE_ONLY}
                          for profile in PROFILES])

    def test_cpu_guard_is_unchanged(self):
        deploy = self.make_deploy()
        with self.assertRaisesRegex(ValueError, "cpu"):
            module.run(deploy, self.tmp / "run", cpu=41, hmmsearch="hmmsearch",
                       runner=FakeHmmsearch())

    def test_main_exits_cleanly_on_a_legacy_manifest(self):
        deploy = self.make_deploy(with_layer_column=False)
        code = module.main(["--deploy-dir", str(deploy), "--run-dir", str(self.tmp / "run"),
                            "--cpu", "1", "--hmmsearch", "hmmsearch"],
                           runner=FakeHmmsearch())
        self.assertEqual(code, 0)


class IndependenceTests(TempFixture):
    def test_functional_calibration_status_is_never_derived_from_the_layer(self):
        manifest = self.make_deploy() / "inputs" / "profile_manifest.tsv"
        rows, _ = module.read_profile_manifest(manifest)
        by_id = {row["profile_id"]: row for row in rows}
        # A validated sequence layer is NOT functionally calibrated...
        self.assertEqual(by_id[TRAINED_ID]["functional_calibration_status"],
                         "candidate_gate_passed_not_promoted")
        # ...while the input field is carried through verbatim, and an absent
        # value is reported as unknown rather than being inferred.
        self.assertEqual(by_id[CALIBRATED_ID]["functional_calibration_status"],
                         "calibrated_candidate_model")

    def test_absent_functional_calibration_status_is_not_inferred(self):
        columns = [column for column in PROFILE_COLUMNS
                   if column != "functional_calibration_status"]
        profiles = [dict(PROFILES[1])]
        deploy = self.make_deploy(profiles=profiles, columns=columns)
        rows, _ = module.read_profile_manifest(deploy / "inputs" / "profile_manifest.tsv")
        self.assertEqual(rows[0]["functional_calibration_status"], "unknown")


if __name__ == "__main__":
    unittest.main()
