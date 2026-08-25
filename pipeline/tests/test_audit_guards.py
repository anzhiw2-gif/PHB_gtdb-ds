"""Audit guard regression tests that require no GTDB or external tools."""
import importlib.util
import json
import os
import sys
import tempfile
import unittest
from unittest.mock import patch


SCRIPTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts")


def load_module(name):
    path = os.path.join(SCRIPTS, f"{name}.py")
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class AuditGuardTests(unittest.TestCase):
    def test_manifest_rejects_missing_required_output(self):
        manifest = load_module("run_manifest")
        with tempfile.TemporaryDirectory() as tmp:
            present = os.path.join(tmp, "present.tsv")
            with open(present, "w", encoding="utf-8") as handle:
                handle.write("ok\n")
            with self.assertRaises(manifest.ManifestError):
                manifest.validate_paths([present, os.path.join(tmp, "missing.tsv")], "output")

    def test_manifest_rejects_empty_or_failed_step_records(self):
        manifest = load_module("run_manifest")
        with self.assertRaises(manifest.ManifestError):
            manifest.validate_steps([])
        with self.assertRaises(manifest.ManifestError):
            manifest.validate_steps([{"step": "06_screen", "exit_code": 1}])

    def test_manifest_includes_successful_finalizer_step(self):
        manifest = load_module("run_manifest")
        built = manifest.build_manifest(
            [{"step": "stage", "exit_code": 0}], [], [],
            final_step={"step": "finalize_manifest", "exit_code": 0},
        )
        self.assertEqual(
            [step["step"] for step in built["steps"]],
            ["stage", "finalize_manifest"],
        )

    def test_manifest_strict_provenance_records_source_inputs_environment_and_commands(self):
        manifest = load_module("run_manifest")
        with tempfile.TemporaryDirectory() as tmp:
            source = os.path.join(tmp, "run_pipeline.sh")
            gtdb = os.path.join(tmp, "bac120.faa")
            hmm = os.path.join(tmp, "ePhaZ.hmm")
            for path, contents in ((source, "source\n"), (gtdb, "gtdb\n"), (hmm, "hmm\n")):
                with open(path, "w", encoding="utf-8") as handle:
                    handle.write(contents)
            built = manifest.build_manifest(
                [{"step": "stage", "exit_code": 0, "command": ["python", "stage.py"]}],
                [gtdb], [],
                source_files=[source],
                gtdb_inputs=[gtdb],
                hmm_inputs=[hmm],
                environment={"entrypoint": "test", "version": "1"},
                strict=True,
            )
            source_hash = manifest.sha256(source)
            gtdb_hash = manifest.sha256(gtdb)
            hmm_hash = manifest.sha256(hmm)
        self.assertEqual(built["source_files"][source], source_hash)
        self.assertTrue(built["source_bundle_sha256"])
        self.assertEqual(built["gtdb_inputs"][gtdb], gtdb_hash)
        self.assertEqual(built["hmm_inputs"][hmm], hmm_hash)
        self.assertEqual(built["environment"]["entrypoint"], "test")
        self.assertEqual(built["commands"], [["python", "stage.py"]])

    def test_manifest_strict_provenance_rejects_missing_declarations(self):
        manifest = load_module("run_manifest")
        with tempfile.TemporaryDirectory() as tmp:
            present = os.path.join(tmp, "input.tsv")
            with open(present, "w", encoding="utf-8") as handle:
                handle.write("ok\n")
            with self.assertRaises(manifest.ManifestError):
                manifest.build_manifest(
                    [{"step": "stage", "exit_code": 0, "command": ["stage"]}],
                    [present], [], strict=True,
                )

    def test_manifest_binds_run_context_and_pending_input_contract(self):
        manifest = load_module("run_manifest")
        with tempfile.TemporaryDirectory() as tmp:
            root = os.path.abspath(tmp)
            source = os.path.join(root, "source.py")
            hmm = os.path.join(root, "family.hmm")
            output = os.path.join(root, "output.tsv")
            contract = os.path.join(root, "input_contract.json")
            for path in (source, hmm, output):
                with open(path, "w", encoding="utf-8") as handle:
                    handle.write("ok\n")
            with open(contract, "w", encoding="utf-8") as handle:
                handle.write(json.dumps({
                    "schema_version": "1.0", "status": "pending",
                    "run_id": "run-1", "run_dir": root,
                }) + "\n")
            built = manifest.build_manifest(
                [{"step": "stage", "exit_code": 0, "command": ["stage"]}],
                [output], [output], source_files=[source],
                gtdb_inputs=[os.path.join(root, "missing.tsv")], hmm_inputs=[hmm],
                strict=True, run_id="run-1", run_root=root,
                input_contract=contract, allow_pending_gtdb=True,
            )
        self.assertEqual(built["run_id"], "run-1")
        self.assertEqual(built["run_root"], root)
        self.assertEqual(built["input_contract"]["status"], "pending")
        self.assertIsNone(next(iter(built["gtdb_inputs"].values())))

    def test_manifest_rejects_mismatched_input_contract_context(self):
        manifest = load_module("run_manifest")
        with tempfile.TemporaryDirectory() as tmp:
            root = os.path.abspath(tmp)
            source = os.path.join(root, "source.py")
            hmm = os.path.join(root, "family.hmm")
            output = os.path.join(root, "output.tsv")
            contract = os.path.join(root, "input_contract.json")
            for path in (source, hmm, output):
                with open(path, "w", encoding="utf-8") as handle:
                    handle.write("ok\n")
            with open(contract, "w", encoding="utf-8") as handle:
                json.dump({"schema_version": "1.0", "run_id": "other", "run_dir": root}, handle)
            with self.assertRaises(manifest.ManifestError):
                manifest.build_manifest(
                    [{"step": "stage", "exit_code": 0, "command": ["stage"]}],
                    [output], [output], source_files=[source], gtdb_inputs=[],
                    hmm_inputs=[hmm], strict=True, run_id="run-1", run_root=root,
                    input_contract=contract,
                )

    def test_prediction_manifest_rejects_incomplete_genome_count(self):
        prediction = load_module("05_validate_prediction_manifest")
        with self.assertRaises(prediction.PredictionManifestError):
            prediction.validate_manifest({
                "total_genomes": 10, "predicted_genomes": 9, "failed_genomes": 0,
                "expected_shards": 1, "shards": [],
            })

    def test_prediction_manifest_requires_the_declared_shard_set(self):
        prediction = load_module("05_validate_prediction_manifest")
        with self.assertRaises(prediction.PredictionManifestError):
            prediction.validate_manifest({
                "total_genomes": 2,
                "predicted_genomes": 2,
                "failed_genomes": 0,
                "expected_shards": ["shard_0001.faa", "shard_0002.faa"],
                "shards": [
                    {"name": "shard_0001.faa", "sha256": "a", "genomes": 2},
                ],
            })

    def test_prediction_script_does_not_mask_parallel_or_manifest_failures(self):
        with open(os.path.join(SCRIPTS, "05_predict_proteins.sh"), encoding="utf-8") as handle:
            script = handle.read()
        self.assertNotIn("parallel -j \"$THREADS\" --progress predict_one {} 2> \"$LOG/predict_progress.log\" || true", script)
        self.assertIn("prediction_manifest.json", script)
        self.assertIn("05_validate_prediction_manifest.py", script)

    def test_tree_manifest_marks_changed_tier_input_stale(self):
        tree_manifest = load_module("09i_tree_manifest")
        status = tree_manifest.classify_status(
            n_leaves=12,
            current_tier_count=12,
            recorded_input_sha256="old-sha",
            current_input_sha256="new-sha",
        )
        self.assertEqual(status, "stale_input")

    def test_tree_manifest_preserves_recorded_hash_when_input_changes(self):
        tree_manifest = load_module("09i_tree_manifest")
        hashes = tree_manifest.input_hashes("old-sha", "new-sha")
        self.assertEqual(hashes["recorded_input_sha256"], "old-sha")
        self.assertEqual(hashes["current_input_sha256"], "new-sha")

    def test_cluster_summary_records_hit_locus_and_genome_counts(self):
        clusters = load_module("11_clusters")
        with tempfile.TemporaryDirectory() as tmp:
            out = os.path.join(tmp, "cluster_summary.tsv")
            clusters.write_cluster_summary(
                out,
                marker_hits={("ArchPhaZ_patatin", "PhaC"): 3},
                supporting_loci={("ArchPhaZ_patatin", "PhaC"): {"locus_a", "locus_b"}},
                supporting_genomes={("ArchPhaZ_patatin", "PhaC"): {"GCA_1"}},
            )
            with open(out, encoding="utf-8") as handle:
                lines = handle.read().splitlines()
        self.assertEqual(
            lines[0],
            "hit_family\tmarker_family\tmarker_hits\tsupporting_loci\tsupporting_genomes",
        )
        self.assertEqual(lines[1], "ArchPhaZ_patatin\tPhaC\t3\t2\t1")

    def test_cluster_rejects_missing_declared_marker_hmm(self):
        clusters = load_module("11_clusters")
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(clusters.ClusterInputError):
                clusters.resolve_marker_hmms(tmp, ["PhaC", "PhaE"])

    def test_cluster_skips_hmmsearch_for_empty_faa_and_records_audit(self):
        clusters = load_module("11_clusters")
        with tempfile.TemporaryDirectory() as tmp:
            faa = os.path.join(tmp, "empty.faa")
            hmm = os.path.join(tmp, "PhaC.hmm")
            open(faa, "w", encoding="utf-8").close()
            with open(hmm, "w", encoding="utf-8") as handle:
                handle.write("placeholder\n")
            clusters.HMMSEARCH = sys.executable
            with patch.object(clusters, "run", side_effect=AssertionError("hmmsearch should be skipped")):
                self.assertEqual(clusters.annotate_markers(faa, [hmm], tmp, 1), {})
            with open(os.path.join(tmp, "invalid_fasta_records.tsv"), encoding="utf-8") as handle:
                self.assertIn("empty.faa\t1\tinvalid_fasta_record", handle.read())

    def test_parallel_cluster_resume_clears_stale_rc_before_launch(self):
        with open(os.path.join(SCRIPTS, "11_clusters_parallel.sh"), encoding="utf-8") as handle:
            script = handle.read()
        self.assertIn('rm -f "$root/rc"', script)

    def test_cluster_rejects_incomplete_locus_audit(self):
        clusters = load_module("11_clusters")
        with self.assertRaises(clusters.ClusterInputError):
            clusters.require_complete_locus_audit([
                {"genome": "GCA_1", "locus": "gene_1", "family": "ePhaZ", "status": "missing_genome"},
            ])

    def test_cluster_audit_writer_preserves_non_analyzed_reason(self):
        clusters = load_module("11_clusters")
        with tempfile.TemporaryDirectory() as tmp:
            out = os.path.join(tmp, "cluster_audit.tsv")
            clusters.write_cluster_audit(out, [
                {"genome": "GCA_1", "locus": "gene_1", "family": "ePhaZ", "status": "missing_genome"},
            ])
            with open(out, encoding="utf-8") as handle:
                lines = handle.read().splitlines()
        self.assertEqual(lines[0], "genome\tlocus\tfamily\tstatus")
        self.assertEqual(lines[1], "GCA_1\tgene_1\tePhaZ\tmissing_genome")

    def test_cluster_genome_audit_summarizes_incomplete_reasons(self):
        clusters = load_module("11_clusters")
        with tempfile.TemporaryDirectory() as tmp:
            out = os.path.join(tmp, "cluster_genome_audit.tsv")
            clusters.write_cluster_genome_audit(out, [
                {"genome": "GCA_1", "locus": "gene_1", "family": "ePhaZ", "status": "analyzed"},
                {"genome": "GCA_1", "locus": "gene_2", "family": "ePhaZ", "status": "locus_not_found"},
            ])
            with open(out, encoding="utf-8") as handle:
                lines = handle.read().splitlines()
        self.assertEqual(lines[0], "genome\trequested_loci\tanalyzed_loci\tnot_analyzed_statuses")
        self.assertEqual(lines[1], "GCA_1\t2\t1\tlocus_not_found:1")

    def test_aggregate_validates_tbl_dom_pairs_but_allows_empty_tbl(self):
        aggregate = load_module("06b_aggregate_hits")
        with tempfile.TemporaryDirectory() as tmp:
            tbl = os.path.join(tmp, "ePhaZ__shard_0001.tbl")
            dom = os.path.join(tmp, "ePhaZ__shard_0001.dom")
            with open(tbl, "w", encoding="utf-8"):
                pass
            with open(dom, "w", encoding="utf-8"):
                pass
            self.assertEqual(aggregate.validate_pairs([tbl]), [])

            os.remove(dom)
            with self.assertRaises(aggregate.AggregateError):
                aggregate.validate_pairs([tbl])

    def test_screen_manifest_requires_complete_family_shard_matrix(self):
        screen_manifest = load_module("06_validate_screen_manifest")
        with self.assertRaises(screen_manifest.ScreenManifestError):
            screen_manifest.validate_manifest({
                "families": ["ePhaZ"],
                "shards": [{"name": "shard_0001.faa", "sha256": "input"}],
                "tasks": [],
            })

    def test_screen_manifest_supports_revalidating_existing_manifest(self):
        screen_manifest = load_module("06_validate_screen_manifest")
        manifest = {
            "families": ["ePhaZ"],
            "shards": [{"name": "shard_0001.faa", "sha256": "input"}],
            "tasks": [{
                "family": "ePhaZ", "shard": "shard_0001.faa",
                "input_sha256": "input", "hmm_sha256": "hmm",
                "tbl_sha256": "tbl", "dom_sha256": "dom", "evalue": "1e-5",
            }],
        }
        self.assertIs(screen_manifest.validate_manifest(manifest), manifest)

    def test_filter_manifest_requires_each_declared_input_shard(self):
        filter_manifest = load_module("06a_validate_filter_manifest")
        with self.assertRaises(filter_manifest.FilterManifestError):
            filter_manifest.validate_manifest({
                "source_shards": [{"name": "shard_0001.faa", "sha256": "input"}],
                "filtered_shards": [],
                "max_aa": 100000,
            })

    def test_filter_and_screen_scripts_publish_provenance_manifests(self):
        for script_name, manifest_name in (
            ("06a_filter_shards.sh", "06a_validate_filter_manifest.py"),
            ("06_screen.sh", "06_validate_screen_manifest.py"),
        ):
            with open(os.path.join(SCRIPTS, script_name), encoding="utf-8") as handle:
                script = handle.read()
            self.assertIn(manifest_name, script)
            self.assertIn(".build.$$.tmp", script)

    def test_screen_fails_closed_when_declared_hmm_is_missing(self):
        with open(os.path.join(SCRIPTS, "06_screen.sh"), encoding="utf-8") as handle:
            script = handle.read()
        self.assertIn('echo "[ERROR] $fam HMM 缺失"', script)
        self.assertIn('exit 1', script[script.index('if [ ! -s "$hmm" ]'):])

    def test_execution_chain_has_no_masked_failures_or_warning_continue(self):
        for script_name in (
            "05_predict_proteins.sh",
            "06_screen.sh",
            "08c_tier_rescore.sh",
            "run_pipeline.sh",
        ):
            with open(os.path.join(SCRIPTS, script_name), encoding="utf-8") as handle:
                script = handle.read()
            self.assertNotIn("|| true", script, script_name)
            self.assertNotRegex(script, r"\bcontinue\b", script_name)
            self.assertIn("set -Eeuo pipefail", script, script_name)

    def test_screen_aborts_on_parallel_failure(self):
        with open(os.path.join(SCRIPTS, "06_screen.sh"), encoding="utf-8") as handle:
            script = handle.read()
        self.assertIn('exit "$local_rc"', script)

    def test_tier_rescore_rejects_empty_declared_inputs_and_outputs(self):
        with open(os.path.join(SCRIPTS, "08c_tier_rescore.sh"), encoding="utf-8") as handle:
            script = handle.read()
        self.assertIn('if [ ! -s "$faa" ]', script)
        self.assertIn('if [ ! -s "$ROOT/$hmm" ]', script)
        self.assertIn('require_nonempty "$fam tier2 FASTA" "$tier2_faa"', script)
        self.assertIn('require_nonempty "$fam tier1 FASTA" "$tier1_faa"', script)

    def test_tier_rescore_does_not_skip_declared_family_inputs(self):
        with open(os.path.join(SCRIPTS, "08c_tier_rescore.sh"), encoding="utf-8") as handle:
            script = handle.read()
        self.assertIn(".build.$$.tmp", script)
        self.assertNotIn("[ -f \"$faa\" ] || continue", script)
        self.assertIn('fail "$fam', script)

    def test_pipeline_requires_explicit_opt_in_for_paused_phylogeny(self):
        with open(os.path.join(SCRIPTS, "run_pipeline.sh"), encoding="utf-8") as handle:
            script = handle.read()
        self.assertIn("RUN_PHYLOGENY", script)
        self.assertIn("--run-phylogeny", script)
        self.assertIn('if [ "$RUN_PHYLOGENY" -eq 1 ]; then', script)
        self.assertIn("skip paused phylogeny", script)
        self.assertIn('if [ "$RUN_PHYLOGENY" -eq 1 ]; then', script[script.index("09i_tree_manifest"):])
        self.assertIn("FINAL_OUTPUTS", script)

    def test_pipeline_captures_step_exit_code_before_failing_closed_exit(self):
        with open(os.path.join(SCRIPTS, "run_pipeline.sh"), encoding="utf-8") as handle:
            script = handle.read()
        self.assertIn('set +e\n    "$@" 2>&1 | tee -a "$LOG"', script)
        self.assertIn('local -a pipe_status=("${PIPESTATUS[@]}")', script)
        self.assertIn('step_end "$rc"', script)

    def test_pipeline_finalizer_binds_strict_provenance_inputs_and_source_bundle(self):
        with open(os.path.join(SCRIPTS, "run_pipeline.sh"), encoding="utf-8") as handle:
            script = handle.read()
        for required in (
            "--source-files",
            "--gtdb-inputs",
            "--hmm-inputs",
            "--strict-provenance",
            "--final-step-command",
            '"command": command',
        ):
            self.assertIn(required, script)

    def test_pipeline_provenance_lists_both_hmm_versions_used_by_stages(self):
        with open(os.path.join(SCRIPTS, "run_pipeline.sh"), encoding="utf-8") as handle:
            script = handle.read()
        for path in (
            "data/hmms/ePhaZ.hmm",
            "data/hmms/iPhaZ.hmm",
            "data/hmms/OH.hmm",
            "data/hmms/v2/ePhaZ.hmm",
            "data/hmms/v2/iPhaZ.hmm",
            "data/hmms/v2/OH.hmm",
            "data/hmms/v2/BdhA.hmm",
            "data/hmms/v2/ArchPhaZ_patatin.hmm",
            "data/hmms/v2/ArchPhaZ_hydrolase.hmm",
            "data/hmms/v2/PhaJ.hmm",
            "data/hmms/v2/phasin.hmm",
            "data/hmms/v2/PhaC.hmm",
        ):
            self.assertIn(path, script)
        self.assertIn('"$ROOT/environment.yml"', script)
        self.assertIn('"$ROOT/pipeline/config/params.yaml"', script)


    def test_pipeline_validates_prediction_manifest_and_runs_required_handoffs(self):
        with open(os.path.join(SCRIPTS, "run_pipeline.sh"), encoding="utf-8") as handle:
            script = handle.read()
        for required in (
            "05_validate_prediction_manifest.py",
            "07b_extract_seqs.py",
            "09a_tier1_summary.py",
            "09i_tree_manifest.py",
            "cluster_locus_audit.tsv",
            "cluster_genome_audit.tsv",
        ):
            self.assertIn(required, script)
        self.assertNotIn("stable=0", script)
        self.assertLess(script.index("05_validate_prediction_manifest.py"), script.index("06a_filter_shards.sh"))
        self.assertLess(script.index("07b_extract_seqs.py"), script.index("08_validate.py"))
        self.assertLess(script.index("09a_tier1_summary.py"), script.index("09b_tier1_phylogeny.sh"))
        self.assertLess(script.index("09i_tree_manifest.py"), script.index("10_distribution.py"))


if __name__ == "__main__":
    unittest.main()
