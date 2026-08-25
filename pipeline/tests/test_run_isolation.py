"""Static contract tests for dated, isolated pipeline runs.

These tests intentionally do not require Bash, GTDB, or a conda environment;
they protect the path contract that the HPC entrypoint and stage scripts share.
"""
import os
import re
import unittest


SCRIPTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts")


def read_script(name):
    with open(os.path.join(SCRIPTS, name), encoding="utf-8") as handle:
        return handle.read()


class RunIsolationTests(unittest.TestCase):
    def test_pipeline_defaults_to_a_new_dated_run_root(self):
        script = read_script("run_pipeline.sh")
        self.assertNotIn("REPO_REPO_ROOT", script)
        self.assertNotIn("`n", script)
        self.assertIn("REPO_ROOT=", script)
        self.assertIn("RUN_ROOT=", script)
        self.assertIn('PHB_RUN_ROOT="$RUN_ROOT"', script)
        self.assertIn('RUN_ROOT="$RUNS_ROOT/$RUN_ID"', script)
        self.assertIn("--legacy-root-results", script)
        self.assertIn("--run-dir", script)
        self.assertIn("--run-id", script)
        self.assertIn("git rev-parse --short HEAD", script)
        self.assertIn("maximum length is 64", script)
        self.assertIn("refusing to use symlinked runs directory", script)

    def test_pipeline_exposes_run_root_and_does_not_write_repo_results_by_default(self):
        script = read_script("run_pipeline.sh")
        self.assertIn('export PHB_RUN_ROOT', script)
        self.assertIn('mkdir -p "$RUN_ROOT/data" "$RUN_ROOT/inputs" "$RUN_ROOT/logs" "$RUN_ROOT/results"', script)
        self.assertIn('ln -s', script)
        self.assertIn('"$RUN_ROOT/results/run_manifest.json"', script)
        self.assertIn("write_input_contract", script)
        self.assertIn("hmm_source = hmm_root.resolve()", script)
        self.assertNotIn('MANIFEST_JSONL="$REPO_ROOT/results/', script)
        self.assertNotRegex(script, r"(?m)^\s+python(?:\s|$)")

    def test_stage_scripts_resolve_data_and_results_from_phb_run_root(self):
        for name in (
            "05_predict_proteins.sh",
            "06a_filter_shards.sh",
            "06_screen.sh",
            "08c_tier_rescore.sh",
            "09b_tier1_phylogeny.sh",
            "09g_fasttree.sh",
            "09_phylogeny.sh",
        ):
            script = read_script(name)
            self.assertIn("REPO_ROOT=", script, name)
            self.assertRegex(script, r'RUN_ROOT=.*PHB_RUN_ROOT', name)
            self.assertNotRegex(script, r'ROOT="\$SCRIPT_DIR/\.\./\.\."', name)
            self.assertNotIn("REPO_REPO_ROOT", script, name)
            self.assertNotIn("`n", script, name)

    def test_pipeline_links_read_only_hmm_input_into_run(self):
        script = read_script("run_pipeline.sh")
        self.assertIn('HMM_SOURCE="$REPO_ROOT/data/hmms"', script)
        self.assertIn('HMM_LINK="$RUN_ROOT/data/hmms"', script)
        self.assertIn('ln -s "$HMM_SOURCE" "$HMM_LINK"', script)


if __name__ == "__main__":
    unittest.main()
