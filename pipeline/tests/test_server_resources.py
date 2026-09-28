"""Governance tests for the dynamic server thread limit.

The project rule is: before any server launch, take
``min(40, nproc - one_minute_load - 10)`` and abort when no capacity remains.
The reserve exists so the host keeps roughly ten cores of headroom, and the
hard cap of 40 exists because the shared server must never be oversubscribed.

These tests are offline: they never touch the server, they only exercise the
pure calculation and the static wiring of the shell entrypoints.
"""
import os
import re
import subprocess
import sys
import unittest


TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(TESTS_DIR, "..", "scripts")
PYTHON = sys.executable


def read_script(name):
    with open(os.path.join(SCRIPTS, name), encoding="utf-8") as handle:
        return handle.read()


def load_module():
    sys.path.insert(0, os.path.abspath(SCRIPTS))
    try:
        import server_resources
    finally:
        sys.path.pop(0)
    return server_resources


module = load_module()


class ComputeThreadLimitTests(unittest.TestCase):
    def test_compute_thread_limit_reserves_ten_and_caps_at_forty(self):
        self.assertEqual(module.compute_thread_limit(80, 12.2), 40)
        self.assertEqual(module.compute_thread_limit(32, 10.1), 11)

    def test_compute_thread_limit_aborts_when_no_capacity_remains(self):
        with self.assertRaisesRegex(RuntimeError, "insufficient free capacity"):
            module.compute_thread_limit(16, 6.2)

    def test_compute_thread_limit_uses_floor_of_available_capacity(self):
        # 80 - 29.9 - 10 = 40.1 -> 40, still capped by the hard cap.
        self.assertEqual(module.compute_thread_limit(80, 29.9), 40)
        # 80 - 30.1 - 10 = 39.9 -> floor 39, below the hard cap.
        self.assertEqual(module.compute_thread_limit(80, 30.1), 39)

    def test_compute_thread_limit_accepts_explicit_reserve_and_cap(self):
        self.assertEqual(module.compute_thread_limit(40, 0.0, reserve=10, hard_cap=40), 30)
        self.assertEqual(module.compute_thread_limit(40, 0.0, reserve=5, hard_cap=12), 12)

    def test_compute_thread_limit_rejects_invalid_measurements(self):
        for args in ((0, 0.0), (-1, 0.0), (16, -0.1), (16, 0.0, -1), (16, 0.0, 10, 0)):
            with self.subTest(args=args):
                with self.assertRaises(ValueError):
                    module.compute_thread_limit(*args)

    def test_resource_record_carries_measurement_provenance(self):
        record = module.resource_record(
            total_cores=80, load_1m=12.2, reserve=10, hard_cap=40, threads=40,
            source="nproc+/proc/loadavg",
        )
        self.assertEqual(record["total_cores"], 80)
        self.assertEqual(record["load_1m"], 12.2)
        self.assertEqual(record["reserve"], 10)
        self.assertEqual(record["hard_cap"], 40)
        self.assertEqual(record["threads"], 40)
        self.assertEqual(record["basis"], "min(hard_cap, total_cores - load_1m - reserve)")
        self.assertEqual(record["source"], "nproc+/proc/loadavg")


class MeasurementParsingTests(unittest.TestCase):
    def test_parse_loadavg_reads_the_one_minute_value(self):
        self.assertAlmostEqual(
            module.parse_loadavg("12.20 9.10 8.30 3/1234 5678\n"), 12.2,
        )

    def test_parse_loadavg_rejects_malformed_input(self):
        for text in ("", "not-a-load", "1.0\n"):
            with self.subTest(text=text):
                with self.assertRaises(ValueError):
                    module.parse_loadavg(text)

    def test_read_measurements_fails_loudly_when_sources_are_unreadable(self):
        with self.assertRaises(RuntimeError) as ctx:
            module.read_measurements(
                cpuinfo_path="/nonexistent/cpuinfo", loadavg_path="/nonexistent/loadavg",
            )
        self.assertIn("cannot measure server capacity", str(ctx.exception))


class CliTests(unittest.TestCase):
    def test_cli_prints_an_auditable_json_record(self):
        result = subprocess.run(
            [PYTHON, os.path.join(SCRIPTS, "server_resources.py"),
             "--total-cores", "80", "--load-1m", "12.2"],
            capture_output=True, text=True, check=False,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        import json
        payload = json.loads(result.stdout)
        self.assertEqual(payload["threads"], 40)
        self.assertEqual(payload["total_cores"], 80)
        self.assertEqual(payload["load_1m"], 12.2)

    def test_cli_exits_non_zero_when_no_capacity_remains(self):
        result = subprocess.run(
            [PYTHON, os.path.join(SCRIPTS, "server_resources.py"),
             "--total-cores", "16", "--load-1m", "6.2"],
            capture_output=True, text=True, check=False,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("insufficient free capacity", result.stderr)

    def test_cli_rejects_a_requested_thread_count_above_the_measured_limit(self):
        result = subprocess.run(
            [PYTHON, os.path.join(SCRIPTS, "server_resources.py"),
             "--total-cores", "80", "--load-1m", "50.0", "--requested", "40"],
            capture_output=True, text=True, check=False,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("exceeds the measured limit", result.stderr)

    def test_cli_can_emit_a_resource_record_file(self):
        import json
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            out = os.path.join(tmp, "resources.json")
            result = subprocess.run(
                [PYTHON, os.path.join(SCRIPTS, "server_resources.py"),
                 "--total-cores", "80", "--load-1m", "0.0", "--out", out],
                capture_output=True, text=True, check=False,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            with open(out, encoding="utf-8") as handle:
                payload = json.load(handle)
            self.assertEqual(payload["threads"], 40)
            self.assertIn("measured_at_utc", payload)

    def test_cli_auto_detects_nproc_and_loadavg_when_not_supplied(self):
        script = read_script("server_resources.py")
        self.assertIn("nproc", script)
        self.assertIn("/proc/loadavg", script)


class EntrypointWiringTests(unittest.TestCase):
    """The 70-thread default is a governance defect: 80 - 70 = 10 means the
    driver assumed an idle host instead of measuring one."""

    def test_pipeline_driver_measures_capacity_instead_of_hardcoding_seventy(self):
        script = read_script("run_pipeline.sh")
        self.assertIn("server_resources.py", script)
        self.assertIn("measure_threads", script)
        self.assertIn("capacity:", script)
        for forbidden in ("THREADS_PREDICT=70", "THREADS_SCREEN=70", "THREADS=70"):
            self.assertNotIn(forbidden, script, forbidden)
        # The driver must not carry a stale launch value for either stage.
        self.assertNotRegex(script, r"THREADS_(?:PREDICT|SCREEN)=\d+")

    def test_stage_scripts_default_to_the_forty_thread_ceiling(self):
        for name in ("05_predict_proteins.sh", "06_screen.sh"):
            script = read_script(name)
            self.assertIn("THREADS=40", script, name)
            self.assertNotIn("THREADS=70", script, name)

    def test_stage_scripts_reject_values_above_the_measured_limit_when_server_mode_is_active(self):
        for name in ("05_predict_proteins.sh", "06_screen.sh"):
            script = read_script(name)
            self.assertIn("--server", script, name)
            self.assertIn("server_resources.py", script, name)
            self.assertIn("measured limit", script, name)

    def test_params_document_forty_as_a_ceiling_not_a_launch_value(self):
        with open(os.path.join(SCRIPTS, "..", "config", "params.yaml"), encoding="utf-8") as handle:
            params = handle.read()
        self.assertNotIn("threads: 70", params)
        # prediction + screening are declared as ceilings; phylogeny stays at 40.
        self.assertEqual(len(re.findall(r"(?m)^\s*threads:\s*40\b", params)), 3)
        self.assertIn("ceiling", params)
        self.assertIn("min(40, nproc - loadavg - 10)", params)


if __name__ == "__main__":
    unittest.main()
