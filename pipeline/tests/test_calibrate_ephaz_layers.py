import csv
import json
import os
import tempfile
import unittest
from unittest.mock import patch


from pipeline.scripts import calibrate_ephaz_layers


class LayeredCalibrationTests(unittest.TestCase):
    def test_parse_hits_uses_hmm_coordinate_columns_for_coverage(self):
        with tempfile.TemporaryDirectory() as root:
            tbl = os.path.join(root, "hits.tbl")
            dom = os.path.join(root, "hits.dom")
            with open(tbl, "w", encoding="utf-8") as handle:
                handle.write("TARGET - QUERY - 1e-20 100.0 0.0 1e-20 100.0 0.0\n")
            # HMMER domtblout: qlen=100 (field 5), HMM from/to=11..80
            fields = [
                "TARGET", "-", "120", "QUERY", "-", "100", "1e-20", "100.0", "0.0",
                "1", "1", "1", "1", "1", "1", "1e-20", "1e-20", "100.0", "0.0",
                "11", "80", "1", "70", "1", "70", "0.99", "description",
            ]
            with open(dom, "w", encoding="utf-8") as handle:
                handle.write(" ".join(fields) + "\n")
            hits = calibrate_ephaz_layers.parse_hits(tbl, dom)
            self.assertAlmostEqual(hits["TARGET"]["cov"], 0.7, places=6)

    def _write_controls(self, root):
        controls = os.path.join(root, "controls")
        os.makedirs(controls)
        with open(os.path.join(controls, "controls.tsv"), "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=["accession", "label", "query_group", "reviewed", "protein_name"],
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(
                [
                    {"accession": "CUR1", "label": "positive", "query_group": "e-PhaZ_EC", "reviewed": "true", "protein_name": "experimental core"},
                    {"accession": "BR1", "label": "positive", "query_group": "e-PhaZ_remote", "reviewed": "false", "protein_name": "putative broad"},
                    {"accession": "NEG1", "label": "negative", "query_group": "3HB_dehydrogenase", "reviewed": "true", "protein_name": "negative control"},
                ]
            )
        for name, rows in {
            "positive.faa": [("CUR1|e-PhaZ_EC", "M" * 100), ("BR1|e-PhaZ_remote", "M" * 100)],
            "negative.faa": [("NEG1|negative", "M" * 100)],
        }.items():
            with open(os.path.join(controls, name), "w", encoding="utf-8") as handle:
                for header, sequence in rows:
                    handle.write(f">{header}\n{sequence}\n")
        manifest = os.path.join(root, "layers.tsv")
        with open(manifest, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=["accession", "layer", "evidence", "architecture", "length"],
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(
                [
                    {"accession": "CUR1", "layer": "ePhaZ_curated_core", "evidence": "experimental;pmid:1", "architecture": "典型", "length": "100"},
                    {"accession": "BR1", "layer": "ePhaZ_broad_discovery", "evidence": "annotation_only", "architecture": "pending", "length": "100"},
                ]
            )
        return controls, manifest

    @staticmethod
    def _fake_hmmsearch(command, check, capture_output, text):
        # The probe FASTA is the final argument. Emit a deterministic hit for
        # only the layer-specific positive, leaving the broad set with one FN.
        tbl = command[command.index("--tblout") + 1]
        dom = command[command.index("--domtblout") + 1]
        hmm = os.path.basename(command[-2])
        target = "CUR1" if "curated" in hmm else "BR1"
        with open(tbl, "w", encoding="utf-8") as handle:
            handle.write(f"{target} - query - 1e-20 100.0 0.0 1e-20 100.0 0.0\n")
        fields = [target, "-", "100", "query", "-", "100", "1e-20", "100.0", "0.0", "1", "1", "1e-20", "1e-20", "100.0", "0.0", "1", "90", "1", "90", "1", "90", "0.99", "desc"]
        with open(dom, "w", encoding="utf-8") as handle:
            handle.write(" ".join(fields) + "\n")

    def test_calibrates_both_layers_and_records_evidence_metadata(self):
        with tempfile.TemporaryDirectory() as root:
            controls, manifest = self._write_controls(root)
            curated_hmm = os.path.join(root, "ePhaZ_curated_core.hmm")
            broad_hmm = os.path.join(root, "ePhaZ_broad_discovery.hmm")
            for path in (curated_hmm, broad_hmm):
                with open(path, "w", encoding="utf-8") as handle:
                    handle.write("HMMER3/f [test]\n")
            outdir = os.path.join(root, "out")
            with patch.object(calibrate_ephaz_layers.subprocess, "run", side_effect=self._fake_hmmsearch):
                calibrate_ephaz_layers.main(
                    [
                        "--curated-hmm", curated_hmm,
                        "--broad-hmm", broad_hmm,
                        "--controls", controls,
                        "--layer-manifest", manifest,
                        "--outdir", outdir,
                        "--hmmsearch-bin", "hmmsearch",
                        "--thresholds", "1e-2,1e-20",
                    ]
                )

            with open(os.path.join(outdir, "calibration_summary.tsv"), encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual({row["family"] for row in rows}, {"ePhaZ_curated_core", "ePhaZ_broad_discovery"})
            curated = [row for row in rows if row["family"] == "ePhaZ_curated_core" and row["threshold"] == "1e-02"][0]
            broad = [row for row in rows if row["family"] == "ePhaZ_broad_discovery" and row["threshold"] == "1e-02"][0]
            self.assertEqual((curated["TP"], curated["FN"], curated["FP"]), ("1", "0", "0"))
            self.assertEqual((broad["TP"], broad["FN"], broad["FP"]), ("1", "1", "0"))

            with open(os.path.join(outdir, "calibration_metadata.json"), encoding="utf-8") as handle:
                metadata = json.load(handle)
            self.assertEqual(set(metadata["layers"]), {"ePhaZ_curated_core", "ePhaZ_broad_discovery"})
            evidence = metadata["layers"]["ePhaZ_curated_core"]["positive_controls"][0]
            self.assertEqual(evidence["accession"], "CUR1")
            self.assertEqual(evidence["evidence"], "experimental;pmid:1")

    def test_requires_distinct_explicit_layer_hmms(self):
        with tempfile.TemporaryDirectory() as root:
            controls, manifest = self._write_controls(root)
            hmm = os.path.join(root, "ePhaZ.hmm")
            with open(hmm, "w", encoding="utf-8") as handle:
                handle.write("HMMER3/f [test]\n")
            with self.assertRaises(SystemExit):
                calibrate_ephaz_layers.main(
                    [
                        "--curated-hmm", hmm,
                        "--broad-hmm", hmm,
                        "--controls", controls,
                        "--layer-manifest", manifest,
                        "--outdir", os.path.join(root, "out"),
                    ]
                )

    def test_rejects_duplicate_control_accession_and_missing_core_evidence(self):
        with tempfile.TemporaryDirectory() as root:
            controls, manifest = self._write_controls(root)
            controls_tsv = os.path.join(controls, "controls.tsv")
            with open(controls_tsv, "a", encoding="utf-8") as handle:
                handle.write("CUR1\tpositive\te-PhaZ_EC\ttrue\tduplicate\n")
            for path in (os.path.join(root, "curated.hmm"), os.path.join(root, "broad.hmm")):
                with open(path, "w", encoding="utf-8") as handle:
                    handle.write("HMMER3/f [test]\n")
            with patch.object(calibrate_ephaz_layers.subprocess, "run", side_effect=self._fake_hmmsearch):
                with self.assertRaises(SystemExit):
                    calibrate_ephaz_layers.main(
                        [
                            "--curated-hmm", os.path.join(root, "curated.hmm"),
                            "--broad-hmm", os.path.join(root, "broad.hmm"),
                            "--controls", controls,
                            "--layer-manifest", manifest,
                            "--outdir", os.path.join(root, "out"),
                            "--hmmsearch-bin", "hmmsearch",
                        ]
                    )



if __name__ == "__main__":
    unittest.main()
