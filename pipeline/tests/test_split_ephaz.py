"""Failing contract tests for the two-layer ePhaZ seed splitter.

These tests deliberately use a small synthetic library.  The production
splitter must keep every input accession, make the two FASTA outputs
disjoint, and fail closed when provenance or sequence inputs are incomplete.
"""

import csv
import hashlib
import json
import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "split_ephaz_seeds.py"


class SplitEphazContractTests(unittest.TestCase):
    def _load_splitter(self):
        self.assertTrue(SCRIPT.is_file(), f"missing splitter implementation: {SCRIPT}")
        spec = importlib.util.spec_from_file_location("split_ephaz_seeds", SCRIPT)
        self.assertIsNotNone(spec)
        module = importlib.util.module_from_spec(spec)
        assert spec.loader is not None
        spec.loader.exec_module(module)
        self.assertTrue(
            hasattr(module, "split_ephaz_seeds"),
            "split_ephaz_seeds.py must expose split_ephaz_seeds()",
        )
        return module

    @staticmethod
    def _write_fixture(root):
        seed_fasta = root / "ephaz.faa"
        manifest = root / "manifest.tsv"
        evidence = root / "curated_evidence.tsv"
        outdir = root / "split"

        records = []
        manifest_rows = []
        evidence_rows = []

        # Four long, complete sequences with explicit experimental support.
        curated = {
            "C001": (340, "Bacillus alpha"),
            "C002": (420, "Pseudomonas beta"),
            "C003": (278, "Delftia gamma"),
            "C004": (494, "Talaromyces delta"),
        }
        for accession, (length, organism) in curated.items():
            records.append((accession, length, organism))
            manifest_rows.append(
                {
                    "accession": accession,
                    "family": "ePhaZ",
                    "organism": organism,
                    "protein_name": "PHB depolymerase",
                    "reviewed": "true",
                    "evidence": "experimental",
                    "architecture": "typical",
                    "completeness": "complete",
                }
            )
            evidence_rows.append(
                {
                    "accession": accession,
                    "evidence_level": "experimental",
                    "architecture_status": "typical",
                    "completeness_status": "complete",
                    "reference": "pmid:10000000",
                }
            )

        # A long sequence with annotation-only support and unresolved
        # architecture must remain discoverable, but cannot enter core.
        records.append(("B001", 281, "Remote organism"))
        manifest_rows.append(
            {
                "accession": "B001",
                "family": "ePhaZ",
                "organism": "Remote organism",
                "protein_name": "putative depolymerase",
                "reviewed": "false",
                "evidence": "annotation_only",
                "architecture": "pending",
                "completeness": "complete",
            }
        )
        evidence_rows.append(
            {
                "accession": "B001",
                "evidence_level": "annotation_only",
                "architecture_status": "pending",
                "completeness_status": "complete",
                "reference": "",
            }
        )

        # Exactly 61 short records.  S001 intentionally has strong evidence;
        # the length gate must still prevent it from entering curated_core.
        for index in range(1, 62):
            accession = f"S{index:03d}"
            length = 100 + ((index * 7) % 100)
            organism = f"Short organism {index}"
            records.append((accession, length, organism))
            manifest_rows.append(
                {
                    "accession": accession,
                    "family": "ePhaZ",
                    "organism": organism,
                    "protein_name": "short putative depolymerase",
                    "reviewed": "false",
                    "evidence": "annotation_only",
                    "architecture": "pending",
                    "completeness": "unknown",
                }
            )
            if accession == "S001":
                evidence_rows.append(
                    {
                        "accession": accession,
                        "evidence_level": "experimental",
                        "architecture_status": "typical",
                        "completeness_status": "complete",
                        "reference": "pmid:10000001",
                    }
                )

        with seed_fasta.open("w", encoding="utf-8", newline="\n") as handle:
            for accession, length, organism in records:
                # A deterministic valid protein alphabet is sufficient for a
                # partition contract test; no biological claim is implied.
                sequence = ("ACDEFGHIKLMNPQRSTVWY" * ((length // 20) + 1))[:length]
                handle.write(f">{accession}|e-PhaZ|{organism}\n{sequence}\n")

        columns = list(manifest_rows[0])
        with manifest.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(manifest_rows)

        evidence_columns = list(evidence_rows[0])
        with evidence.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=evidence_columns, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(evidence_rows)

        return seed_fasta, manifest, evidence, outdir, records

    def _run(self, root):
        module = self._load_splitter()
        seed_fasta, manifest, evidence, outdir, records = self._write_fixture(root)
        result = module.split_ephaz_seeds(
            seed_fasta=seed_fasta,
            manifest=manifest,
            curated_evidence=evidence,
            outdir=outdir,
            short_review_tsv=outdir / "ePhaZ_short_sequence_review.tsv",
        )
        return module, result, seed_fasta, manifest, evidence, outdir, records

    @staticmethod
    def _accessions(path):
        return {
            line[1:].split("|", 1)[0].strip()
            for line in path.read_text(encoding="utf-8").splitlines()
            if line.startswith(">")
        }

    def test_layers_are_disjoint_and_cover_every_input_accession(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _, _, _, _, _, outdir, records = self._run(root)
            core = self._accessions(outdir / "ePhaZ_curated_core.faa")
            broad = self._accessions(outdir / "ePhaZ_broad_discovery.faa")
            original = {accession for accession, _, _ in records}
            self.assertTrue(core)
            self.assertTrue(broad)
            self.assertTrue(core.isdisjoint(broad))
            self.assertEqual(core | broad, original)
            self.assertEqual(len(original), 66)

    def test_all_61_short_sequences_are_retained_and_reviewed(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _, _, _, _, _, outdir, _ = self._run(root)
            core = self._accessions(outdir / "ePhaZ_curated_core.faa")
            broad = self._accessions(outdir / "ePhaZ_broad_discovery.faa")
            short_ids = {f"S{index:03d}" for index in range(1, 62)}
            self.assertTrue(short_ids.isdisjoint(core))
            self.assertTrue(short_ids.issubset(broad))

            with (outdir / "ePhaZ_short_sequence_review.tsv").open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual({row["accession"] for row in rows}, short_ids)
            self.assertEqual(len(rows), 61)
            self.assertTrue(all(row.get("layer") == "ePhaZ_broad_discovery" for row in rows))
            self.assertTrue(all(row.get("review_status") in {"review_required", "architecture_pending"} for row in rows))

    def test_curated_core_requires_evidence_architecture_and_length(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _, _, _, _, _, outdir, _ = self._run(root)
            core = self._accessions(outdir / "ePhaZ_curated_core.faa")
            broad = self._accessions(outdir / "ePhaZ_broad_discovery.faa")
            self.assertEqual(core, {"C001", "C002", "C003", "C004"})
            self.assertIn("B001", broad)
            self.assertNotIn("B001", core)

            with (outdir / "ePhaZ_layer_manifest.tsv").open(encoding="utf-8", newline="") as handle:
                rows = {row["accession"]: row for row in csv.DictReader(handle, delimiter="\t")}
            for accession in core:
                self.assertEqual(rows[accession]["layer"], "ePhaZ_curated_core")
                self.assertEqual(rows[accession]["evidence_level"], "experimental")
                self.assertEqual(rows[accession]["architecture_status"], "typical")
                self.assertEqual(rows[accession]["completeness_status"], "complete")
                self.assertGreaterEqual(int(rows[accession]["length"]), 200)
            self.assertEqual(rows["B001"]["layer"], "ePhaZ_broad_discovery")

    def test_missing_manifest_evidence_or_sequence_fails_closed(self):
        module = self._load_splitter()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            seed_fasta, manifest, evidence, outdir, _ = self._write_fixture(root)
            for kwargs in (
                {"seed_fasta": root / "missing.faa", "manifest": manifest, "curated_evidence": evidence},
                {"seed_fasta": seed_fasta, "manifest": root / "missing.tsv", "curated_evidence": evidence},
                {"seed_fasta": seed_fasta, "manifest": manifest, "curated_evidence": root / "missing_evidence.tsv"},
            ):
                with self.assertRaises(Exception):
                    module.split_ephaz_seeds(outdir=outdir, **kwargs)
                self.assertFalse(outdir.exists(), "failed validation must not leave partial output")

    def test_duplicate_or_invalid_sequence_fails_closed(self):
        module = self._load_splitter()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            seed_fasta, manifest, evidence, outdir, _ = self._write_fixture(root)
            original = seed_fasta.read_text(encoding="utf-8")
            seed_fasta.write_text(original + ">C001|duplicate\nACDE\n", encoding="utf-8")
            with self.assertRaises(Exception):
                module.split_ephaz_seeds(
                    seed_fasta=seed_fasta,
                    manifest=manifest,
                    curated_evidence=evidence,
                    outdir=outdir,
                )
            self.assertFalse(outdir.exists())

            seed_fasta.write_text(original.replace("ACDE", "ACD1", 1), encoding="utf-8")
            with self.assertRaises(Exception):
                module.split_ephaz_seeds(
                    seed_fasta=seed_fasta,
                    manifest=manifest,
                    curated_evidence=evidence,
                    outdir=outdir,
                )
            self.assertFalse(outdir.exists())

    def test_sha256_manifest_covers_layer_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _, _, _, _, _, outdir, _ = self._run(root)
            checksum = outdir / "sha256.tsv"
            self.assertTrue(checksum.is_file())
            with checksum.open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            by_path = {row["path"]: row["sha256"] for row in rows}
            expected = {
                "ePhaZ_curated_core.faa",
                "ePhaZ_broad_discovery.faa",
                "ePhaZ_layer_manifest.tsv",
                "ePhaZ_short_sequence_review.tsv",
            }
            self.assertTrue(expected.issubset(by_path))
            for filename in expected:
                self.assertEqual(
                    by_path[filename],
                    hashlib.sha256((outdir / filename).read_bytes()).hexdigest(),
                )

    def test_split_manifest_records_input_hashes_and_criteria(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _, _, seed_fasta, manifest, evidence, outdir, _ = self._run(root)
            provenance = outdir / "ePhaZ_split_manifest.json"
            self.assertTrue(provenance.is_file())
            data = json.loads(provenance.read_text(encoding="utf-8"))
            for key, path in (("seed_fasta", seed_fasta), ("manifest", manifest), ("curated_evidence", evidence)):
                self.assertEqual(data["inputs"][key]["sha256"], hashlib.sha256(path.read_bytes()).hexdigest())
            self.assertEqual(data["criteria"]["minimum_length_for_core"], 200)


if __name__ == "__main__":
    unittest.main()
