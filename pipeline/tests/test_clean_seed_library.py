"""Tests for the immutable, accession-aware seed library cleaner."""
import csv
import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "clean_seed_library.py"
SPEC = importlib.util.spec_from_file_location("clean_seed_library", SCRIPT)
clean_seed_library = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(clean_seed_library)


EXCLUDED = clean_seed_library.EXCLUDED_ACCESSIONS
FAMILIES = clean_seed_library.REQUIRED_FAMILIES


class CleanSeedLibraryTests(unittest.TestCase):
    def _fixture(self, root):
        seed_dir = root / "seeds"
        seed_dir.mkdir()
        (seed_dir / "OH.faa").write_text(
            ">Q79F77|nylC\nAAAA\n>Q1EPR4|nylC\nCCCC\n>Q1EPR5|nylC\nGGGG\n>P07061|nylB\nTTTT\n>P07062|nylB\nACAC\n>POS1|hydrolase\nACGT\n",
            encoding="utf-8",
        )
        bdh_records = [(acc, "ACGT") for acc in sorted(EXCLUDED - {"Q79F77", "Q1EPR4", "Q1EPR5", "P07061", "P07062"})]
        bdh_records.append(("KEEP1", "TTTT"))
        (seed_dir / "BdhA.faa").write_text(
            "".join(f">{acc}|BDH1\n{sequence}\n" for acc, sequence in bdh_records),
            encoding="utf-8",
        )
        for family in FAMILIES:
            if family in {"OH", "BdhA"}:
                continue
            (seed_dir / f"{family}.faa").write_text(
                f">{family.upper()}_POS|{family}\nACGT\n", encoding="utf-8"
            )
        (seed_dir / "ArchPhaZ_hydrolase.faa").write_text(
            ">DERIVED_POS|ArchPhaZ_hydrolase\nACGT\n", encoding="utf-8"
        )
        manifest = root / "v2_manifest.tsv"
        rows = [
            ("Q79F77", "OH"), ("Q1EPR4", "OH"), ("Q1EPR5", "OH"),
            ("P07061", "OH"), ("P07062", "OH"), ("POS1", "OH"),
            *((acc, "BdhA") for acc in sorted(EXCLUDED - {"Q79F77", "Q1EPR4", "Q1EPR5", "P07061", "P07062"})),
            ("KEEP1", "BdhA"),
        ]
        rows.extend((f"{family.upper()}_POS", family) for family in FAMILIES if family not in {"OH", "BdhA"})
        with manifest.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(("accession", "family", "evidence"))
            writer.writerows((acc, family, "fixture") for acc, family in rows)
        negative = root / "negative.faa"
        negative.write_text(
            "".join(f">{acc}|negative\nACGT\n" for acc in sorted(EXCLUDED)),
            encoding="utf-8",
        )
        return seed_dir, manifest, negative

    def test_excludes_accessions_preserves_negative_and_records_hashes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            seed_dir, manifest, negative = self._fixture(root)
            output = root / "runs" / "run1" / "inputs" / "seeds_clean"
            result = clean_seed_library.clean_library(seed_dir, manifest, negative, output)

            self.assertEqual(result["excluded_accessions_count"], 18)
            self.assertNotIn("Q79F77", (output / "OH.faa").read_text(encoding="utf-8"))
            self.assertIn(">POS1|hydrolase", (output / "OH.faa").read_text(encoding="utf-8"))
            self.assertNotIn("P29147", (output / "BdhA.faa").read_text(encoding="utf-8"))
            self.assertIn("KEEP1", (output / "BdhA.faa").read_text(encoding="utf-8"))
            self.assertEqual(
                set(line.split("\t")[0] for line in (output / "excluded_accessions.tsv").read_text(encoding="utf-8").splitlines()[1:]),
                EXCLUDED,
            )
            self.assertEqual(
                set(line[1:].split("|", 1)[0] for line in (output / "negative.faa").read_text(encoding="utf-8").splitlines() if line.startswith(">")),
                EXCLUDED,
            )
            hashes = (output / "sha256.tsv").read_text(encoding="utf-8")
            self.assertIn("cleaning_manifest.json", hashes)
            manifest_data = json.loads((output / "cleaning_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest_data["excluded_accessions_count"], 18)
            self.assertEqual(manifest_data["outputs"]["OH.faa"]["sha256"], hashlib.sha256((output / "OH.faa").read_bytes()).hexdigest())
            self.assertEqual(manifest_data["negative_fasta"]["excluded_intersection_count"], 18)
            self.assertTrue((output / "v2_manifest.tsv").is_file())
            self.assertEqual(clean_seed_library.verify_sha256_tsv(output), [])

    def test_missing_accession_or_input_fails_without_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            seed_dir, manifest, negative = self._fixture(root)
            manifest.write_text(manifest.read_text(encoding="utf-8") + "\nMISSING\tOH\tfixture\n", encoding="utf-8")
            output = root / "runs" / "run1" / "inputs" / "seeds_clean"
            with self.assertRaises(clean_seed_library.CleaningError):
                clean_seed_library.clean_library(seed_dir, manifest, negative, output)
            self.assertFalse(output.exists())
            with self.assertRaises(clean_seed_library.CleaningError):
                clean_seed_library.clean_library(seed_dir / "missing", manifest, negative, root / "runs" / "run2" / "inputs" / "seeds_clean")

    def test_refuses_output_inside_input_seed_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            seed_dir, manifest, negative = self._fixture(root)
            with self.assertRaises(clean_seed_library.CleaningError):
                clean_seed_library.clean_library(seed_dir, manifest, negative, seed_dir / "cleaned")

    def test_output_must_be_under_runs_runid_inputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            seed_dir, manifest, negative = self._fixture(root)
            with self.assertRaises(clean_seed_library.CleaningError):
                clean_seed_library.clean_library(seed_dir, manifest, negative, root / "run" / "inputs" / "seeds_clean")
            with self.assertRaises(clean_seed_library.CleaningError):
                clean_seed_library.clean_library(seed_dir, manifest, negative, root / "runs" / ".." / "escape" / "inputs" / "seeds_clean")

    def test_all_nine_families_are_required_and_cleaned_nonempty(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            seed_dir, manifest, negative = self._fixture(root)
            (seed_dir / "PhaC.faa").unlink()
            with self.assertRaises(clean_seed_library.CleaningError):
                clean_seed_library.clean_library(seed_dir, manifest, negative, root / "runs" / "r1" / "inputs" / "seeds_clean")

            second = root / "second"
            second.mkdir()
            seed_dir, manifest, negative = self._fixture(second)
            (seed_dir / "PhaC.faa").write_text(">PhaC_POS|PhaC\n", encoding="utf-8")
            with self.assertRaises(clean_seed_library.CleaningError):
                clean_seed_library.clean_library(seed_dir, manifest, negative, second / "runs" / "r1" / "inputs" / "seeds_clean")

    def test_manifest_family_sets_must_match_each_fasta(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            seed_dir, manifest, negative = self._fixture(root)
            text = manifest.read_text(encoding="utf-8")
            text = text.replace("PHAC_POS\tPhaC", "PHAC_POS\tPhaJ")
            manifest.write_text(text, encoding="utf-8")
            with self.assertRaises(clean_seed_library.CleaningError):
                clean_seed_library.clean_library(seed_dir, manifest, negative, root / "runs" / "r1" / "inputs" / "seeds_clean")

    def test_explicit_manifest_output_name_is_safe_and_recorded(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            seed_dir, manifest, negative = self._fixture(root)
            output = root / "runs" / "r2" / "inputs" / "seeds_clean"
            result = clean_seed_library.clean_library(
                seed_dir, manifest, negative, output, manifest_output_name="clean.tsv"
            )
            self.assertTrue((output / "clean.tsv").is_file())
            self.assertEqual(result["manifest_output"], "clean.tsv")
            with self.assertRaises(clean_seed_library.CleaningError):
                clean_seed_library.clean_library(
                    seed_dir, manifest, negative,
                    root / "runs" / "r3" / "inputs" / "seeds_clean",
                    manifest_output_name="../escape.tsv",
                )

    def test_negative_must_contain_exactly_the_18_exclusions(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            seed_dir, manifest, negative = self._fixture(root)
            negative.write_text(negative.read_text(encoding="utf-8") + ">EXTRA|negative\nACGT\n", encoding="utf-8")
            with self.assertRaises(clean_seed_library.CleaningError):
                clean_seed_library.clean_library(seed_dir, manifest, negative, root / "runs" / "r1" / "inputs" / "seeds_clean")

    def test_duplicate_negative_accession_fails_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            seed_dir, manifest, negative = self._fixture(root)
            negative.write_text(negative.read_text(encoding="utf-8") + ">P29147|duplicate\nACGT\n", encoding="utf-8")
            with self.assertRaises(clean_seed_library.CleaningError):
                clean_seed_library.clean_library(seed_dir, manifest, negative, root / "runs" / "r1" / "inputs" / "seeds_clean")

    def test_checksum_verification_reports_tampering(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            seed_dir, manifest, negative = self._fixture(root)
            output = root / "runs" / "r1" / "inputs" / "seeds_clean"
            clean_seed_library.clean_library(seed_dir, manifest, negative, output)
            (output / "OH.faa").write_text((output / "OH.faa").read_text(encoding="utf-8") + "\n", encoding="utf-8")
            self.assertTrue(clean_seed_library.verify_sha256_tsv(output))

    def test_empty_fasta_record_fails_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            seed_dir, manifest, negative = self._fixture(root)
            (seed_dir / "PhaC.faa").write_text(">PHAC_POS|PhaC\n", encoding="utf-8")
            output = root / "runs" / "r-empty" / "inputs" / "seeds_clean"
            with self.assertRaises(clean_seed_library.CleaningError):
                clean_seed_library.clean_library(seed_dir, manifest, negative, output)

    @unittest.skipUnless(hasattr(Path, "symlink_to"), "symlink support unavailable")
    def test_output_path_symlink_ancestor_fails_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            seed_dir, manifest, negative = self._fixture(root)
            outside = root / "outside"
            outside.mkdir()
            runs = root / "runs"
            try:
                runs.symlink_to(outside, target_is_directory=True)
            except OSError as error:
                self.skipTest(f"symlink creation unavailable: {error}")
            output = runs / "r-link" / "inputs" / "seeds_clean"
            with self.assertRaises(clean_seed_library.CleaningError):
                clean_seed_library.clean_library(seed_dir, manifest, negative, output)

    def test_duplicate_accession_in_training_fasta_fails_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            seed_dir, manifest, negative = self._fixture(root)
            with (seed_dir / "BdhA.faa").open("a", encoding="utf-8") as handle:
                handle.write(">KEEP1|duplicate\nAAAA\n")
            with self.assertRaises(clean_seed_library.CleaningError):
                clean_seed_library.clean_library(seed_dir, manifest, negative, root / "runs" / "run1" / "inputs" / "seeds_clean")


if __name__ == "__main__":
    unittest.main()
