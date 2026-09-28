import tempfile
import unittest
from pathlib import Path

from pipeline.scripts.prepare_phaded_competition_review import prepare


class CompetitionReviewInputTests(unittest.TestCase):
    def test_requires_exact_four_buckets(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            priority = root / "priority.tsv"
            dossier = root / "dossier.tsv"
            fasta = root / "candidates.faa"
            reference = root / "reference.faa"
            priority.write_text("accession\tlayer\nA\tePhaZ_competition_review\n", encoding="utf-8")
            dossier.write_text("accession\tprovisional_label\nA\t\n", encoding="utf-8")
            fasta.write_text(">A\nMKT\n", encoding="ascii")
            reference.write_text(">R\nMKT\n", encoding="ascii")
            with self.assertRaisesRegex(ValueError, "unexpected bucket counts"):
                prepare(priority, dossier, fasta, reference, root / "out")

    def test_rejects_duplicate_bucket_membership(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            priority = root / "priority.tsv"
            dossier = root / "dossier.tsv"
            fasta = root / "candidates.faa"
            reference = root / "reference.faa"
            priority.write_text("accession\tlayer\nA\tePhaZ_competition_review\n", encoding="utf-8")
            dossier.write_text("accession\tprovisional_label\nA\tA_retain_high_priority_architecture\n", encoding="utf-8")
            fasta.write_text(">A\nMKT\n", encoding="ascii")
            reference.write_text(">R\nMKT\n", encoding="ascii")
            with self.assertRaisesRegex(ValueError, "unexpected bucket counts"):
                prepare(priority, dossier, fasta, reference, root / "out")


if __name__ == "__main__":
    unittest.main()
