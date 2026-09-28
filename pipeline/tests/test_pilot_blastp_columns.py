"""Regression tests for the pilot blastp artifact column mapping (test-first).

Why this test exists
--------------------
The with-lipase structural-validation pilot ran

    blastp -outfmt "6 qseqid pident length evalue bitscore" -evalue 1e-3 ...

(5 fields: index 0 qseqid, 1 pident, 2 length, **3 evalue**, **4 bitscore**) but
counted hits with ``float(r[4])``, i.e. it applied E-value thresholds to the
*bitscore* column.  Every threshold therefore returned 0 and the run recorded the
false claim "0/1,001 sequences reach E<1e-5", which was then used as evidence that
"sequence similarity cannot narrow the deferred layer".  The same column mix-up was
made once before in this task (handoff §11.13 point 6), so the mapping is now pinned
by a machine-checked invariant instead of being re-read by eye.

The invariant is unambiguous because the run passed ``-evalue 1e-3``: every reported
row must have an E-value <= 1e-3, and a bitscore column can never satisfy that.
Nothing here judges biology; it only asserts that the fields are the fields.

Artifacts live under ``runs/`` (git-ignored), so every test skips cleanly when the
run directory is not present.
"""

import bisect
import json
import unittest
from pathlib import Path

RUN = (
    Path(__file__).resolve().parents[2]
    / "runs"
    / "20260919_phaded_pool_external_evidence_01"
)
PILOT_TSV = RUN / "logs" / "pilot_server" / "pilot_vs_8ynv.tsv"
PILOT_FAA = RUN / "logs" / "pilot_server" / "pilot_sequences.faa"
RECORD = RUN / "logs" / "structure_verification_recon_pilot_20260920.json"
DEFERRED_TSV = (
    Path(__file__).resolve().parents[2]
    / "runs"
    / "20260919_phaded_with_lipase_deferred_archive_01"
    / "results"
    / "pool_external_with_lipase_deferred.tsv"
)

# the blastp command's own cutoff, and the single subject it searched (8YNV chain A)
BLASTP_EVALUE_CUTOFF = 1e-3

# values recomputed from the artifact on 2026-09-20 (see plan/verify_pilot_blastp_columns.py)
EXPECTED = {
    "rows": 987,
    "evalue_min": 1.16e-122,
    "evalue_lt_1e-5": 954,
    "evalue_lt_1e-10": 649,
    "evalue_lt_1e-30": 2,
    "evalue_lt_1e-50": 1,
    "evalue_lt_1e-100": 1,
    "pident_median": 24.254,
    "pident_max": 54.138,
}


def _rows():
    if not PILOT_TSV.is_file():
        return []
    return [
        line.split("\t")
        for line in PILOT_TSV.read_text().splitlines()
        if line.strip()
    ]


class BlastpColumnMappingTest(unittest.TestCase):
    """Index 3 is the E-value and index 4 is the bitscore - pinned, not assumed."""

    def setUp(self):
        self.rows = _rows()
        if not self.rows:
            self.skipTest(f"pilot artifact absent: {PILOT_TSV}")

    def test_every_row_has_five_fields(self):
        self.assertEqual({len(r) for r in self.rows}, {5})

    def test_index_3_is_the_evalue(self):
        """The -evalue 1e-3 cutoff must bound index 3 on every row."""
        values = [float(r[3]) for r in self.rows]
        self.assertLessEqual(max(values), BLASTP_EVALUE_CUTOFF)

    def test_index_4_is_the_bitscore_and_cannot_be_the_evalue(self):
        values = [float(r[4]) for r in self.rows]
        self.assertGreaterEqual(min(values), 10.0)
        # this is the bug: applying the E-value cutoff to the bitscore column yields 0
        self.assertFalse(all(v <= BLASTP_EVALUE_CUTOFF for v in values))
        self.assertEqual(
            sum(1 for v in values if v <= BLASTP_EVALUE_CUTOFF),
            0,
            "the bitscore column reproduces the erroneous 0-hit count exactly",
        )

    def test_corrected_statistics_from_the_same_artifact(self):
        evalues = sorted(float(r[3]) for r in self.rows)
        pident = sorted(float(r[1]) for r in self.rows)
        n = len(self.rows)

        def e_lt(t):
            return bisect.bisect_left(evalues, t)

        self.assertEqual(n, EXPECTED["rows"])
        self.assertEqual(e_lt(1e-5), EXPECTED["evalue_lt_1e-5"])
        self.assertEqual(e_lt(1e-10), EXPECTED["evalue_lt_1e-10"])
        self.assertEqual(e_lt(1e-30), EXPECTED["evalue_lt_1e-30"])
        self.assertEqual(e_lt(1e-50), EXPECTED["evalue_lt_1e-50"])
        self.assertEqual(e_lt(1e-100), EXPECTED["evalue_lt_1e-100"])
        self.assertAlmostEqual(evalues[0] / EXPECTED["evalue_min"], 1.0, places=3)
        self.assertAlmostEqual(pident[n // 2], EXPECTED["pident_median"], places=3)
        self.assertAlmostEqual(pident[-1], EXPECTED["pident_max"], places=3)


class CorrectionRecordTest(unittest.TestCase):
    """The machine-readable record must carry the correction, not the false claim."""

    def setUp(self):
        if not PILOT_TSV.is_file() or not RECORD.is_file():
            self.skipTest(f"artifact or record absent under {RUN}")
        self.record = json.loads(RECORD.read_text())
        self.rows = _rows()

    def test_erroneous_block_is_marked_superseded(self):
        block = self.record["pilot"]["blastp_vs_8ynv"]
        self.assertIs(block.get("superseded"), True)
        self.assertEqual(block.get("hits_evalue_lt_1e-5"), 0)  # kept as history

    def test_corrections_block_names_the_bug(self):
        corrections = self.record.get("corrections")
        self.assertTrue(corrections, "record must carry a corrections block")
        entry = next(
            c for c in corrections if "blastp_vs_8ynv" in str(c.get("field", ""))
        )
        self.assertIn("r[4]", entry["root_cause"])
        self.assertIn("evalue", entry["root_cause"])
        self.assertTrue(entry.get("script"))

    def test_corrected_values_match_the_artifact(self):
        corrections = self.record["corrections"]
        corrected = next(
            c["corrected"]
            for c in corrections
            if "blastp_vs_8ynv" in str(c.get("field", ""))
        )
        evalues = sorted(float(r[3]) for r in self.rows)
        self.assertEqual(corrected["hits_evalue_lt_1e-5"], EXPECTED["evalue_lt_1e-5"])
        self.assertEqual(
            corrected["hits_evalue_lt_1e-5"],
            bisect.bisect_left(evalues, 1e-5),
            "the recorded correction must equal the value recomputed from the artifact",
        )
        self.assertIn("subject_search_space", corrected)

    def test_conclusion_is_restated_not_silently_kept(self):
        """The old conclusion stays, but the corrected interpretation is recorded."""
        citations = json.dumps(self.record.get("corrections"), ensure_ascii=False)
        self.assertIn("representative", citations.lower())


class PilotSampleRepresentativenessTest(unittest.TestCase):
    """A pilot may only speak about the population if it was drawn from it."""

    def test_sample_is_mid_population_not_top_tail(self):
        if not PILOT_FAA.is_file() or not DEFERRED_TSV.is_file():
            self.skipTest("pilot sample or deferred population absent")
        sample_ids = {
            line[1:].split()[0]
            for line in PILOT_FAA.read_text().splitlines()
            if line.startswith(">")
        }
        sample_rows = []
        with DEFERRED_TSV.open() as fh:
            header = fh.readline().rstrip("\n").split("\t")
            i_e = header.index("discovery_best_evalue")
            i_p = header.index("protein_id")
            for row_index, line in enumerate(fh):
                f = line.rstrip("\n").split("\t")
                if f[i_p] in sample_ids:
                    sample_rows.append((row_index, float(f[i_e])))

        self.assertEqual(len(sample_rows), len(sample_ids))
        median_e = sorted(e for _i, e in sample_rows)[len(sample_rows) // 2]

        total = 0
        below = 0
        with DEFERRED_TSV.open() as fh:
            fh.readline()
            for line in fh:
                f = line.rstrip("\n").split("\t")
                total += 1
                try:
                    if float(f[i_e]) < median_e:
                        below += 1
                except ValueError:
                    pass

        rank_ratio = below / total
        self.assertGreater(rank_ratio, 0.25)
        self.assertLess(rank_ratio, 0.75)

    def test_sample_design_is_an_even_stride_covering_the_population(self):
        """A systematic stride is only unbiased if no region of the file is skipped."""
        if not PILOT_FAA.is_file() or not DEFERRED_TSV.is_file():
            self.skipTest("pilot sample or deferred population absent")
        sample_ids = {
            line[1:].split()[0]
            for line in PILOT_FAA.read_text().splitlines()
            if line.startswith(">")
        }
        indices = []
        total = 0
        with DEFERRED_TSV.open() as fh:
            header = fh.readline().rstrip("\n").split("\t")
            i_p = header.index("protein_id")
            for row_index, line in enumerate(fh):
                total += 1
                if line.rstrip("\n").split("\t")[i_p] in sample_ids:
                    indices.append(row_index)

        self.assertEqual(indices[0], 0)
        strides = {b - a for a, b in zip(indices, indices[1:])}
        self.assertEqual(len(strides), 1, f"stride is not constant: {sorted(strides)[:5]}")
        stride = strides.pop()
        # the stride must span the population, i.e. no tail region is left unsampled
        self.assertLessEqual(abs(stride * len(indices) - total), stride)
        self.assertGreater(indices[-1], total * 0.95)


if __name__ == "__main__":
    unittest.main()
