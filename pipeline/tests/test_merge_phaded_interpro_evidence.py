import importlib.util
import unittest

from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "merge_phaded_interpro_evidence.py"


def load_module():
    spec = importlib.util.spec_from_file_location("merge_phaded_interpro_evidence", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class InterProMergeTests(unittest.TestCase):
    def test_aggregates_rows_and_marks_supported(self):
        module = load_module()
        result = module.summarize_interpro(
            [{"accession": "A"}, {"accession": "B"}],
            [
                {"accession": "A", "analysis": "Pfam", "signature": "PF10503", "interpro_id": "IPR010126", "date": "14-09-2026", "description": "Esterase PHB depolymerase", "interpro_desc": "Esterase, PHB depolymerase"},
                {"accession": "A", "analysis": "NCBIfam", "signature": "TIGR01840", "interpro_id": "IPR010126", "date": "14-09-2026", "description": "PHB depolymerase family esterase", "interpro_desc": "Esterase, PHB depolymerase"},
            ],
            {},
        )
        self.assertEqual(result[0]["interpro_status"], "interpro_supported")
        self.assertEqual(result[0]["interpro_hit_count"], "2")
        self.assertEqual(result[0]["interpro_secondary_state"], "explicit_pha_related_domain")
        self.assertEqual(result[1]["interpro_status"], "interpro_not_reported")

    def test_tool_exclusion_is_distinct_from_not_reported(self):
        module = load_module()
        result = module.summarize_interpro(
            [{"accession": "excluded"}, {"accession": "missing"}],
            [],
            {"excluded": "invalid_amino_acid"},
        )
        self.assertEqual(result[0]["interpro_status"], "interpro_tool_input_excluded")
        self.assertEqual(result[1]["interpro_status"], "interpro_not_reported")

    def test_duplicate_accessions_fail_closed(self):
        module = load_module()
        with self.assertRaises(ValueError):
            module.summarize_interpro([{"accession": "A"}, {"accession": "A"}], [], {})


if __name__ == "__main__":
    unittest.main()
