import importlib.util
import unittest

from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "summarize_phaded_interpro.py"


def load_module():
    spec = importlib.util.spec_from_file_location("summarize_phaded_interpro", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class SummarizePhaDEDInterProTests(unittest.TestCase):
    def test_summary_prioritizes_synthase_like_annotation_over_generic_pha_domain(self):
        module = load_module()
        candidates = [{"accession": "A"}]
        hits = [{"accession": "A", "analysis": "PANTHER", "signature": "PTHR36837", "description": "PHA synthase", "interpro_id": "IPR051321", "interpro_desc": "Polyhydroxyalkanoate synthase"}]
        row = module.summarize(candidates, hits)[0]
        self.assertEqual(row["interpro_secondary_state"], "pha_synthase_like_alternative")

    def test_summary_distinguishes_explicit_generic_and_not_submitted_evidence(self):
        module = load_module()
        candidates = [{"accession": "A"}, {"accession": "B"}, {"accession": "C"}]
        hits = [
            {"accession": "A", "analysis": "Pfam", "signature": "PF10503", "description": "Esterase PHB depolymerase", "interpro_id": "IPR010126", "interpro_desc": "Esterase, PHB depolymerase"},
            {"accession": "B", "analysis": "Pfam", "signature": "PF00756", "description": "Putative esterase", "interpro_id": "", "interpro_desc": ""},
        ]
        rows = module.summarize(candidates, hits, {"C": "invalid_amino_acid"})
        by_accession = {row["accession"]: row for row in rows}
        self.assertEqual(by_accession["A"]["interpro_secondary_state"], "explicit_pha_related_domain")
        self.assertEqual(by_accession["B"]["interpro_secondary_state"], "generic_hydrolase_or_esterase")
        self.assertEqual(by_accession["C"]["interpro_secondary_state"], "not_submitted_tool_input_limit")


if __name__ == "__main__":
    unittest.main()
