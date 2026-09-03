import importlib.util
import json
from pathlib import Path


def load_module():
    path = Path(__file__).parents[1] / "scripts" / "build_ephaz_independent_panels.py"
    spec = importlib.util.spec_from_file_location("build_ephaz_independent_panels", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_eval(path: Path, value: str):
    path.write_text(json.dumps({"value": value}), encoding="utf-8")


def prepare_protected_sources(tmp_path: Path):
    bridge = tmp_path / "runs" / "20260830_ephaz_bridge_curation_02"
    (bridge / "results").mkdir(parents=True)
    (bridge / "inputs").mkdir(parents=True)
    (tmp_path / "pipeline" / "seeds" / "controls").mkdir(parents=True)
    (bridge / "results" / "ephaz_bridge_candidate.faa").write_text(
        ">Q51871|bridge\nMBBBB\n>Q5SLU4|bridge\nMCCCCC\n", encoding="utf-8"
    )
    (bridge / "inputs" / "ephaz_curated_core.faa").write_text(
        ">B2NHN2|core\nMNNNN\n", encoding="utf-8"
    )
    (tmp_path / "pipeline" / "seeds" / "controls" / "positive.faa").write_text(
        ">P12625|control\nMPPPP\n", encoding="utf-8"
    )


def test_builder_separates_panels_and_preserves_evidence(tmp_path):
    module = load_module()
    prepare_protected_sources(tmp_path)
    run = tmp_path / "run"
    inputs = run / "inputs"
    inputs.mkdir(parents=True)
    (run / "logs").mkdir()
    (run / "results").mkdir()
    write_eval(inputs / "POS_eval.json", ">sp|POS|test positive\nMKTAA\n")
    write_eval(inputs / "NEG_eval.json", ">sp|NEG|test negative\nMKKLL\n")
    manifest = inputs / "candidate_manifest.tsv"
    manifest.write_text(
        "\t".join(module.MANIFEST_FIELDS)
        + "\n"
        + "\t".join(
            [
                "POS",
                "independent_experimental_positive",
                "PHB",
                "Testus positive",
                "5",
                "true",
                "PMID:1",
                "doi:10/test",
                "uniprot_fasta",
                "POS_eval.json",
                "complete secreted PHB depolymerase",
            ]
        )
        + "\n"
        + "\t".join(
            [
                "NEG",
                "intracellular_non_ephaz_negative",
                "PHB",
                "Testus negative",
                "5",
                "false",
                "PMID:2",
                "doi:10/test2",
                "uniprot_fasta",
                "NEG_eval.json",
                "no signal peptide; intracellular",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    summary = module.build_panels(run, manifest)
    assert summary["counts"] == {
        "independent_experimental_positive": 1,
            "mcl_pha_experimental_positive": 0,
            "mcl_pha_non_phb_negative": 0,
        "expanded_ephaz_near_neighbor_negative": 1,
        "intracellular_non_ephaz_negative": 1,
        "annotation_only_near_neighbor_negative": 0,
        "fragment_or_incomplete_negative": 0,
    }
    assert (run / "results" / "independent_experimental_positive.faa").read_text().startswith(">POS|")
    rows = (run / "results" / "ephaz_panel_evidence.tsv").read_text(encoding="utf-8").splitlines()
    assert "PMID:1" in rows[1]
    assert json.loads((run / "input_contract.json").read_text(encoding="utf-8"))["protected_accessions"]


def test_builder_rejects_bridge_overlap(tmp_path):
    module = load_module()
    prepare_protected_sources(tmp_path)
    run = tmp_path / "run"
    inputs = run / "inputs"
    inputs.mkdir(parents=True)
    (run / "logs").mkdir()
    (run / "results").mkdir()
    write_eval(inputs / "bridge_eval.json", ">sp|Q51871|bridge\nMKTAA\n")
    manifest = inputs / "candidate_manifest.tsv"
    manifest.write_text(
        "\t".join(module.MANIFEST_FIELDS)
        + "\n"
        + "\t".join(
            [
                "Q51871",
                "independent_experimental_positive",
                "PHB",
                "Bridge species",
                "5",
                "false",
                "PMID:1",
                "",
                "uniprot_fasta",
                "bridge_eval.json",
                "bridge sequence",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    try:
        module.build_panels(run, manifest)
    except module.PanelBuildError as exc:
        assert "protected" in str(exc).lower()
    else:
        raise AssertionError("bridge overlap must fail closed")


def test_builder_rejects_missing_raw_response(tmp_path):
    module = load_module()
    prepare_protected_sources(tmp_path)
    run = tmp_path / "run"
    inputs = run / "inputs"
    inputs.mkdir(parents=True)
    (run / "logs").mkdir()
    (run / "results").mkdir()
    manifest = inputs / "candidate_manifest.tsv"
    manifest.write_text(
        "\t".join(module.MANIFEST_FIELDS)
        + "\n"
        + "\t".join(
            [
                "MISS",
                "annotation_only_near_neighbor_negative",
                "NA",
                "Missingus",
                "5",
                "false",
                "",
                "",
                "uniprot_fasta",
                "missing_eval.json",
                "annotation only",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    try:
        module.build_panels(run, manifest)
    except module.PanelBuildError as exc:
        assert "missing" in str(exc).lower()
    else:
        raise AssertionError("missing raw response must fail closed")


def test_sequence_hashes_use_first_header_field_and_keep_all_records(tmp_path):
    module = load_module()
    fasta = tmp_path / "protected.faa"
    fasta.write_text(">Q51871|bridge|x\nMKTAA\n>Q5SLU4|bridge|x\nMKKLL\n", encoding="utf-8")
    hashes = module._sequence_hashes(fasta)
    assert set(hashes) == {"Q51871", "Q5SLU4"}


def test_builder_rejects_substring_accession_and_multirecord_response(tmp_path):
    module = load_module()
    prepare_protected_sources(tmp_path)
    run = tmp_path / "run"
    inputs = run / "inputs"
    inputs.mkdir(parents=True)
    (run / "logs").mkdir()
    (run / "results").mkdir()
    write_eval(inputs / "BAD_eval.json", ">sp|O247190|wrong\nMKTAA\n>sp|SECOND|extra\nMKKLL\n")
    manifest = inputs / "candidate_manifest.tsv"
    manifest.write_text(
        "\t".join(module.MANIFEST_FIELDS) + "\n" + "\t".join([
            "O24719", "annotation_only_near_neighbor_negative", "NA", "Testus", "5", "false", "", "",
            "uniprot_fasta", "BAD_eval.json", "annotation only",
        ]) + "\n", encoding="utf-8")
    try:
        module.build_panels(run, manifest)
    except module.PanelBuildError as exc:
        assert "mismatch" in str(exc).lower() or "multiple" in str(exc).lower()
    else:
        raise AssertionError("non-exact or multi-record response must fail closed")


def test_builder_rejects_response_path_traversal(tmp_path):
    module = load_module()
    prepare_protected_sources(tmp_path)
    run = tmp_path / "run"
    inputs = run / "inputs"
    inputs.mkdir(parents=True)
    (run / "logs").mkdir()
    (run / "results").mkdir()
    outside = tmp_path / "outside.json"
    write_eval(outside, ">sp|X|outside\nMKTAA\n")
    manifest = inputs / "candidate_manifest.tsv"
    manifest.write_text(
        "\t".join(module.MANIFEST_FIELDS) + "\n" + "\t".join([
            "X", "annotation_only_near_neighbor_negative", "NA", "Testus", "5", "false", "", "",
            "uniprot_fasta", "../outside.json", "annotation only",
        ]) + "\n", encoding="utf-8")
    try:
        module.build_panels(run, manifest)
    except module.PanelBuildError as exc:
        assert "path" in str(exc).lower() or "outside" in str(exc).lower()
    else:
        raise AssertionError("response path traversal must fail closed")
