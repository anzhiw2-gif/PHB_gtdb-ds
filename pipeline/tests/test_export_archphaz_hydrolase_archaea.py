import importlib.util
import json
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "export_archphaz_hydrolase_archaea.py"


def load_module():
    spec = importlib.util.spec_from_file_location("export_archphaz_hydrolase_archaea", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_export_filters_archaea_and_preserves_exact_fasta_records(tmp_path):
    module = load_module()
    table = tmp_path / "tier1_genome_family.tsv"
    fasta = tmp_path / "ArchPhaZ_hydrolase_tier1.faa"
    out_faa = tmp_path / "archaea.faa"
    out_tsv = tmp_path / "archaea.tsv"
    provenance = tmp_path / "provenance.json"
    table.write_text(
        "genome\tfamily\tcopies\tgtdb_taxonomy\tphylum\tclass\n"
        "GCA_A\tArchPhaZ_hydrolase\t1\td__Archaea;p__Euryarchaeota;c__Methanobacteria\tEuryarchaeota\tMethanobacteria\n"
        "GCA_B\tArchPhaZ_hydrolase\t1\td__Bacteria;p__Pseudomonadota;c__Gamma\tPseudomonadota\tGamma\n"
        "GCA_A\tePhaZ\t1\td__Archaea;p__Euryarchaeota;c__Methanobacteria\tEuryarchaeota\tMethanobacteria\n",
        encoding="utf-8",
    )
    fasta.write_text(
        ">GCA_A|locus_1\nMARGIN\n>GCA_B|locus_2\nBACT\n",
        encoding="utf-8",
    )

    summary = module.export_archaea(table, fasta, out_faa, out_tsv, provenance)

    assert summary["records"] == 1
    assert summary["genomes"] == 1
    assert out_faa.read_text(encoding="utf-8") == ">GCA_A|locus_1\nMARGIN\n"
    assert out_tsv.read_text(encoding="utf-8").splitlines() == [
        "genome\tfamily\tlocus\tcopies\tgtdb_taxonomy\tphylum\tclass",
        "GCA_A\tArchPhaZ_hydrolase\tlocus_1\t1\td__Archaea;p__Euryarchaeota;c__Methanobacteria\tEuryarchaeota\tMethanobacteria",
    ]
    payload = json.loads(provenance.read_text(encoding="utf-8"))
    assert payload["filter"]["domain"] == "d__Archaea"
    assert payload["records"] == 1


def test_export_rejects_missing_table_match(tmp_path):
    module = load_module()
    table = tmp_path / "tier1.tsv"
    fasta = tmp_path / "tier1.faa"
    table.write_text(
        "genome\tfamily\tcopies\tgtdb_taxonomy\tphylum\tclass\n"
        "GCA_A\tArchPhaZ_hydrolase\t1\td__Archaea;p__X;c__Y\tX\tY\n",
        encoding="utf-8",
    )
    fasta.write_text(">GCA_B|locus_1\nSEQ\n", encoding="utf-8")
    try:
        module.export_archaea(table, fasta, tmp_path / "o.faa", tmp_path / "o.tsv", tmp_path / "p.json")
    except ValueError as exc:
        assert "missing FASTA records" in str(exc)
    else:
        raise AssertionError("expected missing FASTA record validation failure")
