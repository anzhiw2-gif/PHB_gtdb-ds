from pathlib import Path

from pipeline.scripts.parallel_extract_sequences import extract_shard


def test_extract_shard_writes_only_requested_records(tmp_path):
    shard = tmp_path / "shard_0001.faa"
    shard.write_text(">GCA_1|p1\nMKT\n>GCA_2|p2\nAAA\n", encoding="utf-8")
    out = tmp_path / "out.faa"
    assert extract_shard(shard, {"GCA_2|p2"}, out) == 1
    assert out.read_text(encoding="utf-8") == ">GCA_2|p2\nAAA\n"
