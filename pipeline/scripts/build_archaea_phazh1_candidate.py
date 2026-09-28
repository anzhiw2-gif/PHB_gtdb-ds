#!/usr/bin/env python3
"""Build the archaeal PhaZh1-like candidate HMM exactly once per run."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
from pathlib import Path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def build_commands(seed: Path, alignment: Path, hmm: Path, mafft: str = "mafft", hmmbuild: str = "hmmbuild") -> tuple[list[str], list[str]]:
    return ([mafft, "--auto", str(seed)], [hmmbuild, "--amino", str(hmm), str(alignment)])


def build(seed: Path, outdir: Path, *, mafft: str = "mafft", hmmbuild: str = "hmmbuild", logdir: Path | None = None) -> dict[str, object]:
    if not seed.is_file() or seed.stat().st_size == 0:
        raise ValueError(f"seed panel missing or empty: {seed}")
    outdir.mkdir(parents=True, exist_ok=True)
    logdir = logdir or outdir.parent / "logs"
    logdir.mkdir(parents=True, exist_ok=True)
    alignment = outdir / "archaea_PhaZh1_like_candidate.sto"
    hmm = outdir / "archaea_PhaZh1_like_candidate.hmm"
    mafft_cmd, hmmbuild_cmd = build_commands(seed, alignment, hmm, mafft, hmmbuild)
    with alignment.open("w", encoding="utf-8", newline="\n") as handle:
        result = subprocess.run(mafft_cmd, check=True, stdout=handle, stderr=subprocess.PIPE, text=True)
    (logdir / "mafft.stderr.log").write_text(result.stderr or "", encoding="utf-8")
    with (logdir / "hmmbuild.log").open("w", encoding="utf-8", newline="\n") as log:
        subprocess.run(hmmbuild_cmd, check=True, stdout=log, stderr=subprocess.STDOUT, text=True)
    if not hmm.is_file() or hmm.stat().st_size == 0:
        raise RuntimeError("hmmbuild produced no HMM")
    metadata = {
        "status": "candidate-only",
        "seed_sha256": sha256_file(seed),
        "alignment_sha256": sha256_file(alignment),
        "hmm_sha256": sha256_file(hmm),
        "commands": {"mafft": mafft_cmd, "hmmbuild": hmmbuild_cmd},
    }
    (outdir / "build_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    return metadata


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--mafft", default=shutil.which("mafft") or "mafft")
    parser.add_argument("--hmmbuild", default=shutil.which("hmmbuild") or "hmmbuild")
    args = parser.parse_args()
    print(json.dumps(build(args.seed, args.outdir, mafft=args.mafft, hmmbuild=args.hmmbuild), indent=2))


if __name__ == "__main__":
    main()
