#!/usr/bin/env python3
"""Audit the git HISTORY for machine-identity strings, read-only.

The F16 gate checks the WORKING TREE. That gate was itself broken until this
session (``core.quotePath`` made ``git ls-files -z`` return a single element, so
it was a no-op), which means machine-identity strings may have been committed
before the fix landed and would still be reachable in history even though the tree
is clean now.

Rewriting history is a destructive operation needing explicit authorisation, so
this module **measures** instead: it finds every unique blob that violates the
same patterns the gate uses, names the commits carrying each one, and stops there.
It writes nothing unless asked, and changes nothing either way.

Method, in three subprocess calls regardless of repository size: enumerate the
object graph once, ask ``cat-file --batch-check`` for every blob's size, then
stream the text-sized blobs through ``cat-file --batch`` and scan each unique blob
once rather than once per commit.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from pipeline.tests.test_public_repo_safety import IDENTITY_PATTERNS  # noqa: E402

#: The working-tree gate skips its own module, because that module necessarily
#: holds the violating forms as test FIXTURES (it has to, in order to test that
#: they are detected). The history audit must apply the same exemption or it
#: reports the fixture holder as a credential leak - which is exactly what the
#: first version did, producing a false "PEM private key in history" finding.
#: Exempted blobs are COUNTED and reported, never silently dropped.
FIXTURE_HOLDER = "pipeline/tests/test_public_repo_safety.py"


def git_bytes(*arguments: str, workdir: Path = ROOT, stdin: str | None = None) -> bytes:
    """Run git and return raw stdout bytes.

    ``cat-file --batch`` reports sizes in BYTES, so its stream must be sliced as
    bytes: decoding first would shift every offset after the first multi-byte
    character and desynchronise the parser.
    """
    result = subprocess.run(
        ["git", "-c", "core.quotePath=false", *arguments],
        cwd=workdir, capture_output=True, check=False,
        input=None if stdin is None else stdin.encode("utf-8"),
    )
    if result.returncode != 0:
        raise ValueError(
            f"git {' '.join(arguments)} failed ({result.returncode}): "
            f"{result.stderr.decode('utf-8', 'replace').strip()}"
        )
    return result.stdout


def git(*arguments: str, workdir: Path = ROOT, stdin: str | None = None) -> str:
    return git_bytes(*arguments, workdir=workdir, stdin=stdin).decode("utf-8", "replace")


def history_blobs(workdir: Path = ROOT) -> tuple[dict[str, str], int]:
    """{blob_sha: a path it appears at} for every blob reachable from any ref.

    Also returns the number of commits seen, so an empty result cannot be mistaken
    for a clean history.
    """
    seen: dict[str, str] = {}
    commits = 0
    for line in git("rev-list", "--objects", "--all", workdir=workdir).splitlines():
        if not line.strip():
            continue
        parts = line.split(" ", 1)
        if len(parts) == 1:
            commits += 1
            continue
        seen.setdefault(parts[0], parts[1])
    return seen, commits


def scan(workdir: Path = ROOT, *, max_blob_bytes: int = 8 * 1024 * 1024) -> dict:
    blobs, commits = history_blobs(workdir)
    if not blobs:
        raise ValueError("no blobs found; refusing to report a clean history from an empty scan")

    sizes: dict[str, int] = {}
    for line in git(
        "cat-file", "--batch-check",
        workdir=workdir,
        stdin="".join(f"{sha}\n" for sha in blobs),
    ).splitlines():
        parts = line.split()
        if len(parts) == 3 and parts[1] == "blob":
            sizes[parts[0]] = int(parts[2])

    scannable = sorted(
        sha for sha, size in sizes.items() if size <= max_blob_bytes
    )
    skipped = sorted(
        f"{sha[:12]} {blobs[sha]} ({sizes[sha]} bytes)" for sha in sizes
        if sizes[sha] > max_blob_bytes
    )

    violations: list[dict[str, str]] = []
    fixture_holder: list[dict[str, str]] = []
    scanned = 0
    binary = 0
    stream = git_bytes(
        "cat-file", "--batch",
        workdir=workdir,
        stdin="".join(f"{sha}\n" for sha in scannable),
    )
    # Batch output is: b"<sha> blob <size>\n<content>\n", repeated. Parse by byte
    # offsets rather than splitting on newlines, because content is arbitrary text.
    position = 0
    while position < len(stream):
        header_end = stream.find(b"\n", position)
        if header_end < 0:
            break
        header = stream[position:header_end].split()
        if len(header) != 3 or header[1] != b"blob":
            raise ValueError(
                f"unexpected batch header: {stream[position:header_end]!r}"
            )
        sha = header[0].decode("ascii")
        size = int(header[2])
        raw = stream[header_end + 1:header_end + 1 + size]
        position = header_end + 1 + size + 1
        if b"\x00" in raw:
            binary += 1
            continue
        scanned += 1
        content = raw.decode("utf-8", "replace")
        hits = [index for index, pattern in enumerate(IDENTITY_PATTERNS, start=1)
                if re.search(pattern, content)]
        if not hits:
            continue
        path = blobs.get(sha, "")
        if path.replace("\\", "/") == FIXTURE_HOLDER:
            fixture_holder.append({
                "blob": sha,
                "path": path,
                "pattern_indexes": ",".join(str(index) for index in hits),
                "why": "test fixture holder: the gate skips its own module for the same reason",
            })
            continue
        violations.append({
            "blob": sha,
            "path": path,
            "pattern_indexes": ",".join(str(index) for index in hits),
        })

    for violation in violations:
        containing = git(
            "log", "--all", "--oneline", f"--find-object={violation['blob']}",
            workdir=workdir,
        ).strip().splitlines()
        violation["commit_count"] = str(len(containing))
        violation["commit_shas"] = " ".join(line.split(" ", 1)[0] for line in containing)

    return {
        "commits_examined": commits,
        "unique_blobs": len(blobs),
        "blobs_sized": len(sizes),
        "text_blobs_scanned": scanned,
        "binary_blobs_skipped": binary,
        "oversize_blobs_skipped": len(skipped),
        "oversize_blobs_detail": skipped[:20],
        "violating_blobs": len(violations),
        "violations": violations,
        "fixture_holder_blobs_exempted": len(fixture_holder),
        "fixture_holder_detail": fixture_holder,
        "pattern_count": len(IDENTITY_PATTERNS),
        "verdict": (
            "history is CLEAN under the same patterns the working-tree gate uses"
            if not violations else
            "history CONTAINS identity strings; a clean working tree is not sufficient, "
            "and rewriting history needs explicit authorisation"
        ),
        "boundary": (
            "Read-only measurement. Nothing is rewritten and no ref is moved. This module "
            "makes no judgement about whether a rewrite is warranted."
        ),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", type=Path, default=ROOT)
    parser.add_argument("--out", type=Path, default=None,
                        help="optional JSON report path; nothing is written without it")
    args = parser.parse_args(argv)

    report = scan(args.workdir)
    if args.out is not None:
        if args.out.exists():
            raise ValueError(f"refusing to overwrite an existing report: {args.out}")
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(
            json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
            encoding="utf-8", newline="\n",
        )
    print(json.dumps({
        "commits_examined": report["commits_examined"],
        "unique_blobs": report["unique_blobs"],
        "text_blobs_scanned": report["text_blobs_scanned"],
        "binary_blobs_skipped": report["binary_blobs_skipped"],
        "oversize_blobs_skipped": report["oversize_blobs_skipped"],
        "violating_blobs": report["violating_blobs"],
        "verdict": report["verdict"],
    }, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
