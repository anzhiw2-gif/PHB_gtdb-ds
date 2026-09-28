#!/usr/bin/env python3
"""Prove an adapter revision is additive by diffing it against the shipped output.

Both files come from the SAME twelve frozen inputs; the only difference is the
adapter revision. So the expected result is: identical accessions in identical
order, identical values in every shared column, and new columns present only in
the revision. Anything else means the revision changed behaviour rather than
adding a column, and the caller needs to know that before deciding to rebuild.

Read-only. Streams both files row by row so the 75 MB inputs do not have to be
held in memory twice.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
from pathlib import Path


def read_header(path: Path) -> list[str]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        return next(reader)


def rows(path: Path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            yield row


def compare(shipped: Path, revision: Path, *, max_examples: int = 10) -> dict:
    for path in (shipped, revision):
        if not path.is_file() or path.is_symlink():
            raise ValueError(f"input is not a regular file: {path}")

    shipped_columns = read_header(shipped)
    revision_columns = read_header(revision)
    shared = [column for column in shipped_columns if column in revision_columns]
    added = [column for column in revision_columns if column not in shipped_columns]
    removed = [column for column in shipped_columns if column not in revision_columns]
    if "accession" not in shared:
        raise ValueError("both files must carry an accession column to be aligned")

    mismatched_values: list[str] = []
    count = 0
    extra_shipped = 0
    extra_revision = 0
    # zip_longest, not zip: with plain zip the longer iterator is pulled one item
    # past the end and that item is DISCARDED, so a row-count difference between
    # the two files would go unnoticed - which is exactly the defect the first
    # version of this module had.
    sentinel = object()
    for shipped_row, revision_row in itertools.zip_longest(rows(shipped), rows(revision),
                                                           fillvalue=sentinel):
        if shipped_row is sentinel:
            extra_revision += 1
            continue
        if revision_row is sentinel:
            extra_shipped += 1
            continue
        count += 1
        if shipped_row["accession"] != revision_row["accession"]:
            mismatched_values.append(
                f"row {count}: accession {shipped_row['accession']!r} vs "
                f"{revision_row['accession']!r}"
            )
        for column in shared:
            if column == "accession":
                continue
            if shipped_row.get(column) != revision_row.get(column):
                if len(mismatched_values) < max_examples:
                    mismatched_values.append(
                        f"row {count} {shipped_row['accession']} {column}: "
                        f"{shipped_row.get(column)!r} vs {revision_row.get(column)!r}"
                    )

    populated = {}
    for column in added:
        with revision.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            total = 0
            non_empty = 0
            for row in reader:
                total += 1
                if (row.get(column) or "").strip():
                    non_empty += 1
        populated[column] = {"rows": total, "non_empty": non_empty}

    same_rows = extra_shipped == 0 and extra_revision == 0
    verdict = (
        "ADDITIVE: identical rows and values; the revision only adds column(s)"
        if same_rows and not mismatched_values and not removed
        else "NOT ADDITIVE: the revision changes existing output, see the examples"
    )
    return {
        "shipped_rows": count + extra_shipped,
        "revision_rows": count + extra_revision,
        "rows_compared": count,
        "shipped_columns": len(shipped_columns),
        "revision_columns": len(revision_columns),
        "shared_columns": len(shared),
        "added_columns": added,
        "removed_columns": removed,
        "added_column_fill": populated,
        "mismatched_value_count": len(mismatched_values),
        "mismatch_examples": mismatched_values[:max_examples],
        "verdict": verdict,
        "boundary": (
            "Read-only comparison of two adapter outputs over the same frozen inputs. "
            "It says nothing about whether rebuilding the shipped artifact is warranted."
        ),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shipped", required=True, type=Path)
    parser.add_argument("--revision", required=True, type=Path)
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args(argv)

    report = compare(args.shipped, args.revision)
    if args.out is not None:
        if args.out.exists():
            raise ValueError(f"refusing to overwrite an existing report: {args.out}")
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(
            json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
            encoding="utf-8", newline="\n",
        )
    print(json.dumps({
        "rows_compared": report["rows_compared"],
        "shipped_columns": report["shipped_columns"],
        "revision_columns": report["revision_columns"],
        "added_columns": report["added_columns"],
        "removed_columns": report["removed_columns"],
        "added_column_fill": report["added_column_fill"],
        "mismatched_value_count": report["mismatched_value_count"],
        "verdict": report["verdict"],
    }, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
