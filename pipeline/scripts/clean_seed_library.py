#!/usr/bin/env python3
"""Create an immutable cleaned seed-library input bundle for one dated run.

The tool never edits its inputs.  It removes the fixed non-depolymerase
accessions from every family FASTA and manifest record, while copying the
verified negative control FASTA into the resulting explicit output directory.
"""
import argparse
import csv
import hashlib
import json
import re
import shutil
import tempfile
from collections import Counter
from pathlib import Path


EXCLUDED_ACCESSIONS = frozenset({
    "P29147", "Q02337", "Q02338", "Q80XN0", "P86198", "Q5ZJZ5", "D4A1J4",
    "Q561X9", "Q8JZV9", "Q9BUT1", "C1C4R8", "Q3KPT7", "Q3T046",
    "Q79F77", "Q1EPR4", "Q1EPR5", "P07061", "P07062",
})
REQUIRED_FAMILIES = (
    "ePhaZ", "iPhaZ", "OH", "BdhA", "ArchPhaZ_patatin", "PhaJ", "PhaC", "phasin",
)
DERIVED_FAMILIES = ("ArchPhaZ_hydrolase",)
# ArchPhaZ_patatin contains a small derived expansion whose accession rows
# predate the v2 manifest.  Keep those records, but report them explicitly.
EXPANDED_FAMILIES = ("ArchPhaZ_patatin",)
SAFE_RUN_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


class CleaningError(ValueError):
    """Raised when a cleanup input contract is incomplete or unsafe."""


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def accession_from_header(header):
    accession = header[1:].split("|", 1)[0].strip()
    if not accession:
        raise CleaningError("FASTA header has no accession: " + header.rstrip())
    return accession


def read_fasta(path):
    records = []
    header = None
    sequence = []
    for raw_line in Path(path).read_text(encoding="utf-8").splitlines():
        if raw_line.startswith(">"):
            if header is not None:
                if not sequence:
                    raise CleaningError("FASTA record has an empty sequence: " + header)
                records.append((accession_from_header(header), header, sequence))
            header = raw_line
            sequence = []
        elif not raw_line.strip():
            continue
        elif header is None:
            raise CleaningError("FASTA sequence occurs before a header in " + str(path))
        else:
            sequence.append(raw_line.strip())
    if header is not None:
        if not sequence:
            raise CleaningError("FASTA record has an empty sequence: " + header)
        records.append((accession_from_header(header), header, sequence))
    if not records:
        raise CleaningError("FASTA contains no records: " + str(path))
    accessions = [record[0] for record in records]
    duplicates = sorted({accession for accession in accessions if accessions.count(accession) > 1})
    if duplicates:
        raise CleaningError("duplicate accession(s) in FASTA: " + ",".join(duplicates))
    return records


def write_fasta(path, records):
    with Path(path).open("w", encoding="utf-8", newline="\n") as handle:
        for _, header, sequence in records:
            handle.write(header + "\n")
            handle.write("\n".join(sequence) + "\n")


def file_metadata(path):
    path = Path(path)
    return {"path": path.name, "bytes": path.stat().st_size, "sha256": sha256_file(path)}


def _assert_output_is_safe(seed_dir, output_dir):
    if output_dir.exists():
        raise CleaningError("output directory already exists: " + str(output_dir))
    try:
        output_dir.relative_to(seed_dir)
    except ValueError:
        pass
    else:
        raise CleaningError("output directory must not be inside seed directory")

    # Every cleanup is an auditable run input.  Requiring this exact shape
    # prevents accidental writes to a shared/legacy directory or traversal
    # paths such as ``runs/../...``.
    if output_dir.name != "seeds_clean" or output_dir.parent.name != "inputs":
        raise CleaningError("output directory must be runs/<run_id>/inputs/seeds_clean")
    runs_dir = output_dir.parent.parent.parent
    run_dir = output_dir.parent.parent
    if runs_dir.name != "runs" or not SAFE_RUN_ID.fullmatch(run_dir.name):
        raise CleaningError("output directory must be runs/<safe_run_id>/inputs/seeds_clean")

    current = output_dir
    while True:
        if current.is_symlink():
            raise CleaningError("output path cannot contain symlink components: " + str(current))
        parent = current.parent
        if parent == current:
            break
        current = parent


def _assert_regular_nonempty(path, label):
    if path.is_symlink() or not path.is_file():
        raise CleaningError(f"{label} must be a regular file: {path}")
    if path.stat().st_size <= 0:
        raise CleaningError(f"{label} must be nonempty: {path}")


def verify_sha256_tsv(output_dir):
    """Return checksum/provenance errors for an output bundle.

    ``sha256.tsv`` intentionally lists every file except itself, avoiding a
    recursive self-hash.  The verifier also checks that the table has no
    duplicate or unsafe paths and covers the directory exactly.
    """
    output_dir = Path(output_dir)
    checksum_path = output_dir / "sha256.tsv"
    errors = []
    if not checksum_path.is_file():
        return ["missing sha256.tsv"]
    try:
        with checksum_path.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
    except (OSError, csv.Error) as error:
        return [f"cannot read sha256.tsv: {error}"]
    if not rows or set(rows[0]) != {"path", "bytes", "sha256"}:
        errors.append("sha256.tsv requires path, bytes, sha256 columns")
    listed = set()
    for row in rows:
        name = (row.get("path") or "").strip()
        if not name or Path(name).name != name or name in listed or name == "sha256.tsv":
            errors.append(f"invalid or duplicate checksum path: {name}")
            continue
        listed.add(name)
        target = output_dir / name
        if not target.is_file():
            errors.append(f"checksum target missing: {name}")
            continue
        try:
            expected_bytes = int(row["bytes"])
        except (KeyError, TypeError, ValueError):
            errors.append(f"invalid byte count for {name}")
            continue
        expected_hash = (row.get("sha256") or "").strip()
        if not re.fullmatch(r"[0-9a-f]{64}", expected_hash):
            errors.append(f"invalid sha256 for {name}")
            continue
        actual_hash = sha256_file(target)
        if target.stat().st_size != expected_bytes or actual_hash != expected_hash:
            errors.append(f"checksum mismatch: {name}")
    actual = {path.name for path in output_dir.iterdir() if path.is_file()} - {"sha256.tsv"}
    errors.extend(f"unlisted output: {name}" for name in sorted(actual - listed))
    errors.extend(f"stale checksum entry: {name}" for name in sorted(listed - actual))
    provenance_path = output_dir / "cleaning_manifest.json"
    if provenance_path.is_file():
        try:
            provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
            if provenance.get("checksum") != {
                "file": "sha256.tsv", "scope": "all output files except sha256.tsv"
            }:
                errors.append("cleaning_manifest checksum scope is invalid")
            for name, metadata in provenance.get("outputs", {}).items():
                row_path = output_dir / name
                if name not in listed or not row_path.is_file():
                    errors.append(f"provenance output missing from checksum table: {name}")
                elif metadata.get("sha256") != next(
                    (row.get("sha256") for row in rows if row.get("path") == name), None
                ):
                    errors.append(f"provenance/checksum mismatch: {name}")
        except (OSError, ValueError, TypeError) as error:
            errors.append(f"invalid cleaning_manifest.json: {error}")
    else:
        errors.append("missing cleaning_manifest.json")
    return errors


def clean_library(seed_dir, manifest_path, negative_fasta, output_dir, manifest_output_name=None):
    """Clean a library into a new output directory and return its provenance."""
    seed_dir = Path(seed_dir)
    manifest_path = Path(manifest_path)
    negative_fasta = Path(negative_fasta)
    output_dir = Path(output_dir)
    for input_path in (seed_dir, manifest_path, negative_fasta, output_dir):
        current = input_path
        while True:
            if current.is_symlink():
                raise CleaningError("input path cannot contain symlink components: " + str(current))
            parent = current.parent
            if parent == current:
                break
            current = parent
    seed_dir = seed_dir.resolve()
    manifest_path = manifest_path.resolve()
    negative_fasta = negative_fasta.resolve()
    output_dir = output_dir.resolve()
    if not seed_dir.is_dir():
        raise CleaningError("seed directory is missing: " + str(seed_dir))
    for required in (manifest_path, negative_fasta):
        if not required.is_file():
            raise CleaningError("required input is missing: " + str(required))
    _assert_output_is_safe(seed_dir, output_dir)
    fasta_paths = []
    expected_families = (*REQUIRED_FAMILIES, *DERIVED_FAMILIES)
    for family in expected_families:
        path = seed_dir / f"{family}.faa"
        _assert_regular_nonempty(path, f"seed FASTA {family}")
        fasta_paths.append(path)
    extras = sorted(path.name for path in seed_dir.glob("*.faa") if path.name not in {p.name for p in fasta_paths})
    if extras:
        raise CleaningError("unexpected family FASTA(s): " + ",".join(extras))

    family_records = {path.stem: read_fasta(path) for path in fasta_paths}
    retained_by_family = {
        family: [record for record in records if record[0] not in EXCLUDED_ACCESSIONS]
        for family, records in family_records.items()
    }
    empty_cleaned = sorted(family for family, records in retained_by_family.items() if not records)
    if empty_cleaned:
        raise CleaningError("cleaned family FASTA would be empty: " + ",".join(empty_cleaned))
    training_accession_counts = Counter(record[0] for records in family_records.values() for record in records)
    training_accessions = set(training_accession_counts)
    missing_training = sorted(EXCLUDED_ACCESSIONS - training_accessions)
    if missing_training:
        raise CleaningError("excluded accessions missing from seed FASTA: " + ",".join(missing_training))
    duplicate_excluded = sorted(accession for accession in EXCLUDED_ACCESSIONS if training_accession_counts[accession] != 1)
    if duplicate_excluded:
        raise CleaningError("excluded accession must occur once in training FASTA: " + ",".join(duplicate_excluded))

    with manifest_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not {"accession", "family"}.issubset(reader.fieldnames):
            raise CleaningError("manifest requires accession and family columns")
        manifest_rows = list(reader)
        fieldnames = reader.fieldnames
    if any(not isinstance(row.get("accession"), str) or not row["accession"].strip() for row in manifest_rows):
        raise CleaningError("manifest contains an empty accession")
    if any(not isinstance(row.get("family"), str) or not row["family"].strip() for row in manifest_rows):
        raise CleaningError("manifest contains an empty family")
    manifest_accessions = {row["accession"].strip() for row in manifest_rows}
    if len(manifest_accessions) != len(manifest_rows):
        raise CleaningError("duplicate accession(s) in manifest")
    manifest_by_family = {}
    manifest_by_accession = {}
    for row in manifest_rows:
        family = row["family"].strip()
        if not family:
            raise CleaningError("manifest contains an empty family")
        manifest_by_family.setdefault(family, set()).add(row["accession"].strip())
        manifest_by_accession[row["accession"].strip()] = row
    if set(manifest_by_family) != set(REQUIRED_FAMILIES):
        raise CleaningError("manifest families must exactly match required family set")
    fasta_families_by_accession = {}
    for family, records in family_records.items():
        for accession, header, _ in records:
            fasta_families_by_accession.setdefault(accession, set()).add(family)
    curated_accessions = {
        record[0]
        for family in REQUIRED_FAMILIES
        for record in family_records[family]
        if family not in EXPANDED_FAMILIES
    }
    missing_fasta_manifest = sorted(curated_accessions - manifest_accessions)
    if missing_fasta_manifest:
        raise CleaningError("FASTA accessions missing from manifest: " + ",".join(missing_fasta_manifest[:10]))
    unmanifested_accessions = sorted(
        training_accessions - manifest_accessions
    )
    for accession, families in fasta_families_by_accession.items():
        if families.intersection((*DERIVED_FAMILIES, *EXPANDED_FAMILIES)):
            continue
        row = manifest_by_accession[accession]
        declared = {row["family"].strip()}
        extra = row.get("families", "") or ""
        declared.update(part.strip() for part in re.split(r"[,;]", extra) if part.strip())
        if not families.intersection(declared):
            raise CleaningError("manifest/FASTA family mismatch for accession: " + accession)
    duplicate_excluded_manifest = sorted(
        accession for accession in EXCLUDED_ACCESSIONS if sum(row["accession"].strip() == accession for row in manifest_rows) != 1
    )
    if duplicate_excluded_manifest:
        raise CleaningError("excluded accession must occur once in manifest: " + ",".join(duplicate_excluded_manifest))
    missing_manifest = sorted(EXCLUDED_ACCESSIONS - manifest_accessions)
    if missing_manifest:
        raise CleaningError("excluded accessions missing from manifest: " + ",".join(missing_manifest))
    missing_sequence = sorted(manifest_accessions - training_accessions)
    if missing_sequence:
        raise CleaningError("manifest accessions missing from seed FASTA: " + ",".join(missing_sequence))

    _assert_regular_nonempty(negative_fasta, "negative FASTA")
    negative_records = read_fasta(negative_fasta)
    negative_accessions = {record[0] for record in negative_records}
    if len(negative_records) != len(EXCLUDED_ACCESSIONS) or negative_accessions != EXCLUDED_ACCESSIONS:
        missing_negative = sorted(EXCLUDED_ACCESSIONS - negative_accessions)
        extra_negative = sorted(negative_accessions - EXCLUDED_ACCESSIONS)
        details = []
        if missing_negative:
            details.append("missing=" + ",".join(missing_negative))
        if extra_negative:
            details.append("extra=" + ",".join(extra_negative))
        raise CleaningError("negative FASTA must contain exactly the 18 exclusions (" + "; ".join(details) + ")")

    manifest_output_name = manifest_output_name or manifest_path.name
    manifest_name_path = Path(manifest_output_name)
    if (not manifest_output_name or manifest_name_path.is_absolute()
            or manifest_name_path.name != manifest_output_name
            or manifest_output_name in {".", ".."}):
        raise CleaningError("manifest output name must be a simple basename")
    if manifest_output_name in {"sha256.tsv", "cleaning_manifest.json", "excluded_accessions.tsv", "negative.faa"}:
        raise CleaningError("manifest output name is reserved: " + manifest_output_name)

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    temporary_dir = Path(tempfile.mkdtemp(prefix=output_dir.name + ".tmp-", dir=output_dir.parent))
    try:
        excluded = []
        for family_name, records in family_records.items():
            retained = retained_by_family[family_name]
            write_fasta(temporary_dir / f"{family_name}.faa", retained)
            for accession, _, _ in records:
                if accession in EXCLUDED_ACCESSIONS:
                    excluded.append((accession, family_name))
        retained_rows = [row for row in manifest_rows if row["accession"].strip() not in EXCLUDED_ACCESSIONS]
        with (temporary_dir / manifest_output_name).open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(retained_rows)
        shutil.copyfile(negative_fasta, temporary_dir / "negative.faa")
        with (temporary_dir / "excluded_accessions.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(("accession", "family_fasta", "reason"))
            for accession, family_name in sorted(excluded):
                writer.writerow((accession, family_name, "known_non_depolymerase"))

        outputs = {path.name: file_metadata(path) for path in sorted(temporary_dir.iterdir()) if path.is_file()}
        provenance = {
            "schema_version": 2,
            "excluded_accessions_count": len(EXCLUDED_ACCESSIONS),
            "excluded_accessions": sorted(EXCLUDED_ACCESSIONS),
            "inputs": {
                "seed_dir": str(seed_dir),
                "seed_fastas": {path.name: file_metadata(path) for path in fasta_paths},
                "manifest": {"path": str(manifest_path), "bytes": manifest_path.stat().st_size, "sha256": sha256_file(manifest_path)},
                "negative_fasta": {"path": str(negative_fasta), "bytes": negative_fasta.stat().st_size, "sha256": sha256_file(negative_fasta)},
            },
            "negative_fasta": {"output": "negative.faa", "accession_count": len(negative_accessions), "excluded_intersection_count": len(EXCLUDED_ACCESSIONS & negative_accessions)},
            "manifest_output": manifest_output_name,
            "unmanifested_fasta_accessions": unmanifested_accessions,
            "run_id": output_dir.parent.parent.name,
            "checksum": {"file": "sha256.tsv", "scope": "all output files except sha256.tsv"},
            "outputs": outputs,
        }
        (temporary_dir / "cleaning_manifest.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        with (temporary_dir / "sha256.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(("path", "bytes", "sha256"))
            for path in sorted(item for item in temporary_dir.iterdir() if item.is_file()):
                if path.name == "sha256.tsv":
                    continue
                metadata = file_metadata(path)
                writer.writerow((metadata["path"], metadata["bytes"], metadata["sha256"]))
        temporary_dir.rename(output_dir)
    except Exception:
        shutil.rmtree(temporary_dir, ignore_errors=True)
        raise
    return json.loads((output_dir / "cleaning_manifest.json").read_text(encoding="utf-8"))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed-dir", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--negative-fasta", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--manifest-output-name", default=None,
                        help="output manifest basename (default: input manifest basename)")
    args = parser.parse_args(argv)
    try:
        result = clean_library(
            args.seed_dir, args.manifest, args.negative_fasta, args.output_dir,
            manifest_output_name=args.manifest_output_name,
        )
    except CleaningError as error:
        parser.error(str(error))
    print("cleaned {} excluded accessions -> {}".format(result["excluded_accessions_count"], args.output_dir))


if __name__ == "__main__":
    main()
