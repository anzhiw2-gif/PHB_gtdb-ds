#!/usr/bin/env python3
"""Task 9 — one HMMER command shape pinned to the *full* database size ``Z``.

Why this module exists
----------------------
HMMER E-values are ``E = Z * P``: the E-value scales linearly with the number
of target sequences ``Z`` that the search is told it looked at.  A sharded
library search that omits ``-Z`` therefore silently uses HMMER's default,
``Z = <number of sequences in this shard>``.  Each shard then reports E-values
on its *own* scale, and those E-values can never be compared with each other or
reported as a single full-library scale.

The frozen GTDB scan 13 (``runs/20260901_formal_frozen_scan_13``) was executed
exactly that way (see the Task 9 audit in
``docs/T141_20260928_phaded_evidence_model_redesign_status.md``), which is why
the requirement in
``docs/superpowers/specs/2026-09-28-phaded-evidence-model-redesign.md``
("HMMER 分片尺度") is: 正式入口必须记录完整目标序列数 ``Z``; every shard must be
searched with the *same* full-database ``-Z``; the manifest must record
``database_size_Z``, the per-shard sequence counts, and the command.

This module is the single, dependency-free source of truth for that command
shape.  ``pipeline/scripts/06_screen.sh`` builds the same flag set in Bash, and
``pipeline/tests/test_hmmer_database_size.py`` asserts both stay aligned.

Contract
--------
* ``-Z <total_targets>`` is **always** present and equals ``str(total_targets)``.
* ``total_targets`` must be a positive ``int`` (``bool`` is rejected); anything
  else raises ``ValueError`` — there is no default and no silent fallback.
* ``--domZ`` is **never** added: domain E-values are not rescaled anywhere in
  the pipeline, so adding it would create a second, undocumented scale.  The
  plan requires a separate verified decision before that changes.
* The result is a ``list[str]`` of literal argv elements, safe to hand to
  ``subprocess`` without a shell: paths are never interpolated, quoted, or
  word-split.

The rescaling convention itself lives in
``pipeline/scripts/parse_phaded_cys_targeted_recall.py::rescale_evalue`` and is
reused rather than reinvented here.
"""

from __future__ import annotations

import shlex

PROGRAM = "hmmsearch"
DATABASE_SIZE_FLAG = "-Z"
# Recorded verbatim in the screen manifest so a later reader can see the exact
# command shape that produced the E-values.
COMMAND_TEMPLATE = (
    'hmmsearch --tblout "$out" --domtblout "${out%.tbl}.dom" '
    '-Z "$DATABASE_SIZE_Z" -E "$EVAL" --cpu "$HMM_CPU" "$hmm" "$shard"'
)


def _require_non_empty_text(value, field):
    """Return ``value`` as text, raising ``ValueError`` when it is unusable."""
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{field} must be a non-empty string, got {value!r}")
    return value


def validate_total_targets(total_targets):
    """Return the full-database target count as a positive ``int``.

    ``bool`` is rejected explicitly: ``True`` is an ``int`` in Python but is
    never a meaningful sequence count.  Missing evidence must raise here rather
    than fall back to a shard-local default.
    """
    if isinstance(total_targets, bool) or not isinstance(total_targets, int):
        raise ValueError(
            f"total_targets must be a positive integer, got {total_targets!r}"
        )
    if total_targets < 1:
        raise ValueError(
            f"total_targets must be a positive integer, got {total_targets!r}"
        )
    return total_targets


def validate_cpu(cpu):
    """Return ``cpu`` as a positive ``int`` (HMMER's ``--cpu``)."""
    if isinstance(cpu, bool) or not isinstance(cpu, int) or cpu < 1:
        raise ValueError(f"cpu must be a positive integer, got {cpu!r}")
    return cpu


def build_hmmsearch_command(
    hmm,
    shard,
    *,
    total_targets,
    evalue="1e-5",
    tblout=None,
    domtblout=None,
    cpu=1,
):
    """Build the full-database-scaled ``hmmsearch`` argv for one shard.

    Parameters
    ----------
    hmm:
        Path to the profile HMM.
    shard:
        Path to the shard FASTA that is being searched.
    total_targets:
        The **complete** library target sequence count ``Z`` (every shard gets
        this same value).  Required; must be a positive ``int``.
    evalue:
        Reporting threshold passed as ``-E``.
    tblout, domtblout:
        Optional ``--tblout`` / ``--domtblout`` paths.
    cpu:
        Positive ``int`` passed as ``--cpu``.

    Returns
    -------
    list[str]
        A literal argv list: ``subprocess.run(command, ...)`` needs no shell.
        ``--domZ`` is never included.
    """
    hmm = _require_non_empty_text(hmm, "hmm")
    shard = _require_non_empty_text(shard, "shard")
    evalue = _require_non_empty_text(evalue, "evalue")
    total_targets = validate_total_targets(total_targets)
    cpu = validate_cpu(cpu)

    argv = [PROGRAM]
    if tblout is not None:
        argv += ["--tblout", _require_non_empty_text(tblout, "tblout")]
    if domtblout is not None:
        argv += ["--domtblout", _require_non_empty_text(domtblout, "domtblout")]
    # The whole point of Task 9: one full-library Z for every shard.  Never
    # --domZ: no downstream consumer rescales domain E-values.
    argv += [DATABASE_SIZE_FLAG, str(total_targets)]
    argv += ["-E", evalue, "--cpu", str(cpu), hmm, shard]
    return argv


def command_template() -> str:
    """Return the manifest-facing template matching ``build_hmmsearch_command``."""
    return COMMAND_TEMPLATE


def format_command_for_log(argv) -> str:
    """Render argv for a human-readable log line only (never for execution)."""
    return " ".join(shlex.quote(str(item)) for item in argv)
