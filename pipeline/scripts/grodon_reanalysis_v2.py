#!/usr/bin/env python3
"""gRodon re-analysis v2: corrected labels, units, overlap, and dependence.

PLACEMENT NOTE (read before running this file)
----------------------------------------------
This module is the tested implementation of plan **Task 12** of
``docs/superpowers/plans/2026-09-28-phaded-evidence-model-redesign.md``.

The plan's file map puts the real code in
``deploy/20260928_phaded_grodon_reanalysis_v2_01/scripts/`` once execution is
authorized.  Creating that dated deploy directory requires separate
authorization, so the failure-test-first implementation lives here, in
``pipeline/scripts/grodon_reanalysis_v2.py``.  When execution is authorized the
dated deploy copies this module (byte-identical, hash recorded in the deploy
manifest); it must not be re-derived by hand.

This module **supersedes** the two legacy scripts of the frozen deploy
``deploy/20260920_phaded_grodon_growth_01/``::

    scripts/build_grodon_manifest.py   -> build_manifest()
    scripts/grodon_group_stats.py      -> genus_level_* functions

The frozen deploy is historical evidence and is never modified.  Only the
scientific corrections are ported: same-genus 1:1 balancing with a
deterministic seed, existing-FASTA filtering, genus-level deltas, BH
adjustment, and the between-group comparison all keep their proven logic.

What v2 corrects (design document, "下游统计"):

1. **Labels and overlap.**  Group names are "candidate gene carrier"
   (``candidate_gene_carrier``) and "candidate not detected under the defined
   search/quality conditions"
   (``candidate_not_detected_under_defined_search``).  No phenotype word is
   emitted.  The final positive set and the control set are asserted to be
   disjoint and the builder refuses to construct a manifest that would overlap
   them.
2. **Units.**  gRodon returns the predicted minimum doubling time in hours
   (``res$d`` -> ``doubling_time_h``).  ``growth_rate_per_h()`` is the single
   tested conversion (``ln(2) / d``); ``resolve_growth_rate()`` decides which
   field is authoritative so the conversion is applied exactly once and never
   twice.
3. **Dependence.**  The **genus-level mean difference is the primary
   estimand**.  Bootstrap and permutation resample *genera* (clusters), never
   individual genomes.  Median, sign test, per-genome and
   intracellular/extracellular comparisons are secondary and labelled as such.
   "No difference" becomes "no difference detected under the current design"
   unless an equivalence margin is preregistered.
4. **Reuse audit.**  The 713 reused predictions of the frozen run enter the
   primary analysis only when tool, tool version, model settings, temperature
   handling, CDS preparation, ribosomal-marker procedure and output formula are
   all comparable; the rest stay in a sensitivity table with a named reason.
5. **Selection coverage.**  Included / no-same-genus-control / missing-FASTA /
   failed-prediction genomes are stratified by taxonomy, candidate family,
   group and available quality metrics.
6. **The 66-genome gap.**  The 4,507 eligible positives to 4,441 manifest
   positives difference is explained by named reason buckets whose counts sum
   to exactly 66 and are computed from **accession-level sets**, never by
   subtracting aggregates.

Boundary: this module is offline and tool-free.  It never imports or executes
gRodon/R, never touches the network, and never writes into ``runs/``,
``results/`` or any ``deploy/`` directory.  Growth rates are codon-usage
*predictions* of maximum growth potential, and every label is candidate-only:
carrying a candidate gene is not a validated PHB/PHA degradation phenotype.
"""

from __future__ import annotations

import math
import sys
from collections import Counter, defaultdict
from typing import Any, Callable, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

__all__ = [
    "CANDIDATE_CARRIER",
    "CANDIDATE_NOT_DETECTED",
    "CONTROL_GROUP",
    "EQUIVALENCE_MARGIN_ABSENT",
    "COMPARABILITY_FIELDS",
    "FORBIDDEN_LEGACY_LABELS",
    "GROWTH_RATE_UNIT",
    "GROWTH_RATE_FORMULA",
    "LEGACY_PHENOTYPE_LABELS",
    "MANIFEST_COLUMNS",
    "PRIMARY_ESTIMAND",
    "RESAMPLING_UNIT",
    "SECONDARY_ANALYSES",
    "SELECTION_STATUSES",
    "assert_sets_are_disjoint",
    "audit_reused_predictions",
    "audit_reused_predictions_table",
    "carrier_genomes",
    "control_genomes",
    "cluster_bootstrap_mean_difference",
    "equivalence_statement",
    "explain_manifest_difference",
    "growth_rate_formula",
    "growth_rate_per_h",
    "genus_level_estimates",
    "genus_level_deltas",
    "genus_level_permutation_test",
    "intracellular_vs_extracellular",
    "resolve_growth_rate",
    "resolve_growth_rate_series",
    "resolve_growth_rate_source",
    "sets_are_disjoint",
    "summarize_selection_coverage",
]

# --------------------------------------------------------------------- labels

#: genome carries at least one candidate gene under the defined search
CANDIDATE_CARRIER = "candidate_gene_carrier"
#: no candidate was detected under the defined search/quality conditions
CANDIDATE_NOT_DETECTED = "candidate_not_detected_under_defined_search"
#: value of the ``group`` column for a control row
CONTROL_GROUP = "control"

DETECTION_STATUSES = (CANDIDATE_CARRIER, CANDIDATE_NOT_DETECTED)

#: historical phenotype vocabulary that must never be emitted again.  The
#: frozen 2026-09-20 deploy used these as group names; they assert a phenotype
#: that a candidate classification does not establish.
LEGACY_PHENOTYPE_LABELS = ("degrader", "non-degrader")
FORBIDDEN_LEGACY_LABELS = (
    "degrader",
    "non-degrader",
    "non degrader",
    "non_degrader",
)

#: exact column order of a v2 manifest row
MANIFEST_COLUMNS = (
    "genome_id",
    "candidate_detection_status",
    "candidate_count",
    "candidate_family",
    "group",
    "pair_id",
    "control_basis",
    "domain",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
    "genome_path",
)

# ---------------------------------------------------------------------- units

#: unit of every growth-rate number this module produces or consumes
GROWTH_RATE_UNIT = "1/h"

#: doubling time field names accepted from a gRodon prediction table.  ``d`` is
#: the field name of gRodon's own result object (``res$d``, minimum doubling
#: time in hours); the frozen deploy copied it into ``doubling_time_h``.
DOUBLING_TIME_FIELDS = ("doubling_time_h", "d", "doubling_time", "min_doubling_time_h")

#: growth-rate field names accepted from a table that already carries a rate
GROWTH_RATE_FIELDS = (
    "growth_rate_per_h",
    "growth_rate",
    "growth_rate_1_per_h",
    "mu_per_h",
)

#: fields that state the unit of the growth-rate column
GROWTH_RATE_UNIT_FIELDS = ("growth_rate_unit", "rate_unit", "growth_rate_units")

_ACCEPTED_RATE_UNITS = frozenset(
    {
        "1/h",
        "h^-1",
        "h-1",
        "h^{-1}",
        "per_h",
        "perhour",
        "per_hour",
        "hour^-1",
        "hours^-1",
        "/h",
    }
)

#: relative tolerance for accepting a table that carries both a doubling time
#: and a pre-computed rate as the *same* quantity
FORMULA_TOLERANCE = 5e-2


def _canonical_unit(unit: str) -> str:
    return unit.strip().lower().replace(" ", "").replace("**", "^")


def growth_rate_per_h(doubling_time_h: float) -> float:
    """Convert a gRodon *minimum doubling time in hours* to a growth rate.

    ``mu = ln(2) / d`` where ``d`` is the predicted minimum doubling time in
    hours and ``mu`` is the predicted maximum growth rate in ``1/h``.  This is
    the single conversion in this module: call it once per genome, on a
    doubling time, and never on a value that is already a rate.
    """
    value = float(doubling_time_h)
    if not math.isfinite(value) or value <= 0:
        raise ValueError("doubling time must be positive")
    return math.log(2.0) / value


def growth_rate_formula(doubling_time_field: str = "doubling_time_h") -> str:
    """Canonical string describing how a growth-rate column was produced.

    Records that carry a growth rate must state the formula that produced it,
    so a reused table built with a different convention (``1/d``, minutes
    instead of hours, an already-converted rate converted again) is detectable.
    """
    return f"growth_rate_per_h=ln2/{doubling_time_field}"


GROWTH_RATE_FORMULA = growth_rate_formula()


def _clean_field(value: Any) -> str | None:
    if value is None:
        return None
    if isinstance(value, float) and math.isnan(value):
        return None
    text = str(value).strip()
    return text or None


def _numeric_field(
    record: Mapping[str, Any],
    names: Sequence[str],
    *,
    label: str,
) -> tuple[str, float] | None:
    """First present, numeric field among ``names``; sign is checked by the caller."""
    for name in names:
        if name not in record:
            continue
        raw = record[name]
        if raw is None:
            continue
        if isinstance(raw, float) and not math.isfinite(raw):
            raise ValueError(f"{label} field {name!r} is not finite: {raw!r}")
        text = str(raw).strip()
        if not text:
            continue
        try:
            value = float(text)
        except (TypeError, ValueError):
            raise ValueError(f"{label} field {name!r} is not numeric: {text!r}") from None
        if not math.isfinite(value):
            raise ValueError(f"{label} field {name!r} is not finite: {text!r}")
        return name, value
    return None


def _positive(value: tuple[str, float], *, label: str, message: str) -> tuple[str, float]:
    if value[1] <= 0:
        raise ValueError(f"{message}: {value[0]}={value[1]!r}")
    return value


def _declared_rate_unit(record: Mapping[str, Any]) -> tuple[str, str] | None:
    for name in GROWTH_RATE_UNIT_FIELDS:
        raw = _clean_field(record.get(name))
        if raw is None:
            continue
        if _canonical_unit(raw) not in _ACCEPTED_RATE_UNITS:
            raise ValueError(
                f"{name} must be a per-hour unit, got {raw!r}; expected one of "
                f"{sorted(_ACCEPTED_RATE_UNITS)}"
            )
        return name, raw
    return None


def resolve_growth_rate_source(record: Mapping[str, Any]) -> str:
    """Return which field of ``record`` is authoritative for the growth rate."""
    rate = _numeric_field(record, GROWTH_RATE_FIELDS, label="growth rate")
    doubling = _numeric_field(record, DOUBLING_TIME_FIELDS, label="doubling time")
    if rate is not None:
        return _positive(rate, label="growth rate", message="growth rate must be positive")[0]
    if doubling is not None:
        return _positive(
            doubling, label="doubling time", message="doubling time must be positive"
        )[0]
    raise ValueError(
        "record carries no growth rate: expected one of "
        f"{list(GROWTH_RATE_FIELDS)} or one of {list(DOUBLING_TIME_FIELDS)}"
    )


def resolve_growth_rate(record: Mapping[str, Any]) -> float:
    """Return the growth rate in ``1/h`` for one prediction record, once.

    Resolution rule (documented once, tested in both directions):

    * a growth-rate field is authoritative and is returned **unchanged** -- the
      value must already be a per-hour rate (checked against the declared unit
      when the table states one), so no second conversion is applied;
    * otherwise a doubling-time field is converted exactly once with
      :func:`growth_rate_per_h`;
    * if the record carries *both*, the pre-computed rate wins (prefer the
      recorded value, never re-derive it) and the doubling time is used only as
      a cross-check: a disagreement larger than
      :data:`FORMULA_TOLERANCE` relative error raises ``ValueError`` naming
      ``output_formula``, which is the signature of a rate stored in the wrong
      unit or produced by a different formula.
    """
    _declared_rate_unit(record)
    rate = _numeric_field(record, GROWTH_RATE_FIELDS, label="growth rate")
    doubling = _numeric_field(record, DOUBLING_TIME_FIELDS, label="doubling time")

    if rate is not None:
        _positive(rate, label="growth rate", message="growth rate must be positive")
        if doubling is not None:
            # Convert the doubling time once, purely to cross-check the formula.
            value = _positive(
                doubling, label="doubling time", message="doubling time must be positive"
            )
            derived = growth_rate_per_h(value[1])
            relative = abs(rate[1] - derived) / derived
            if relative > FORMULA_TOLERANCE:
                raise ValueError(
                    "output_formula mismatch: growth rate field "
                    f"{rate[0]!r}={rate[1]!r} is not {growth_rate_formula(doubling[0])} "
                    f"= {derived!r} (relative error {relative:.3g} > {FORMULA_TOLERANCE:g}); "
                    "the record must not be converted a second time and must not "
                    "enter the primary analysis until the formula is verified"
                )
        return rate[1]

    if doubling is not None:
        return growth_rate_per_h(
            _positive(
                doubling, label="doubling time", message="doubling time must be positive"
            )[1]
        )

    raise ValueError(
        "record carries no growth rate: expected one of "
        f"{list(GROWTH_RATE_FIELDS)} or one of {list(DOUBLING_TIME_FIELDS)}"
    )


def resolve_growth_rate_series(
    table: pd.DataFrame,
    canonical_name: str = "growth_rate_per_h",
) -> pd.Series:
    """Vectorised :func:`resolve_growth_rate` over a prediction table.

    Returns a ``1/h`` series indexed like ``table`` (by ``genome_id`` when that
    column exists).  Exactly one conversion is applied per row: a pre-computed
    rate is carried through unchanged and a doubling time is converted once.
    """
    if not isinstance(table, pd.DataFrame):
        raise ValueError("table must be a pandas DataFrame")
    named = table
    if "genome_id" in table.columns:
        named = table.set_index(pd.Index(table["genome_id"].astype(str), name="genome_id"))

    rate_column = next((c for c in GROWTH_RATE_FIELDS if c in named.columns), None)
    doubling_column = next((c for c in DOUBLING_TIME_FIELDS if c in named.columns), None)
    if rate_column is None and doubling_column is None:
        raise ValueError(
            "table carries no growth rate: expected one of "
            f"{list(GROWTH_RATE_FIELDS)} or one of {list(DOUBLING_TIME_FIELDS)}"
        )

    rate = (
        pd.to_numeric(named[rate_column], errors="coerce")
        if rate_column is not None
        else None
    )
    doubling = (
        pd.to_numeric(named[doubling_column], errors="coerce")
        if doubling_column is not None
        else None
    )

    if rate is not None and doubling is not None:
        derived = math.log(2.0) / doubling
        relative = (rate - derived).abs() / derived
        bad = relative > FORMULA_TOLERANCE
        if bool(np.any(bad.fillna(False).to_numpy())):
            offenders = ", ".join(sorted(named.index[bad.fillna(False)]))[:200]
            raise ValueError(
                "output_formula mismatch for genomes "
                f"{offenders}: the pre-computed rate is not "
                f"{growth_rate_formula(doubling_column)}"
            )
        values = rate
    elif rate is not None:
        values = rate
    else:
        if bool(np.any((doubling <= 0).to_numpy())):
            raise ValueError("doubling time must be positive")
        values = math.log(2.0) / doubling

    if bool(values.isna().to_numpy().any()):
        missing = sorted(named.index[values.isna()])
        raise ValueError(f"missing growth rate for genomes {missing[:20]}")
    values = values.astype(float)
    values.name = canonical_name
    return values


# ------------------------------------------------------------------- manifest


def _as_int(value: Any, field: str) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        raise ValueError(f"{field} must be an integer, got {value!r}") from None


def _row_genome_id(key: str, row: Mapping[str, Any]) -> str:
    raw = row.get("genome_id")
    text = "" if raw is None else str(raw).strip()
    return text or str(key)


def _validate_carriers(carriers: Mapping[str, Mapping[str, Any]]) -> dict[str, dict[str, Any]]:
    if carriers is None:
        raise ValueError("carriers mapping must not be None")
    out: dict[str, dict[str, Any]] = {}
    for key, row in carriers.items():
        if row is None:
            raise ValueError(f"carrier row for {key!r} is empty")
        genome_id = _row_genome_id(key, row)
        if genome_id in out:
            raise ValueError(f"duplicate carrier genome_id: {genome_id}")
        out[genome_id] = dict(row)
    return out


def _control_pool(
    tax_by_gid: Mapping[str, Mapping[str, Any]],
    carriers: Mapping[str, Mapping[str, Any]],
    exclusion: set[str],
) -> dict[str, list[str]]:
    """Genomes eligible to be a control, grouped by genus.

    A genome is eligible when it is bacterial, has a usable genus label, is not
    a candidate carrier and is not in the exclusion set.
    """
    pool: dict[str, list[str]] = defaultdict(list)
    for genome_id, tax_row in tax_by_gid.items():
        tax_row = tax_row or {}
        if str(tax_row.get("domain", "")) != "Bacteria":
            continue
        genus = str(tax_row.get("genus") or "").strip()
        if genus in ("", "Unknown", "unclassified", "nan", "None"):
            continue
        if genome_id in carriers or genome_id in exclusion:
            continue
        pool[genus].append(genome_id)
    return {genus: sorted(ids) for genus, ids in pool.items()}


def _tax_row(tax_by_gid: Mapping[str, Mapping[str, Any]], genome_id: str) -> Mapping[str, Any]:
    return tax_by_gid.get(genome_id) or {}


def _manifest_row(
    genome_id: str,
    tax_row: Mapping[str, Any],
    *,
    status: str,
    candidate_count: int,
    candidate_family: str,
    group: str,
    pair_id: str,
    control_basis: str,
    genome_path: str,
) -> dict[str, Any]:
    row = {
        "genome_id": genome_id,
        "candidate_detection_status": status,
        "candidate_count": int(candidate_count),
        "candidate_family": candidate_family,
        "group": group,
        "pair_id": pair_id,
        "control_basis": control_basis,
        "genome_path": genome_path,
    }
    for rank in ("domain", "phylum", "class", "order", "family", "genus", "species"):
        row[rank] = _clean_field(tax_row.get(rank)) or ""
    return {column: row[column] for column in MANIFEST_COLUMNS}


def carrier_genomes(rows: Iterable[Mapping[str, Any]]) -> set[str]:
    """Genomes in the final positive set (candidate carriers) of a manifest."""
    return {
        str(row["genome_id"])
        for row in rows
        if row.get("candidate_detection_status") == CANDIDATE_CARRIER
    }


def control_genomes(rows: Iterable[Mapping[str, Any]]) -> set[str]:
    """Genomes in the control set (no candidate detected) of a manifest."""
    return {
        str(row["genome_id"])
        for row in rows
        if row.get("candidate_detection_status") == CANDIDATE_NOT_DETECTED
    }


def sets_are_disjoint(rows: Iterable[Mapping[str, Any]]) -> bool:
    """True when no genome is both a carrier and a control in ``rows``."""
    rows = list(rows)
    return not (carrier_genomes(rows) & control_genomes(rows))


def assert_sets_are_disjoint(rows: Iterable[Mapping[str, Any]]) -> None:
    """Raise ``ValueError('positive/control overlap ...')`` on any overlap."""
    rows = list(rows)
    positives = carrier_genomes(rows)
    controls = control_genomes(rows)
    overlap = positives & controls
    if overlap:
        raise ValueError(
            "positive/control overlap: "
            f"{len(overlap)} genome(s) appear in both the final candidate-carrier "
            f"set and the control set: {sorted(overlap)[:20]}"
        )


def build_manifest(
    carriers: Mapping[str, Mapping[str, Any]],
    exclusion: set[str],
    tax_by_gid: Mapping[str, Mapping[str, Any]],
    seed: int,
    max_per_genus: int,
    file_exists: Callable[[str], bool] | None = None,
    control_pool: Iterable[str] | None = None,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Build the same-genus matched manifest for the v2 candidate-carrier design.

    Design (ported from the frozen deploy, with the v2 labels):

    * the final positive set is exactly ``carriers`` (genomes carrying a
      candidate gene under the defined search);
    * the control candidates are the bacterial genomes of ``tax_by_gid`` that
      are neither carriers nor in ``exclusion``.  Pass ``control_pool`` to
      declare that candidate set explicitly (the real pipeline materialises it
      once and reuses it); when omitted it is derived from ``tax_by_gid``;
    * ``max_per_genus`` caps the matched pairs per genus (``<= 0`` means no
      limit, as in the frozen deploy) and balance is 1:1 within each genus,
      selected deterministically from ``seed``;
    * only genomes whose FASTA ``file_exists`` are kept, and genera without a
      usable same-genus control are skipped and reported, never filled with a
      genome drawn from outside the genus.

    Overlap rule, asserted not assumed.  The invariant is that the final
    positive set (``carriers``) and the final control set are disjoint: no
    genome may be labelled both a candidate carrier and a control.

    1. Contradictory inputs are rejected before a row is built.  A genome that
       is both a carrier and named in the control-pool ledger, or both a
       carrier and named in the ``exclusion`` ledger, cannot be classified, so
       the call raises ``ValueError`` whose message contains the literal text
       ``positive/control overlap``;
    2. carriers are removed from the control candidate set by construction, so
       a carrier can never be *selected* as a control;
    3. :func:`assert_sets_are_disjoint` re-checks the finished rows and raises
       the same ``ValueError`` if any genome ended up in both sets, which makes
       the disjointness a checked post-condition of every returned manifest
       rather than a property assumed from the code path.

    Consequences worth stating explicitly: the final positive set and the
    control set are provably disjoint for every returned manifest; a genus with
    no usable same-genus control contributes no rows at all (so no genome is
    ever silently re-labelled to fill a pair); and a caller that hands in a
    control pool contaminated with positives gets a hard error instead of a
    control arm built from its own cases.

    Returns ``(rows, stats)``.  Every row is a dict with exactly
    :data:`MANIFEST_COLUMNS` keys, and ``candidate_detection_status`` is either
    :data:`CANDIDATE_CARRIER` or :data:`CANDIDATE_NOT_DETECTED`.
    """
    if file_exists is None:
        def file_exists(_path: str) -> bool:  # pragma: no cover - trivial default
            return True

    cap = _as_int(max_per_genus, "max_per_genus")
    if cap < 0:
        raise ValueError(f"max_per_genus must be >= 0, got {cap}")
    exclusion_set = {str(item) for item in (exclusion or set())}

    carrier_map = _validate_carriers(carriers)
    positive_set = set(carrier_map)

    # Guard 1: the final positive set and the control candidate set must be
    # disjoint.  Both ledgers below name accessions that are, or feed, the
    # control arm; a carrier appearing in either one is a contradiction, and
    # resolving it silently is how a positive genome ends up as its own
    # control.
    overlap = sorted(positive_set & exclusion_set)
    if overlap:
        raise ValueError(
            "positive/control overlap: "
            f"{len(overlap)} genome(s) are declared both as candidate carriers "
            "(final positive set) and in the pool-exclusion ledger used to "
            f"build the control set: {overlap[:20]}"
        )
    if control_pool is not None:
        declared_pool = {str(item) for item in control_pool}
        overlap = sorted(positive_set & declared_pool)
        if overlap:
            raise ValueError(
                "positive/control overlap: "
                f"{len(overlap)} genome(s) are declared both as candidate "
                "carriers (final positive set) and inside the declared control "
                f"pool: {overlap[:20]}"
            )

    pool = _control_pool(tax_by_gid, carrier_map, exclusion_set)
    rng = np.random.default_rng(int(seed))
    by_genus_carriers: dict[str, list[str]] = defaultdict(list)
    for genome_id in sorted(positive_set):
        genus = str(_tax_row(tax_by_gid, genome_id).get("genus") or "").strip()
        if genus in ("", "Unknown"):
            genus = ""
        by_genus_carriers[genus].append(genome_id)

    rows: list[dict[str, Any]] = []
    cache: dict[str, bool] = {}

    def exists(path: str) -> bool:
        if path not in cache:
            cache[path] = bool(file_exists(path))
        return cache[path]

    def genome_path(genome_id: str) -> str:
        return f"{genome_id}_genomic.fna.gz"

    stats: dict[str, Any] = {
        "schema_version": 2,
        "seed": int(seed),
        "max_per_genus": cap,
        "regime": "matched_same_genus_1to1",
        "carrier_genomes_input": len(positive_set),
        "exclusion_genomes": len(exclusion_set),
        "taxonomy_bacterial_genomes": sum(
            1
            for tax_row in tax_by_gid.values()
            if str((tax_row or {}).get("domain", "")) == "Bacteria"
        ),
        "control_pool_genomes": sum(len(ids) for ids in pool.values()),
        "control_pool_genera": len(pool),
    }

    carrier_missing_fasta = 0
    control_missing_fasta = 0
    no_control = 0
    skipped_genera = 0
    matched_pairs = 0

    for genus in sorted(by_genus_carriers):
        carrier_ids = by_genus_carriers[genus]
        candidates = list(pool.get(genus, []))
        if genus == "" or not candidates:
            no_control += len(carrier_ids)
            skipped_genera += 1
            continue
        order = rng.permutation(len(candidates))
        usable_controls = [candidates[i] for i in order if exists(genome_path(candidates[i]))]
        control_missing_fasta += len(candidates) - len(usable_controls)
        usable_carriers = [
            genome_id
            for genome_id in sorted(carrier_ids)
            if exists(str(carrier_map[genome_id].get("genome_path") or genome_path(genome_id)))
        ]
        carrier_missing_fasta += len(carrier_ids) - len(usable_carriers)
        if not usable_carriers or not usable_controls:
            no_control += len(carrier_ids)
            skipped_genera += 1
            continue
        pair_cap = sys.maxsize if cap <= 0 else cap
        n_pairs = min(len(usable_carriers), len(usable_controls), pair_cap)
        matched_pairs += n_pairs
        no_control += len(carrier_ids) - n_pairs
        for index in range(n_pairs):
            genome_id = usable_carriers[index]
            control_id = usable_controls[index]
            pair_id = f"{genus}::{index + 1}"
            info = carrier_map[genome_id]
            rows.append(
                _manifest_row(
                    genome_id,
                    _tax_row(tax_by_gid, genome_id),
                    status=CANDIDATE_CARRIER,
                    candidate_count=_as_int(info.get("candidate_count", 0), "candidate_count"),
                    candidate_family=_clean_field(info.get("candidate_family")) or "",
                    group=_clean_field(info.get("group")) or "",
                    pair_id=pair_id,
                    control_basis="same_genus_no_candidate_detected",
                    genome_path=str(info.get("genome_path") or genome_path(genome_id)),
                )
            )
            rows.append(
                _manifest_row(
                    control_id,
                    _tax_row(tax_by_gid, control_id),
                    status=CANDIDATE_NOT_DETECTED,
                    candidate_count=0,
                    candidate_family="",
                    group=CONTROL_GROUP,
                    pair_id=pair_id,
                    control_basis="same_genus_no_candidate_detected",
                    genome_path=genome_path(control_id),
                )
            )

    stats.update(
        manifest_rows=len(rows),
        manifest_carrier=len(carrier_genomes(rows)),
        manifest_control=len(control_genomes(rows)),
        matched_pairs=matched_pairs,
        manifest_genera=len({row["genus"] for row in rows}),
        genera_with_a_control=len(
            {row["genus"] for row in rows if row["candidate_detection_status"] == CANDIDATE_NOT_DETECTED}
        ),
        control_genomes_no_same_genus_control=no_control,
        carrier_genomes_without_a_same_genus_control=no_control,
        genera_skipped_no_control=skipped_genera,
        carrier_genomes_missing_fasta=carrier_missing_fasta,
        control_genomes_missing_fasta=control_missing_fasta,
        overlap_detected=False,
    )
    assert_sets_are_disjoint(rows)
    return rows, stats


# ----------------------------------------------------------------- statistics

PRIMARY_ESTIMAND = "genus_level_mean_difference"
RESAMPLING_UNIT = "genus"
SECONDARY_ANALYSES = (
    "genus_median_difference",
    "genus_sign_test",
    "genome_level_mann_whitney_u",
    "intracellular_vs_extracellular",
)
EQUIVALENCE_MARGIN_ABSENT = "(not applicable: no equivalence margin preregistered)"


def _require_rng(rng: np.random.Generator | None) -> np.random.Generator:
    if rng is None:
        raise ValueError("an explicit numpy Generator is required for reproducibility")
    if not isinstance(rng, np.random.Generator):
        raise ValueError("rng must be a numpy.random.Generator")
    return rng


def _require_columns(table: pd.DataFrame, needed: Sequence[str]) -> None:
    if not isinstance(table, pd.DataFrame):
        raise ValueError("table must be a pandas DataFrame")
    missing = [column for column in needed if column not in table.columns]
    if missing:
        raise ValueError(f"table is missing required column(s): {missing}")


def _detection_series(table: pd.DataFrame) -> pd.Series:
    return table["candidate_detection_status"].astype(str).str.strip()


def genus_level_deltas(
    table: pd.DataFrame,
    *,
    status_column: str = "candidate_detection_status",
    growth_column: str = "growth_rate_per_h",
    genus_column: str = "genus",
    carrier_value: str = CANDIDATE_CARRIER,
    control_value: str = CANDIDATE_NOT_DETECTED,
    id_column: str = "genome_id",
) -> list[dict[str, Any]]:
    """One balanced genus-level delta per genus that has both arms.

    The analysis unit is the **genus**: carriers and controls are matched 1:1
    within the genus (largest-index arms truncated by genome_id, the same
    deterministic rule as the frozen deploy), and the delta is
    ``mean(carriers) - mean(controls)`` in ``1/h``.  A mask built from a
    filtered sub-frame is accepted: selection is by column values, never by
    index alignment.
    """
    _require_columns(
        table, [status_column, growth_column, genus_column, id_column]
    )
    frame = table.copy()
    frame[growth_column] = pd.to_numeric(frame[growth_column], errors="coerce")
    statuses = frame[status_column].astype(str).str.strip()

    out: list[dict[str, Any]] = []
    for genus in sorted(set(frame[genus_column].astype(str))):
        if genus in ("", "nan", "None"):
            continue
        arm_pos = frame[
            (frame[genus_column].astype(str) == genus) & (statuses == carrier_value)
        ]
        arm_neg = frame[
            (frame[genus_column].astype(str) == genus) & (statuses == control_value)
        ]
        take = min(len(arm_pos), len(arm_neg))
        if take == 0:
            continue
        pos = arm_pos.sort_values(id_column).head(take)[growth_column]
        neg = arm_neg.sort_values(id_column).head(take)[growth_column]
        if bool(pos.isna().any()) or bool(neg.isna().any()):
            continue
        out.append(
            {
                "genus": genus,
                "n_carrier": int(take),
                "n_control": int(take),
                "mean_growth_carrier": float(pos.mean()),
                "mean_growth_control": float(neg.mean()),
                "median_growth_carrier": float(pos.median()),
                "median_growth_control": float(neg.median()),
                "delta": float(pos.mean() - neg.mean()),
                "effect_unit": GROWTH_RATE_UNIT,
            }
        )
    return out


def _compartment_deltas(
    table: pd.DataFrame,
    compartment: str,
    *,
    carrier_value: str = CANDIDATE_CARRIER,
    control_value: str = CANDIDATE_NOT_DETECTED,
) -> list[dict[str, Any]]:
    """Genus deltas restricted to one *pure* compartment.

    The compartment is a property of the carrier genome (``group``); control
    rows carry ``group == 'control'`` and are always required for the matched
    arm, so they are never filtered out by the compartment label.  Carriers
    whose ``group`` is not the pure compartment value (for example multi-
    candidate ``mixed`` genomes) are excluded.
    """
    statuses = table["candidate_detection_status"].astype(str).str.strip()
    groups = table["group"].astype(str).str.strip()
    keep = ((statuses == carrier_value) & (groups == compartment)) | (
        statuses == control_value
    )
    return genus_level_deltas(table.loc[keep])


def cluster_bootstrap_mean_difference(
    table: pd.DataFrame,
    rng: np.random.Generator | None,
    *,
    n_boot: int = 2000,
    alpha: float = 0.05,
    deltas: Sequence[float] | None = None,
    genus_keys: Sequence[str] | None = None,
) -> dict[str, Any]:
    """Cluster bootstrap of the **mean genus-level difference**.

    Genera -- not genomes -- are the resampling unit: each draw takes
    ``n_genera`` genus labels with replacement and averages their genus deltas.
    Four carriers spread over two genera therefore contribute two independent
    observations, not four.  ``resampled_units`` reports the genus label drawn
    by every resample so the resampling unit is auditable, never inferred.
    """
    generator = _require_rng(rng)
    if not isinstance(n_boot, int) or n_boot < 1:
        raise ValueError("n_boot must be a positive integer")
    if not 0.0 < alpha < 1.0:
        raise ValueError("alpha must be in (0, 1)")
    if deltas is None:
        rows = genus_level_deltas(table)
        deltas = [row["delta"] for row in rows]
        if genus_keys is None:
            genus_keys = [str(row["genus"]) for row in rows]
    values = [float(value) for value in deltas]
    n_genera = len(values)
    if n_genera == 0:
        return {
            "estimate": float("nan"),
            "ci_low": float("nan"),
            "ci_high": float("nan"),
            "n_genera": 0,
            "n_resamples": 0,
            "resample_size": 0,
            "resampling_unit": RESAMPLING_UNIT,
            "alpha": alpha,
            "resampled_units": [],
        }
    array = np.asarray(values, dtype=float)
    keys = None if genus_keys is None else [str(key) for key in genus_keys]
    if keys is not None and len(keys) != n_genera:
        raise ValueError(
            f"genus_keys has {len(keys)} entries but there are {n_genera} genus deltas"
        )
    draws = np.empty(n_boot, dtype=float)
    units: list[str] = []
    for index in range(n_boot):
        picked = generator.integers(0, n_genera, size=n_genera)
        draws[index] = array[picked].mean()
        if keys is not None:
            units.extend(keys[int(i)] for i in picked)
    low, high = np.percentile(draws, [100.0 * alpha / 2.0, 100.0 * (1.0 - alpha / 2.0)])
    return {
        "estimate": float(array.mean()),
        "ci_low": float(low),
        "ci_high": float(high),
        "n_genera": int(n_genera),
        "n_resamples": int(n_boot),
        "resample_size": int(n_genera),
        "resampling_unit": RESAMPLING_UNIT,
        "alpha": float(alpha),
        "resampled_units": units,
    }


def genus_level_permutation_test(
    table: pd.DataFrame,
    rng: np.random.Generator | None,
    *,
    n_permutations: int = 10000,
    deltas: Sequence[float] | None = None,
) -> dict[str, Any]:
    """Within-genus sign-flip permutation test on the means of the genus deltas.

    A paired genus design is permuted by flipping the sign of whole genus
    deltas, which is the exact null of "carrier status is exchangeable within a
    genus".  Genomes are never shuffled: the genus stays the resampling unit.
    """
    generator = _require_rng(rng)
    if not isinstance(n_permutations, int) or n_permutations < 1:
        raise ValueError("n_permutations must be a positive integer")
    rows = genus_level_deltas(table)
    values = [float(value) for value in (deltas if deltas is not None else [r["delta"] for r in rows])]
    keys = [str(row["genus"]) for row in rows]
    n_genera = len(values)
    if n_genera == 0:
        return {
            "observed": float("nan"),
            "p_value": float("nan"),
            "n_genera": 0,
            "n_permutations": 0,
            "resampling_unit": RESAMPLING_UNIT,
            "resampled_units": [],
        }
    array = np.asarray(values, dtype=float)
    observed = float(array.mean())
    null = np.empty(n_permutations, dtype=float)
    units: list[str] = []
    for index in range(n_permutations):
        signs = generator.integers(0, 2, size=n_genera) * 2 - 1
        null[index] = float((array * signs).mean())
        units.extend(keys)
    exceed = int(np.sum(np.abs(null) >= abs(observed)))
    return {
        "observed": observed,
        "p_value": float((exceed + 1) / (n_permutations + 1)),
        "n_genera": int(n_genera),
        "n_permutations": int(n_permutations),
        "resampling_unit": RESAMPLING_UNIT,
        "resampled_units": units,
    }


def equivalence_statement(
    interval: tuple[float, float],
    *,
    mean_delta: float,
    equivalence_margin: float | None,
    alpha: float = 0.05,
) -> dict[str, Any]:
    """Bounded difference wording, allowed only with a preregistered margin.

    Without ``equivalence_margin`` the only defensible statement is
    "no difference detected under the current design"; an equivalence claim
    requires the margin to have been preregistered *before* the data.
    """
    low, high = float(interval[0]), float(interval[1])
    if equivalence_margin is None:
        return {
            "basis": "not_detected",
            "statement": (
                "no difference detected under the current design"
                if low <= 0.0 <= high
                else "a difference was detected under the current design"
            ),
            "difference_detected": not (low <= 0.0 <= high),
            "ci_low": low,
            "ci_high": high,
            "mean_delta": float(mean_delta),
            "equivalence_margin": None,
            "alpha": float(alpha),
            "effect_unit": GROWTH_RATE_UNIT,
        }
    margin = float(equivalence_margin)
    if not math.isfinite(margin) or margin <= 0:
        raise ValueError("equivalence_margin must be a positive number")
    within = (-margin <= low) and (high <= margin)
    if within:
        statement = (
            f"the genus-level mean difference is bounded within the preregistered "
            f"equivalence margin of +/-{margin:g} {GROWTH_RATE_UNIT} "
            f"({100 * (1 - alpha):g}% CI {low:+.4f} to {high:+.4f})"
        )
    else:
        statement = (
            f"the genus-level mean difference is not bounded within the "
            f"preregistered equivalence margin of +/-{margin:g} {GROWTH_RATE_UNIT} "
            f"({100 * (1 - alpha):g}% CI {low:+.4f} to {high:+.4f})"
        )
    return {
        "basis": "within_equivalence_margin" if within else "outside_equivalence_margin",
        "statement": statement,
        "difference_detected": within,
        "ci_low": low,
        "ci_high": high,
        "mean_delta": float(mean_delta),
        "equivalence_margin": margin,
        "alpha": float(alpha),
        "effect_unit": GROWTH_RATE_UNIT,
    }


def _sign_test_p(positive: int, negative: int) -> float:
    """Two-sided exact sign test (binomial, p = 0.5), no scipy dependency."""
    nonzero = positive + negative
    if nonzero == 0:
        return float("nan")
    k = min(positive, negative)
    tail = sum(math.comb(nonzero, i) for i in range(0, k + 1)) / float(2 ** nonzero)
    return float(min(1.0, 2.0 * tail))


def _mann_whitney_p(a: Sequence[float], b: Sequence[float]) -> float:
    """Two-sided Mann-Whitney U p-value with the normal approximation.

    Used only for the *secondary* genome-level description, where the
    independence assumption does not hold; the genus-level test is primary.
    """
    x = np.asarray(list(a), dtype=float)
    y = np.asarray(list(b), dtype=float)
    n_x, n_y = len(x), len(y)
    if n_x == 0 or n_y == 0:
        return float("nan")
    combined = np.concatenate([x, y])
    order = combined.argsort(kind="mergesort")
    ranks = np.empty(len(combined), dtype=float)
    sorted_values = combined[order]
    index = 0
    while index < len(sorted_values):
        stop = index
        while stop + 1 < len(sorted_values) and sorted_values[stop + 1] == sorted_values[index]:
            stop += 1
        ranks[order[index : stop + 1]] = (index + stop) / 2.0 + 1.0
        index = stop + 1
    u_stat = ranks[:n_x].sum() - n_x * (n_x + 1) / 2.0
    mu = n_x * n_y / 2.0
    sigma = math.sqrt(n_x * n_y * (n_x + n_y + 1) / 12.0)
    if sigma == 0:
        return float("nan")
    z = (u_stat - mu) / sigma
    return float(math.erfc(abs(z) / math.sqrt(2.0)))


def genus_level_estimates(
    table: pd.DataFrame,
    rng: np.random.Generator | None,
    *,
    n_boot: int = 2000,
    n_permutations: int = 10000,
    equivalence_margin: float | None = None,
    alpha: float = 0.05,
) -> dict[str, Any]:
    """Primary genus-level estimate plus the explicitly secondary analyses.

    ``primary_estimand`` is the mean genus-level difference and
    ``resampling_unit`` is ``genus``; ``secondary_analyses`` names every test
    that is *not* the primary estimand (median, sign test, genome-level
    comparison, compartment comparison) so downstream tables and figures can
    label them as secondary.
    """
    _require_columns(
        table,
        ["candidate_detection_status", "growth_rate_per_h", "genus", "genome_id"],
    )
    generator = _require_rng(rng)
    if equivalence_margin is not None:
        margin = float(equivalence_margin)
        if not math.isfinite(margin) or margin <= 0:
            raise ValueError("equivalence_margin must be a positive number")

    delta_rows = genus_level_deltas(table)
    deltas = [row["delta"] for row in delta_rows]
    genus_keys = [str(row["genus"]) for row in delta_rows]
    bootstrap = cluster_bootstrap_mean_difference(
        table,
        generator,
        n_boot=n_boot,
        alpha=alpha,
        deltas=deltas,
        genus_keys=genus_keys,
    )
    permutation = genus_level_permutation_test(
        table, generator, n_permutations=n_permutations, deltas=deltas
    )
    if not deltas:
        return {
            "primary_estimand": PRIMARY_ESTIMAND,
            "resampling_unit": RESAMPLING_UNIT,
            "analysis_role": "primary",
            "secondary_analyses": list(SECONDARY_ANALYSES),
            "n_genera": 0,
            "n_genomes": int(len(table)),
            "mean_delta": float("nan"),
            "median_delta": float("nan"),
            "sign_test_p": float("nan"),
            "genome_level_mann_whitney_p": float("nan"),
            "ci_low": float("nan"),
            "ci_high": float("nan"),
            "permutation_p": float("nan"),
            "difference_detected": False,
            "conclusion_basis": "not_detected",
            "statement": "no difference detected under the current design",
            "equivalence_margin": equivalence_margin,
            "effect_unit": GROWTH_RATE_UNIT,
            "genus_deltas": [],
            "resampled_units": [],
        }

    array = np.asarray(deltas, dtype=float)
    positive = int(np.sum(array > 0))
    negative = int(np.sum(array < 0))
    verdict = equivalence_statement(
        (bootstrap["ci_low"], bootstrap["ci_high"]),
        mean_delta=bootstrap["estimate"],
        equivalence_margin=equivalence_margin,
        alpha=alpha,
    )
    statuses = _detection_series(table)
    group_a = pd.to_numeric(
        table.loc[statuses == CANDIDATE_CARRIER, "growth_rate_per_h"], errors="coerce"
    ).dropna()
    group_b = pd.to_numeric(
        table.loc[statuses == CANDIDATE_NOT_DETECTED, "growth_rate_per_h"], errors="coerce"
    ).dropna()

    return {
        "primary_estimand": PRIMARY_ESTIMAND,
        "resampling_unit": RESAMPLING_UNIT,
        "analysis_role": "primary",
        "secondary_analyses": list(SECONDARY_ANALYSES),
        "n_genera": int(len(deltas)),
        "n_genomes": int(len(table)),
        "mean_delta": float(array.mean()),
        "median_delta": float(np.median(array)),
        "sign_test_p": _sign_test_p(positive, negative),
        "positive_genera": positive,
        "negative_genera": negative,
        "zero_genera": int(np.sum(array == 0)),
        "genome_level_mann_whitney_p": _mann_whitney_p(group_a, group_b),
        "ci_low": bootstrap["ci_low"],
        "ci_high": bootstrap["ci_high"],
        "n_bootstrap": bootstrap["n_resamples"],
        "permutation_p": permutation["p_value"],
        "n_permutations": permutation["n_permutations"],
        "difference_detected": verdict["difference_detected"],
        "conclusion_basis": verdict["basis"],
        "statement": verdict["statement"],
        "equivalence_margin": equivalence_margin,
        "effect_unit": GROWTH_RATE_UNIT,
        "genus_deltas": delta_rows,
        "resampled_units": bootstrap["resampled_units"],
    }


def intracellular_vs_extracellular(
    table: pd.DataFrame,
    rng: np.random.Generator | None,
    *,
    n_permutations: int = 10000,
) -> dict[str, Any]:
    """Secondary compartment comparison of genus-level deltas.

    Only pure-compartment genomes contribute; the resampling unit is still the
    genus.  This comparison is secondary by design (plan Task 12, Step 3).
    """
    generator = _require_rng(rng)
    _require_columns(table, ["group", "candidate_detection_status", "genus", "genome_id"])
    inner = [row["delta"] for row in _compartment_deltas(table, "intracellular")]
    outer = [row["delta"] for row in _compartment_deltas(table, "extracellular")]
    if not inner or not outer:
        return {
            "analysis_role": "secondary",
            "resampling_unit": RESAMPLING_UNIT,
            "n_intracellular_genera": len(inner),
            "n_extracellular_genera": len(outer),
            "mean_intracellular_delta": float(np.mean(inner)) if inner else float("nan"),
            "mean_extracellular_delta": float(np.mean(outer)) if outer else float("nan"),
            "observed": float("nan"),
            "p_value": float("nan"),
            "n_permutations": 0,
            "resampled_units": [],
            "effect_unit": GROWTH_RATE_UNIT,
        }
    a = np.asarray(inner, dtype=float)
    b = np.asarray(outer, dtype=float)
    observed = float(a.mean() - b.mean())
    pool = np.concatenate([a, b])
    n_a = len(a)
    units: list[str] = []
    keys = [f"intracellular:{i}" for i in range(n_a)] + [
        f"extracellular:{i}" for i in range(len(b))
    ]
    exceed = 0
    for _ in range(n_permutations):
        assignment = generator.permutation(len(pool))
        permuted_a = pool[assignment[:n_a]]
        permuted_b = pool[assignment[n_a:]]
        if abs(float(permuted_a.mean() - permuted_b.mean())) >= abs(observed):
            exceed += 1
        units.extend(keys[i] for i in assignment[:n_a])
    return {
        "analysis_role": "secondary",
        "resampling_unit": RESAMPLING_UNIT,
        "n_intracellular_genera": int(n_a),
        "n_extracellular_genera": int(len(b)),
        "mean_intracellular_delta": float(a.mean()),
        "mean_extracellular_delta": float(b.mean()),
        "observed": observed,
        "p_value": float((exceed + 1) / (n_permutations + 1)),
        "n_permutations": int(n_permutations),
        "resampled_units": units,
        "effect_unit": GROWTH_RATE_UNIT,
    }


# ------------------------------------------------------- reuse compatibility

#: every dimension that must match before a reused prediction is comparable
COMPARABILITY_FIELDS = (
    "tool",
    "tool_version",
    "model_settings",
    "temperature_handling",
    "cds_preparation",
    "ribosomal_marker_procedure",
    "output_formula",
)


def audit_reused_predictions(
    records: Iterable[Mapping[str, Any]],
    reference_profile: Mapping[str, Any],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Split reused gRodon predictions into ``(compatible, incompatible)``.

    A record is comparable only when every field of
    :data:`COMPARABILITY_FIELDS` is present and equal to the reference profile,
    i.e. the same tool version, the same model/settings, the same temperature
    handling, the same CDS preparation, the same ribosomal-marker procedure and
    the same output formula.  Each incompatible record carries
    ``mismatch_fields`` (sorted field names) and a ``reason`` that names the
    mismatched field.

    Compatible records may enter the primary analysis; incompatible ones are
    excluded from it and retained for the sensitivity table -- never silently
    dropped and never merged into the primary result.
    """
    if reference_profile is None:
        raise ValueError("reference_profile must not be None")
    missing_reference = [f for f in COMPARABILITY_FIELDS if not _clean_field(reference_profile.get(f))]
    if missing_reference:
        raise ValueError(
            f"reference_profile must state every comparability field; missing {missing_reference}"
        )

    compatible: list[dict[str, Any]] = []
    incompatible: list[dict[str, Any]] = []
    for index, record in enumerate(records):
        row = dict(record or {})
        row.setdefault("genome_id", f"record_{index}")
        mismatch_fields: list[str] = []
        reasons: list[str] = []
        for field in COMPARABILITY_FIELDS:
            expected = str(reference_profile[field]).strip()
            actual = _clean_field(row.get(field))
            if actual is None:
                mismatch_fields.append(field)
                reasons.append(f"{field}: missing (reference {expected!r})")
            elif actual != expected:
                mismatch_fields.append(field)
                reasons.append(f"{field}: {actual!r} != reference {expected!r}")
        entry = dict(row)
        entry["mismatch_fields"] = sorted(mismatch_fields)
        entry["comparable"] = not mismatch_fields
        entry["reason"] = (
            "" if not mismatch_fields else "incomparable reused prediction; " + "; ".join(reasons)
        )
        if mismatch_fields:
            incompatible.append(entry)
        else:
            compatible.append(entry)
    return compatible, incompatible


def audit_reused_predictions_table(
    records: Iterable[Mapping[str, Any]],
    reference_profile: Mapping[str, Any],
) -> dict[str, Any]:
    """Audit payload for the 713 reused predictions of the frozen run.

    Reports the counts, the primary-analysis genome ids, the sensitivity table
    (one row per incompatible record, carrying its prediction payload and the
    named reason) and the duplicate genome ids that must be resolved before the
    counts are quoted.
    """
    records = [dict(record or {}) for record in records]
    compatible, incompatible = audit_reused_predictions(records, reference_profile)
    identifiers = [str(row.get("genome_id", "")) for row in records]
    duplicates = sorted(
        {identifier for identifier, count in Counter(identifiers).items() if count > 1}
    )
    sensitivity = []
    for row in incompatible:
        sensitivity.append(
            {
                "genome_id": str(row.get("genome_id", "")),
                "doubling_time_h": row.get("doubling_time_h", row.get("d", "")),
                "growth_rate_per_h": row.get("growth_rate_per_h", ""),
                "mismatch_fields": row["mismatch_fields"],
                "reason": row["reason"],
                "used_in_primary_analysis": False,
                "used_in_sensitivity_analysis": True,
                "source_table": "reused_predictions",
                "comparability_reference": {
                    field: reference_profile[field] for field in COMPARABILITY_FIELDS
                },
            }
        )
    audit = {
        "schema_version": 1,
        "n_records": len(records),
        "n_unique_genomes": len(set(identifiers)),
        "duplicate_genome_ids": duplicates,
        "n_compatible": len(compatible),
        "n_incompatible": len(incompatible),
        "primary_genome_ids": [str(row.get("genome_id", "")) for row in compatible],
        "sensitivity_genome_ids": [str(row.get("genome_id", "")) for row in incompatible],
        "sensitivity_table": sensitivity,
        "comparability_fields": list(COMPARABILITY_FIELDS),
        "reference_profile": {field: reference_profile[field] for field in COMPARABILITY_FIELDS},
    }
    audit["status"] = (
        "ok" if audit["n_compatible"] + audit["n_incompatible"] == audit["n_records"] else "mismatch"
    )
    return audit


# ------------------------------------------------------- selection coverage

#: every disposition a genome can have in the v2 selection ledger
SELECTION_STATUSES = (
    "included",
    "excluded_no_same_genus_control",
    "excluded_missing_fasta",
    "excluded_failed_prediction",
    "excluded_invalid_growth_rate",
    "excluded_incomparable_reused_prediction",
    "excluded_not_in_balanced_design",
)

#: statuses counted as explicit reason buckets for the manifest difference
REASON_STATUS_TO_BUCKET = {
    "excluded_no_same_genus_control": "no_same_genus_control",
    "excluded_missing_fasta": "missing_fasta",
    "excluded_failed_prediction": "failed_prediction",
    "excluded_invalid_growth_rate": "invalid_growth_rate",
    "excluded_incomparable_reused_prediction": "incomparable_reused_prediction",
    "excluded_not_in_balanced_design": "not_reaching_the_balanced_same_genus_design",
}

COVERAGE_REQUIRED_COLUMNS = (
    "genome_id",
    "taxonomy",
    "candidate_family",
    "group",
    "quality",
    "selection_status",
)


def _coverage_counts(frame: pd.DataFrame) -> dict[str, int]:
    counts = {status: 0 for status in SELECTION_STATUSES}
    for status, count in Counter(frame["selection_status"].astype(str)).items():
        counts[status] = int(count)
    return counts


def _coverage_stratum(frame: pd.DataFrame) -> dict[str, Any]:
    counts = _coverage_counts(frame)
    total = int(len(frame))
    included = counts["included"]
    return {
        "total": total,
        **counts,
        "included_count": included,
        "excluded_count": total - included,
        "included_fraction": round(included / total, 6) if total else 0.0,
    }


def summarize_selection_coverage(table: pd.DataFrame) -> dict[str, Any]:
    """Quantify which genomes entered the analysis and which did not.

    The v2 selection ledger is stratified by taxonomy, candidate family, group
    and the available quality metric, with one named disposition per genome.
    Every count is computed from the accession-level rows; no count is derived
    by subtracting aggregates.
    """
    _require_columns(table, COVERAGE_REQUIRED_COLUMNS)
    unknown = sorted(
        set(table["selection_status"].astype(str)) - set(SELECTION_STATUSES)
    )
    if unknown:
        raise ValueError(
            f"unknown selection_status value(s): {unknown}; allowed values are "
            f"{list(SELECTION_STATUSES)}"
        )

    frame = table.copy()
    frame["selection_status"] = frame["selection_status"].astype(str)
    counts = _coverage_counts(frame)
    total = int(len(frame))
    included = counts["included"]

    def by(column: str) -> dict[str, dict[str, Any]]:
        return {
            str(value): _coverage_stratum(group)
            for value, group in frame.groupby(column, dropna=False, sort=True)
        }

    return {
        "schema_version": 1,
        "total_genomes": total,
        "included_count": included,
        "excluded_count": total - included,
        "status_counts": counts,
        "sum_of_status_counts": sum(counts.values()),
        "no_same_genus_control_count": counts["excluded_no_same_genus_control"],
        "missing_fasta_count": counts["excluded_missing_fasta"],
        "failed_prediction_count": counts["excluded_failed_prediction"],
        "by_taxonomy": by("taxonomy"),
        "by_candidate_family": by("candidate_family"),
        "by_group": by("group"),
        "by_quality": by("quality"),
    }


# ------------------------------------------------ manifest difference (66)

#: fields accepted as an accession key, in priority order
ACCESSION_FIELDS = ("accession", "genome_id", "genome", "gtdb_accession")


def _accession_set(
    value: Iterable[Any] | Mapping[str, Any] | pd.Series | None,
    field: str,
) -> set[str]:
    """Normalise any accession-level input into a set of accession strings.

    Accepted shapes are a mapping (its keys are the accessions), a one-column
    table with a recognised accession column, a pandas Series, or any iterable
    of accession strings.  A bare string is rejected: it would silently explode
    into single characters.
    """
    if value is None:
        raise ValueError(f"{field} must not be None")
    if isinstance(value, (str, bytes)):
        raise ValueError(
            f"{field} must be a collection of accessions, not a single string"
        )
    if isinstance(value, pd.DataFrame):
        for candidate in ACCESSION_FIELDS:
            if candidate in value.columns:
                return {str(item).strip() for item in value[candidate] if str(item).strip()}
        raise ValueError(f"{field} has no accession column among {list(ACCESSION_FIELDS)}")
    if isinstance(value, Mapping):
        return {str(key).strip() for key in value if str(key).strip()}
    return {str(item).strip() for item in value if str(item).strip()}


def explain_manifest_difference(
    *,
    eligible_positive_accessions: Iterable[Any] | Mapping[str, Any],
    manifest_positive_accessions: Iterable[Any] | Mapping[str, Any],
    reason_buckets: Mapping[str, Iterable[Any]] | None,
    expected_eligible_count: int | None = None,
    expected_manifest_count: int | None = None,
    expected_difference: int | None = None,
) -> dict[str, Any]:
    """Explain the eligible-positive to manifest-positive difference by set algebra.

    The 2026-09-21 review records 4,507 eligible positives against 4,441
    manifest positives -- a difference of 66 genomes.  This function computes
    every number from **accession-level sets**:

    * ``difference = len(eligible - manifest)``;
    * each ``reason_buckets[name]`` is a set of accessions;
    * ``reason_bucket_counts`` are set sizes, and their sum is required to equal
      ``difference`` exactly;
    * any eligible genome missing from the manifest and from every bucket is
      listed in ``unexplained_accessions`` and makes ``status`` a
      ``"mismatch"``;
    * an accession claimed by two buckets is a ``"mismatch"`` (bucket overlap),
      so the buckets must be a partition, not a collection of overlapping sets;
    * accessions named in a bucket but absent from the eligible set, and
      manifest positives absent from the eligible set, are reported separately.

    Aggregates are never subtracted from one another: nothing here can produce a
    number that does not come from an accession identity.
    """
    eligible = _accession_set(eligible_positive_accessions, "eligible_positive_accessions")
    manifest = _accession_set(manifest_positive_accessions, "manifest_positive_accessions")
    if not eligible:
        raise ValueError("eligible_positive_accessions must not be empty")

    difference_set = eligible - manifest
    not_eligible = manifest - eligible

    buckets: dict[str, set[str]] = {}
    for name, members in (reason_buckets or {}).items():
        buckets[str(name)] = _accession_set(members, f"reason_buckets[{name!r}]")

    counts = {name: len(members) for name, members in buckets.items()}
    union: set[str] = set()
    overlap: set[str] = set()
    for members in buckets.values():
        overlap |= union & members
        union |= members

    unexplained = difference_set - union
    outside = union - eligible
    sum_of_counts = sum(counts.values())
    difference = len(difference_set)

    counts_ok = True
    count_failures: list[str] = []
    if expected_eligible_count is not None and int(expected_eligible_count) != len(eligible):
        counts_ok = False
        count_failures.append("eligible_positive_count")
    if expected_manifest_count is not None and int(expected_manifest_count) != len(manifest):
        counts_ok = False
        count_failures.append("manifest_positive_count")
    if expected_difference is not None and int(expected_difference) != difference:
        counts_ok = False
        count_failures.append("difference")

    status = (
        "ok"
        if (
            sum_of_counts == difference
            and not unexplained
            and not overlap
            and not outside
            and not not_eligible
            and counts_ok
        )
        else "mismatch"
    )
    return {
        "schema_version": 1,
        "source": "accession_level_sets",
        "status": status,
        "eligible_positive_count": len(eligible),
        "manifest_positive_count": len(manifest),
        "difference": difference,
        "expected_eligible_count": expected_eligible_count,
        "expected_manifest_count": expected_manifest_count,
        "expected_difference": expected_difference,
        "count_expectation_failures": count_failures,
        "reason_bucket_counts": counts,
        "reason_bucket_accessions": {
            name: sorted(members) for name, members in sorted(buckets.items())
        },
        "sum_of_reason_bucket_counts": sum_of_counts,
        "buckets_partition_the_difference": sum_of_counts == difference and not overlap,
        "overlap_count": len(overlap),
        "overlap_accessions": sorted(overlap),
        "unexplained_count": len(unexplained),
        "unexplained_accessions": sorted(unexplained),
        "bucket_members_outside_eligible_count": len(outside),
        "bucket_members_outside_eligible_accessions": sorted(outside),
        "manifest_not_eligible_count": len(not_eligible),
        "manifest_not_eligible_accessions": sorted(not_eligible),
    }
