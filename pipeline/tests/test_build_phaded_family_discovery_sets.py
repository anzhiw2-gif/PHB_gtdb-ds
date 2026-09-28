"""Tests for the per-family discovery-set builder (F3, test-first).

Discovery sets exist to **RECALL**: a family training set plus its provenance is
built so that a recall-only discovery HMM can be trained later.  A discovery HMM
has no discriminating power (78.69% sensitivity on the 563 x1 confounder hits)
and therefore may never filter, drop, demote or exclude a candidate, and may
never produce a family call.

Load-bearing invariants pinned here:

* the family selector is keyed on ``model_layer`` (``reference_query_only``:
  no HMM exists yet, so a discovery set is built), never on the legacy
  ``model_status`` column;
* families that already have a *discriminating* layer
  (``sequence_family_hmm_validated`` / ``calibrated_candidate_model``) are
  excluded from discovery-set construction — they do not need a recall-only
  model, and a recall-only model must never be presented as their model;
* a manifest that predates ``model_layer`` (column absent) is read as
  ``reference_query_only`` / non-discriminating and says so — "trained in the
  legacy column" is never read as "validated";
* the builder only *adds* training sets and provenance rows: it has no code path
  that can exclude, drop or demote any candidate, and every emitted ledger
  accession must appear in the emitted training FASTA;
* ``model_layer`` and ``functional_calibration_status`` are independent.

Fixtures are synthetic (temp directory); no project data is read.
"""

from __future__ import annotations

import csv
import hashlib
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import build_phaded_family_discovery_sets as module  # noqa: E402


REFERENCE_ONLY = "reference_query_only"
DISCOVERY = "discovery_hmm_uncalibrated"
VALIDATED = "sequence_family_hmm_validated"
CALIBRATED = "calibrated_candidate_model"

LAYERED_COLUMNS = [
    "profile_id", "profile_kind", "phaded_superfamily", "phaded_family_id",
    "model_layer", "model_status", "functional_calibration_status",
]
LEGACY_COLUMNS = [
    "profile_id", "profile_kind", "phaded_superfamily", "phaded_family_id",
    "model_status", "model_reason",
]

LEDGER_COLUMNS = [
    "accession", "phaded_family_id", "evidence_status", "pmid", "primary_doi",
]

#: The five already-trained families of the frozen registry (mirrors, never
#: reads, the project values).
TRAINED_FAMILIES = sorted(module.TRAINED_FAMILIES)


def family_rows(prefix: str, count: int, *, offset: int = 0) -> list[str]:
    return [f"{prefix}{index:02d}" for index in range(offset, offset + count)]


def make_profiles(
    *,
    reference_only_families: list[str],
    validated_families: list[str],
    discovery_families: list[str],
    calibrated_families: list[str] | None = None,
) -> list[dict[str, str]]:
    profiles: list[dict[str, str]] = []
    for family in reference_only_families:
        profiles.append({
            "profile_id": f"family_{family}_hash", "profile_kind": "family",
            "phaded_superfamily": f"sf_{family}", "phaded_family_id": family,
            "model_layer": REFERENCE_ONLY, "model_status": "reference_only",
            "functional_calibration_status": "not_function_calibrated",
        })
    for family in validated_families:
        profiles.append({
            "profile_id": f"family_{family}_hash", "profile_kind": "family",
            "phaded_superfamily": f"sf_{family}", "phaded_family_id": family,
            # legacy word is ALSO "trained": the trap this task closes
            "model_layer": VALIDATED, "model_status": "trained",
            "functional_calibration_status": "candidate_gate_passed_not_promoted",
        })
    for family in discovery_families:
        profiles.append({
            "profile_id": f"family_{family}_hash", "profile_kind": "family",
            "phaded_superfamily": f"sf_{family}", "phaded_family_id": family,
            "model_layer": DISCOVERY, "model_status": "trained",
            "functional_calibration_status": "not_function_calibrated",
        })
    for family in calibrated_families or []:
        profiles.append({
            "profile_id": f"family_{family}_hash", "profile_kind": "family",
            "phaded_superfamily": f"sf_{family}", "phaded_family_id": family,
            "model_layer": CALIBRATED, "model_status": "trained",
            "functional_calibration_status": "calibrated_candidate_model",
        })
    # superfamily-level rows have no family id and are never selectable
    profiles.append({
        "profile_id": "superfamily_sf_x_hash", "profile_kind": "superfamily",
        "phaded_superfamily": "sf_x", "phaded_family_id": "",
        "model_layer": REFERENCE_ONLY, "model_status": "reference_only",
        "functional_calibration_status": "not_function_calibrated",
    })
    return profiles


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t",
                                lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})
    return path


class Fixture(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="dsh_family_discovery_f3_"))
        self.ledger = self.tmp / "ledger.tsv"
        self.fasta = self.tmp / "reference.faa"
        self.out_dir = self.tmp / "out"

    def build(self, profiles, columns=LAYERED_COLUMNS, out_dir=None):
        manifest = self.tmp / "profile_manifest.tsv"
        write_tsv(manifest, columns, profiles)
        ledger_rows: list[dict[str, str]] = []
        fasta_lines: list[str] = []
        family_ids = sorted({
            row["phaded_family_id"] for row in profiles if row.get("phaded_family_id")
        })
        for family in family_ids:
            accession = f"ACC_{family}"
            ledger_rows.append({
                "accession": accession, "phaded_family_id": family,
                "evidence_status": "experimental_positive", "pmid": f"PMID_{family}",
                "primary_doi": f"10.0000/{family}",
            })
            fasta_lines.append(f">ref_{accession}|{accession} {family}")
            fasta_lines.append("MKWV")
        write_tsv(self.ledger, LEDGER_COLUMNS, ledger_rows)
        self.fasta.write_text("\n".join(fasta_lines) + "\n", encoding="utf-8", newline="\n")
        return module.build(self.ledger, self.fasta, manifest,
                            out_dir or self.out_dir)

    def bundle(self, profiles, **kwargs):
        """Build and read back the emitted training FASTA + provenance manifest."""
        summary = self.build(profiles, **kwargs)
        manifest_path = Path(summary["provenance_manifest"])
        with manifest_path.open(encoding="utf-8", newline="") as handle:
            provenance = list(csv.DictReader(handle, delimiter="\t"))
        emitted: dict[str, set[str]] = {}
        for row in provenance:
            faa = manifest_path.parent / "training" / row["training_faa"]
            emitted[row["phaded_family_id"]] = {
                line[1:].split()[0]
                for line in faa.read_text(encoding="utf-8").splitlines()
                if line.startswith(">")
            }
        return summary, provenance, emitted


class LayerVocabularyTests(unittest.TestCase):
    def test_layer_vocabulary_is_the_four_governed_layers(self):
        self.assertEqual(module.MODEL_LAYERS,
                         {REFERENCE_ONLY, DISCOVERY, VALIDATED, CALIBRATED})

    def test_only_the_no_hmm_layer_needs_a_discovery_set(self):
        self.assertEqual(module.DISCOVERY_SET_MODEL_LAYERS, {REFERENCE_ONLY})
        self.assertEqual(module.EXCLUDED_DISCOVERY_SET_MODEL_LAYERS,
                         {VALIDATED, CALIBRATED})
        # a recall-only discovery HMM is never a reason to build another
        # discovery set, and never a discriminating family model
        self.assertNotIn(DISCOVERY, module.DISCOVERY_SET_MODEL_LAYERS)
        self.assertNotIn(DISCOVERY, module.EXCLUDED_DISCOVERY_SET_MODEL_LAYERS)

    def test_unrecognized_declared_layer_fails_closed(self):
        with self.assertRaisesRegex(ValueError, "model_layer"):
            module.resolve_model_layer({"profile_id": "p1", "model_layer": "trained"})

    def test_legacy_column_is_not_a_layer_fallback(self):
        self.assertEqual(module.resolve_model_layer({"profile_id": "p1",
                                                     "model_status": "trained"}),
                         REFERENCE_ONLY)


class SelectionTests(Fixture):
    def setUp(self):
        super().setUp()
        self.reference_families = family_rows("DED_hfam_", 33)
        self.validated_families = TRAINED_FAMILIES
        self.profiles = make_profiles(
            reference_only_families=self.reference_families,
            validated_families=self.validated_families,
            discovery_families=["DED_hfam_900"],
            calibrated_families=["DED_hfam_901"],
        )
        self.manifest = write_tsv(self.tmp / "profile_manifest.tsv",
                                  LAYERED_COLUMNS, self.profiles)

    def selected(self):
        return module.select_discovery_set_families(
            module.read_profile_manifest(self.manifest)[0]
        )

    def test_reference_query_only_families_are_selected(self):
        families = self.selected()
        self.assertEqual(len(families), 33)
        self.assertEqual(families, set(self.reference_families))

    def test_validated_family_layer_needs_no_discovery_set(self):
        families = self.selected()
        for family in self.validated_families:
            self.assertNotIn(family, families)

    def test_a_legacy_trained_word_never_selects_a_family(self):
        for family in self.validated_families:
            self.assertIn(family, module.TRAINED_FAMILIES)
        families = self.selected()
        # the five rows whose model_status says "trained" are excluded because
        # their declared layer is validated, not because of the legacy word
        self.assertEqual(len(families), 33)

    def test_selector_excludes_only_the_no_hmm_layer(self):
        """Every profile row is either selected or explicitly classified."""
        rows, report = module.read_profile_manifest(self.manifest)
        selected = self.selected()
        counts = report["model_layer_counts"]
        self.assertEqual(counts[DISCOVERY], 1)
        self.assertEqual(counts[CALIBRATED], 1)
        self.assertEqual(counts[VALIDATED], 5)
        self.assertEqual(len(rows), counts[REFERENCE_ONLY] + 7)
        # the discovery-only family is NOT selected: it already has a model
        self.assertNotIn("DED_hfam_900", selected)
        # ...and it is not in the excluded discriminating set either — it is a
        # third, recall-only class that simply needs no new discovery set
        self.assertNotIn("DED_hfam_900", module.EXCLUDED_DISCOVERY_SET_MODEL_LAYERS)
        self.assertEqual(
            selected,
            {row["phaded_family_id"] for row in rows
             if row["model_layer"] == REFERENCE_ONLY and row["phaded_family_id"]},
        )

    def test_superfamily_level_rows_are_never_selected(self):
        self.assertTrue(all(family.startswith("DED_hfam_") for family in self.selected()))


class LegacyManifestTests(Fixture):
    def test_absent_layer_column_is_read_as_non_discriminating(self):
        profiles = make_profiles(
            reference_only_families=family_rows("DED_hfam_", 28),
            validated_families=TRAINED_FAMILIES,
            discovery_families=[],
        )
        manifest = write_tsv(self.tmp / "legacy_manifest.tsv", LEGACY_COLUMNS, profiles)
        rows, report = module.read_profile_manifest(manifest)
        self.assertFalse(report["model_layer_column_present"])
        self.assertEqual(report["layer_source"], "absent_model_layer_column")
        self.assertEqual(report["fallback_layer"], REFERENCE_ONLY)
        self.assertIn("non-discriminating", report["fallback_note"])
        self.assertIn("never assumed validated", report["fallback_note"])
        self.assertEqual({row["model_layer"] for row in rows}, {REFERENCE_ONLY})
        # a legacy "trained" row is selected like any other no-HMM row: the
        # builder is additive and a legacy manifest must keep parsing
        families = module.select_discovery_set_families(rows)
        self.assertEqual(families, set(family_rows("DED_hfam_", 28)) | set(TRAINED_FAMILIES))

    def test_legacy_and_layered_manifests_select_the_same_recall_set(self):
        """The layer semantics do not change *which* families need a set here;
        both forms read no-HMM families, and neither treats a legacy "trained"
        word as a validated layer (`model_status` is not a layer fallback)."""
        reference_families = family_rows("DED_hfam_", 28)
        legacy_profiles = make_profiles(reference_only_families=reference_families,
                                        validated_families=[],
                                        discovery_families=[])
        legacy_profiles += [{
            "profile_id": f"family_{family}_hash", "profile_kind": "family",
            "phaded_superfamily": f"sf_{family}", "phaded_family_id": family,
            "model_status": "trained",
        } for family in TRAINED_FAMILIES]
        legacy_manifest = write_tsv(self.tmp / "legacy_manifest.tsv",
                                    LEGACY_COLUMNS, legacy_profiles)
        layered_manifest = write_tsv(
            self.tmp / "layered_manifest.tsv", LAYERED_COLUMNS,
            make_profiles(reference_only_families=reference_families,
                          validated_families=TRAINED_FAMILIES,
                          discovery_families=[]),
        )
        legacy_rows, legacy_report = module.read_profile_manifest(legacy_manifest)
        layered_rows, layered_report = module.read_profile_manifest(layered_manifest)
        legacy = module.select_discovery_set_families(legacy_rows)
        layered = module.select_discovery_set_families(layered_rows)
        # a legacy "trained" row falls back to the non-discriminating layer and
        # stays in the recall set (never silently dropped)...
        self.assertEqual(legacy, set(reference_families) | set(TRAINED_FAMILIES))
        self.assertEqual(legacy_report["fallback_layer"], REFERENCE_ONLY)
        # ...while an explicit validated layer removes the family from discovery
        # -set construction, because a recall-only model must never stand in for
        # an already validated family model
        self.assertEqual(layered, set(reference_families))
        self.assertEqual(layered_report["fallback_note"], "")


class BuildTests(Fixture):
    def setUp(self):
        super().setUp()
        self.reference_families = family_rows("DED_hfam_", 33)
        self.profiles = make_profiles(
            reference_only_families=self.reference_families,
            validated_families=TRAINED_FAMILIES,
            discovery_families=["DED_hfam_900"],
            calibrated_families=["DED_hfam_901"],
        )

    def test_build_emits_one_training_set_per_reference_only_family(self):
        summary, provenance, emitted = self.bundle(self.profiles)
        self.assertEqual(summary["reference_only_families"], 33)
        self.assertEqual(summary["training_fasta_files"], 33)
        self.assertEqual({row["phaded_family_id"] for row in provenance},
                         set(self.reference_families))
        self.assertEqual(summary["decision_field"], "model_layer")
        self.assertFalse(summary["legacy_model_status_is_a_decision_input"])
        self.assertEqual(summary["discovery_set_model_layers"], [REFERENCE_ONLY])
        self.assertEqual(summary["excluded_trained_families"], TRAINED_FAMILIES)
        self.assertTrue(summary["recall_only"])
        self.assertFalse(summary["filters_candidates"])

    def test_discovery_only_profile_never_excludes_a_candidate(self):
        """A recall-only model may not cause any ledger accession to be dropped."""
        summary, provenance, emitted = self.bundle(self.profiles)
        for family in self.reference_families:
            self.assertIn(f"ACC_{family}", emitted[family])
        # every emitted training FASTA row is a ledger accession, and the counts
        # are consistent — nothing is silently filtered out
        for row in provenance:
            self.assertEqual(row["n_sequences"], "1")
            self.assertEqual(row["n_experimental_positive"], "1")
            self.assertEqual(row["n_annotation_only"], "0")
            self.assertIn(f"ACC_{row['phaded_family_id']}", emitted[row["phaded_family_id"]])
            faa = Path(summary["provenance_manifest"]).parent / "training" / row["training_faa"]
            self.assertEqual(hashlib.sha256(faa.read_bytes()).hexdigest(),
                             row["training_faa_sha256"])

    def test_builder_has_no_candidate_exclusion_path(self):
        """The builder emits training sets; it must contain no filter/delete path."""
        source = (Path(__file__).resolve().parents[1] / "scripts"
                  / "build_phaded_family_discovery_sets.py").read_text(encoding="utf-8")
        for forbidden in ("primary_disposition", "excluded_input_quality",
                          "absence_statement", "emit_family_call", "sequence_family_call"):
            self.assertNotIn(forbidden, source)
        # no dataset-level removal API: the only writes are training FASTA files
        # plus the provenance manifest
        self.assertNotIn(".unlink(", source)
        self.assertNotIn("shutil", source)
        self.assertIn("recall-only", source)
        # the exclusion decision is a layer decision, and it excludes *model
        # construction*, never a candidate
        self.assertEqual(module.EXCLUDED_DISCOVERY_SET_MODEL_LAYERS,
                         {"sequence_family_hmm_validated", "calibrated_candidate_model"})

    def test_candidate_ledger_ids_survive_the_recall_set_construction(self):
        """Every family the ledger knows keeps its own emitted training set.

        A discovery-only profile elsewhere in the manifest cannot make the
        builder drop a ledger accession: the selector reads the profile layer,
        the ledger is only ever grouped (never filtered).
        """
        summary, provenance, emitted = self.bundle(self.profiles)
        ledger_families = {
            line.split("\t")[1]
            for line in self.ledger.read_text(encoding="utf-8").splitlines()[1:]
            if line.strip()
        }
        # the ledger also holds the validated/discovery/calibrated families, which
        # get no discovery set — and their accessions are simply not part of this
        # additive output, never "excluded candidates"
        self.assertTrue(set(self.reference_families) <= ledger_families)
        for family in self.reference_families:
            self.assertEqual(emitted[family], {f"ACC_{family}"})
        self.assertEqual(summary["filters_candidates"], False)

    def test_manifest_without_layer_column_still_parses(self):
        """A pre-``model_layer`` manifest parses and is never read as validated."""
        legacy_profiles = make_profiles(
            reference_only_families=family_rows("DED_hfam_", 28),
            validated_families=[], discovery_families=[],
        )
        legacy_profiles += [{
            "profile_id": f"family_{family}_hash", "profile_kind": "family",
            "phaded_superfamily": f"sf_{family}", "phaded_family_id": family,
            "model_status": "trained",
        } for family in TRAINED_FAMILIES]
        legacy_manifest = write_tsv(self.tmp / "legacy_manifest.tsv",
                                    LEGACY_COLUMNS, legacy_profiles)
        rows, report = module.read_profile_manifest(legacy_manifest)
        self.assertFalse(report["model_layer_column_present"])
        self.assertEqual(report["layer_source"], "absent_model_layer_column")
        self.assertEqual(report["fallback_layer"], REFERENCE_ONLY)
        self.assertIn("never assumed validated", report["fallback_note"])
        self.assertEqual({row["model_layer"] for row in rows}, {REFERENCE_ONLY})
        for row in rows:
            # the legacy word is carried through verbatim and is never the decision
            self.assertIn(row["model_status"], {"trained", "reference_only"})
            self.assertEqual(row["functional_calibration_status"], "unknown")

    def test_legacy_manifest_builds_and_reports_the_fallback(self):
        """A legacy manifest whose rows are all no-HMM builds, and says why."""
        profiles = make_profiles(reference_only_families=family_rows("DED_hfam_", 33),
                                 validated_families=[], discovery_families=[])
        summary, provenance, _ = self.bundle(profiles, columns=LEGACY_COLUMNS)
        self.assertEqual(summary["reference_only_families"], 33)
        self.assertEqual(summary["layer_source"], "absent_model_layer_column")
        self.assertEqual(summary["fallback_layer"], REFERENCE_ONLY)
        self.assertIn("never assumed validated", summary["fallback_note"])
        self.assertEqual({row["model_layer"] for row in provenance},
                         {module.DISCOVERY_LAYER})
        self.assertEqual({row["recall_only"] for row in provenance}, {"true"})

    def test_legacy_manifest_with_trained_rows_fails_closed_on_the_count(self):
        """No row is silently dropped to make an expected count fit."""
        profiles = make_profiles(reference_only_families=family_rows("DED_hfam_", 33),
                                 validated_families=[], discovery_families=[])
        profiles += [{
            "profile_id": f"family_{family}_hash", "profile_kind": "family",
            "phaded_superfamily": f"sf_{family}", "phaded_family_id": family,
            "model_status": "trained",
        } for family in TRAINED_FAMILIES]
        manifest = write_tsv(self.tmp / "legacy_manifest.tsv", LEGACY_COLUMNS, profiles)
        write_tsv(self.ledger, LEDGER_COLUMNS, [])
        self.fasta.write_text("", encoding="utf-8", newline="\n")
        with self.assertRaisesRegex(ValueError, "reference-only families, got 38"):
            module.build(self.ledger, self.fasta, manifest, self.out_dir)

    def test_functional_calibration_is_carried_through_never_inferred(self):
        summary, provenance, _ = self.bundle(self.profiles)
        for row in provenance:
            self.assertEqual(row["functional_calibration_status"],
                             "not_function_calibrated")
        self.assertNotIn("calibrated_candidate_model",
                         {row["functional_calibration_status"] for row in provenance})


if __name__ == "__main__":
    unittest.main()
