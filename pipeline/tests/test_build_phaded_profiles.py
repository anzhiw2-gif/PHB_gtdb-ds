import csv
import hashlib
import importlib.util
import tempfile
import unittest
from pathlib import Path
from unittest import mock


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "build_phaded_profiles.py"


FIELDS = [
    "reference_id", "accession", "sequence_sha256", "sequence_source",
    "source_database", "source_database_version", "retrieval_date", "organism",
    "taxonomy_id", "phaded_superfamily", "phaded_family_id",
    "reported_localization", "substrate_class", "substrate_detail",
    "experimental_assay", "experimental_result", "evidence_status",
    "positive_negative_control", "catalytic_residues", "lipase_box_or_ahsmg",
    "oxyanion_hole_evidence", "domain_architecture", "primary_doi", "pmid",
    "pmcid", "notes",
]

#: Optional evidence column introduced by the 2026-09-28 evidence-model redesign.
#: The frozen 723-row ledger does not carry it, so every path must also work when it
#: is absent (eligibility is then derived from ``evidence_status``).
ELIGIBILITY_FIELD = "discovery_training_eligible"

DEFAULT_SUPERFAMILY = "extracellular dPHASCL type 1"
WITH_LIPASE_SUPERFAMILY = "intracellular nPHASCL with lipase box"
FAKE_VERSIONS = {"mafft": "fake-mafft-7.525", "hmmer": "fake-hmmer-3.4"}
DISCOVERY_LAYER = "discovery_hmm_uncalibrated"
REFERENCE_LAYER = "reference_query_only"
VALIDATED_LAYERS = {"sequence_family_hmm_validated", "calibrated_candidate_model"}
BLOCKED_FUNCTIONAL_STATUS = "blocked_contradictory_experimental_negative"


def _load_module():
    spec = importlib.util.spec_from_file_location("build_phaded_profiles", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_inputs(root, rows, *, extra_fields=()):
    fasta = root / "reference.faa"
    ledger = root / "ledger.tsv"
    fields = FIELDS + [field for field in extra_fields if field not in FIELDS]
    sequences = {}
    for index, row in enumerate(rows, start=1):
        sequence = row.get("sequence", "M" + "A" * (19 + index))
        accession = row["accession"]
        reference_id = row.get("reference_id", f"ref-{index}")
        sequences[reference_id] = (accession, sequence)
        row["reference_id"] = reference_id
        row["sequence_sha256"] = hashlib.sha256(sequence.encode("ascii")).hexdigest()
    fasta.write_text(
        "".join(f">{ref}|{accession}\n{sequence}\n" for ref, (accession, sequence) in sequences.items()),
        encoding="ascii",
    )
    with ledger.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            output = {field: "" for field in fields}
            output.update(row)
            output.pop("sequence", None)
            writer.writerow(output)
    return ledger, fasta


def blank_sequence_sha256(ledger, accession):
    """Blank one ledger hash so a claimed-eligible row has no sequence hash at all."""
    with ledger.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        rows = list(reader)
    for row in rows:
        if row["accession"] == accession:
            row["sequence_sha256"] = ""
    with ledger.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def row(accession, family="family-1", superfamily=DEFAULT_SUPERFAMILY):
    return {
        "accession": accession,
        "phaded_superfamily": superfamily,
        "phaded_family_id": family,
        "evidence_status": "experimental_positive",
        "experimental_assay": "purified enzyme assay",
        "experimental_result": "positive",
    }


def annotation_row(accession, family="family-1", superfamily=DEFAULT_SUPERFAMILY):
    """A quality-controlled annotation-only reference: no direct experiment."""
    return {
        "accession": accession,
        "phaded_superfamily": superfamily,
        "phaded_family_id": family,
        "evidence_status": "annotation_only",
    }


def negative_row(accession, family="family-1", superfamily=DEFAULT_SUPERFAMILY):
    """A direct experimental negative: contradictory support for its family."""
    return {
        "accession": accession,
        "phaded_superfamily": superfamily,
        "phaded_family_id": family,
        "evidence_status": "experimental_negative",
        "experimental_assay": "purified enzyme assay",
        "experimental_result": "negative",
    }


def family_profile(result):
    return next(item for item in result["profiles"] if item["profile_kind"] == "family")


def superfamily_profile(result):
    return next(item for item in result["profiles"] if item["profile_kind"] == "superfamily")


def read_tsv(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class PhaDEDProfileTests(unittest.TestCase):
    def setUp(self):
        self._temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self._temporary.cleanup)
        self.root = Path(self._temporary.name)
        self._build_count = 0

    def build_with_fake_tools(self, rows, *, extra_fields=(), software_versions=None, mafft_threads=None):
        """Build profiles from synthetic rows with stubbed MAFFT/hmmbuild.

        No real alignment or HMM tool is ever executed: the stub copies the training
        FASTA into the alignment file and writes a fake HMM, exactly like the
        pre-existing tests in this module.
        """
        module = _load_module()
        ledger, fasta = write_inputs(self.root, rows, extra_fields=extra_fields)
        self._build_count += 1
        output_dir = self.root / f"profiles-{self._build_count}"

        def fake_build(training, alignment, hmm, mafft, hmmbuild, **kwargs):
            alignment.write_text(training.read_text(encoding="ascii"), encoding="ascii")
            hmm.write_text("fake-hmm\n", encoding="ascii")

        with mock.patch.object(module, "build_hmm", side_effect=fake_build):
            return module.build_profiles(
                ledger, fasta, output_dir, expected_family_count=None,
                software_versions=FAKE_VERSIONS if software_versions is None else software_versions,
                mafft_threads=mafft_threads,
            )

    def test_build_hmm_binds_mafft_thread_limit(self):
        module = _load_module()
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            training = root / "training.faa"
            alignment = root / "alignment.faa"
            hmm = root / "profile.hmm"
            training.write_text(">P1\nMAAA\n", encoding="ascii")
            calls = []

            def fake_run(command):
                calls.append(command)
                if command[0] == "mafft":
                    return mock.Mock(stdout=">P1\nMAAA\n")
                alignment.write_text("alignment\n", encoding="ascii")
                hmm.write_text("hmm\n", encoding="ascii")
                return mock.Mock(stdout="")

            with mock.patch.object(module, "_run", side_effect=fake_run):
                module.build_hmm(training, alignment, hmm, "mafft", "hmmbuild", mafft_threads=40)
            self.assertIn("--thread", calls[0])
            self.assertEqual(calls[0][calls[0].index("--thread") + 1], "40")

    def test_rejects_mixed_superfamily_family_group(self):
        module = _load_module()
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            ledger, fasta = write_inputs(root, [
                row("P1", family="same-family", superfamily="extracellular dPHASCL type 1"),
                row("P2", family="same-family", superfamily="intracellular nPHASCL without lipase box"),
            ])
            with self.assertRaisesRegex(ValueError, "mixed superfamily"):
                module.build_profiles(ledger, fasta, root / "profiles", expected_family_count=None)

    def test_singleton_family_is_reference_only_with_empty_hmm_path(self):
        module = _load_module()
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            ledger, fasta = write_inputs(root, [row("P1")])
            result = module.build_profiles(ledger, fasta, root / "profiles", expected_family_count=None)
            self.assertEqual(result["family_status"], "reference_only")
            self.assertEqual(result["hmm_path"], "")
            self.assertEqual(result["profiles"][0]["model_status"], "reference_only")

    def test_manifest_has_deterministic_training_accessions_and_required_fields(self):
        module = _load_module()
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            ledger, fasta = write_inputs(root, [row("P3"), row("P1"), row("P2")])
            def fake_build(training, alignment, hmm, mafft, hmmbuild, **kwargs):
                alignment.write_text(training.read_text(encoding="ascii"), encoding="ascii")
                hmm.write_text("fake-hmm\n", encoding="ascii")
            with mock.patch.object(module, "build_hmm", side_effect=fake_build):
                result = module.build_profiles(
                    ledger, fasta, root / "profiles", expected_family_count=None,
                    software_versions={"mafft": "fake-mafft", "hmmer": "fake-hmmer"},
                )
            profile = next(item for item in result["profiles"] if item["profile_kind"] == "family")
            self.assertEqual(profile["training_accessions"], "P1;P2;P3")
            for field in (
                "profile_id", "phaded_superfamily", "phaded_family_id", "training_accessions",
                "training_count", "model_status", "alignment_sha256", "hmm_sha256",
                "mafft_version", "hmmer_version", "hmm_path",
                "model_layer", "functional_calibration_status", "training_sequence_count",
                "experimental_anchor_count", "annotation_only_count", "training_set_sha256",
                "bit_reproducible",
            ):
                self.assertIn(field, profile)

    def test_build_profiles_passes_mafft_threads_to_builder(self):
        module = _load_module()
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            ledger, fasta = write_inputs(root, [row("P1"), row("P2"), row("P3")])
            calls = []

            def fake_build(training, alignment, hmm, mafft, hmmbuild, **kwargs):
                calls.append(kwargs.get("mafft_threads"))
                alignment.write_text(training.read_text(encoding="ascii"), encoding="ascii")
                hmm.write_text("fake-hmm\n", encoding="ascii")

            with mock.patch.object(module, "build_hmm", side_effect=fake_build):
                module.build_profiles(
                    ledger, fasta, root / "profiles", expected_family_count=None,
                    software_versions={"mafft": "fake-mafft", "hmmer": "fake-hmmer"},
                    mafft_threads=40,
                )
            self.assertEqual(calls, [40, 40])

    def test_rejects_duplicate_sequence_hash_in_training_group(self):
        module = _load_module()
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            duplicate = row("P2")
            duplicate["sequence"] = "M" + "A" * 20
            first = row("P1")
            first["sequence"] = duplicate["sequence"]
            ledger, fasta = write_inputs(root, [first, duplicate])
            with self.assertRaisesRegex(ValueError, "duplicate sequence hash"):
                module.build_profiles(ledger, fasta, root / "profiles", expected_family_count=None)

    def test_manifest_has_unique_ids_for_slug_collisions(self):
        module = _load_module()
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            rows = [
                row(f"P{i}", family=family)
                for family, start in (("family-a", 1), ("family_a", 4))
                for i in range(start, start + 3)
            ]
            ledger, fasta = write_inputs(root, rows)
            def fake_build(training, alignment, hmm, mafft, hmmbuild, **kwargs):
                alignment.write_text(training.read_text(encoding="ascii"), encoding="ascii")
                hmm.write_text("fake-hmm\n", encoding="ascii")
            with mock.patch.object(module, "build_hmm", side_effect=fake_build):
                result = module.build_profiles(
                    ledger, fasta, root / "profiles", expected_family_count=None,
                    software_versions={"mafft": "fake-mafft", "hmmer": "fake-hmmer"},
                )
            with (root / "profiles" / "profile_manifest.tsv").open(encoding="utf-8", newline="") as handle:
                manifest = list(csv.DictReader(handle, delimiter="\t"))
            ids = [item["profile_id"] for item in manifest]
            paths = [item["hmm_path"] for item in manifest]
            self.assertEqual(len(ids), len(set(ids)))
            self.assertEqual(len(paths), len(set(paths)))
            self.assertEqual(len(result["profiles"]), 3)

    def test_annotation_only_family_builds_discovery_hmm_when_three_eligible_sequences_exist(self):
        rows = [annotation_row("A1"), annotation_row("A2"), annotation_row("A3")]
        result = self.build_with_fake_tools(rows)
        family = family_profile(result)
        self.assertEqual(family["model_layer"], "discovery_hmm_uncalibrated")
        self.assertEqual(family["experimental_anchor_count"], 0)
        self.assertEqual(family["annotation_only_count"], 3)

    def test_two_sequence_family_remains_reference_query_only(self):
        result = self.build_with_fake_tools([annotation_row("A1"), annotation_row("A2")])
        self.assertEqual(family_profile(result)["model_layer"], "reference_query_only")

    def test_experimental_positives_are_counted_as_anchors_in_the_discovery_layer(self):
        """The functional-calibration gate must no longer decide model construction."""
        result = self.build_with_fake_tools([row("P1"), row("P2"), row("P3")])
        family = family_profile(result)
        self.assertEqual(family["model_layer"], DISCOVERY_LAYER)
        self.assertEqual(family["training_sequence_count"], 3)
        self.assertEqual(family["experimental_anchor_count"], 3)
        self.assertEqual(family["annotation_only_count"], 0)
        self.assertEqual(family["functional_calibration_status"], "not_function_calibrated")

    def test_mixed_annotation_and_experimental_references_build_one_discovery_hmm(self):
        rows = [row("P1"), annotation_row("A1"), annotation_row("A2")]
        family = family_profile(self.build_with_fake_tools(rows))
        self.assertEqual(family["model_layer"], DISCOVERY_LAYER)
        self.assertEqual(family["experimental_anchor_count"], 1)
        self.assertEqual(family["annotation_only_count"], 2)
        self.assertEqual(family["training_accessions"], "A1;A2;P1")

    def test_duplicate_sequence_hashes_count_as_one_training_sequence(self):
        module = _load_module()
        records = [
            {"accession": "A3", "evidence_status": "annotation_only", "sequence_sha256": "a" * 64},
            {"accession": "A2", "evidence_status": "annotation_only", "sequence_sha256": "b" * 64},
            {"accession": "A1", "evidence_status": "annotation_only", "sequence_sha256": "b" * 64},
        ]
        self.assertEqual(len(module.select_discovery_training_rows(records)), 3)
        unique = module.unique_discovery_training_rows(records)
        self.assertEqual(len(unique), 2)
        self.assertEqual([item["accession"] for item in unique], ["A3", "A1"])

    def test_three_rows_sharing_two_hashes_are_below_the_training_threshold(self):
        module = _load_module()
        records = [
            {"accession": "A1", "evidence_status": "annotation_only", "sequence_sha256": "b" * 64},
            {"accession": "A2", "evidence_status": "annotation_only", "sequence_sha256": "b" * 64},
            {"accession": "A3", "evidence_status": "annotation_only", "sequence_sha256": "a" * 64},
        ]
        self.assertEqual(module.MINIMUM_DISCOVERY_TRAINING_SEQUENCES, 3)
        self.assertLess(
            len(module.unique_discovery_training_rows(records)),
            module.MINIMUM_DISCOVERY_TRAINING_SEQUENCES,
        )

    def test_training_selection_is_independent_of_ledger_order(self):
        module = _load_module()
        rows = [annotation_row("A2"), annotation_row("A1"), annotation_row("A3")]
        for index, item in enumerate(rows):
            item["sequence_sha256"] = hashlib.sha256(f"sequence-{index}".encode("ascii")).hexdigest()
        forward = [
            (item["sequence_sha256"], item["accession"])
            for item in module.select_discovery_training_rows(rows)
        ]
        backward = [
            (item["sequence_sha256"], item["accession"])
            for item in module.select_discovery_training_rows(list(reversed(rows)))
        ]
        self.assertEqual(forward, backward)
        self.assertEqual(forward, sorted(forward))

    def test_eligibility_falls_back_to_evidence_status_when_the_column_is_absent(self):
        module = _load_module()
        self.assertTrue(module.discovery_training_eligible(annotation_row("A1")))
        self.assertTrue(module.discovery_training_eligible(row("P1")))
        self.assertFalse(module.discovery_training_eligible(negative_row("N1")))
        self.assertFalse(module.discovery_training_eligible(
            {"accession": "C1", "evidence_status": "challenge_control"}
        ))
        self.assertFalse(module.discovery_training_eligible(
            {"accession": "R1", "evidence_status": "pending_review"}
        ))

    def test_explicit_eligibility_column_overrides_the_fallback(self):
        module = _load_module()
        declared_ineligible = annotation_row("A1")
        declared_ineligible[ELIGIBILITY_FIELD] = "false"
        self.assertFalse(module.discovery_training_eligible(declared_ineligible))
        declared_eligible = {"accession": "R1", "evidence_status": "pending_review", ELIGIBILITY_FIELD: "true"}
        self.assertTrue(module.discovery_training_eligible(declared_eligible))
        with self.assertRaisesRegex(ValueError, ELIGIBILITY_FIELD):
            module.discovery_training_eligible(
                {"accession": "X1", "evidence_status": "annotation_only", ELIGIBILITY_FIELD: "maybe"}
            )

    def test_declared_ineligibility_keeps_a_row_out_of_the_training_set(self):
        rows = [annotation_row("A1"), annotation_row("A2"), annotation_row("A3")]
        rows[0][ELIGIBILITY_FIELD] = "false"
        family = family_profile(self.build_with_fake_tools(rows, extra_fields=(ELIGIBILITY_FIELD,)))
        self.assertEqual(family["model_layer"], REFERENCE_LAYER)
        self.assertEqual(family["training_sequence_count"], 2)
        self.assertEqual(family["training_accessions"], "A2;A3")

    def test_eligible_row_without_sequence_hash_is_rejected(self):
        module = _load_module()
        ledger, fasta = write_inputs(self.root, [row("P1"), row("P2"), row("P3")])
        blank_sequence_sha256(ledger, "P2")
        with self.assertRaisesRegex(ValueError, "P2: discovery training eligibility requires a sequence_sha256"):
            module.build_profiles(ledger, fasta, self.root / "profiles-hashless", expected_family_count=None)

    def test_contradictory_negative_blocks_functional_calibration_and_stays_in_the_audit(self):
        rows = [annotation_row("A1"), annotation_row("A2"), annotation_row("A3"), negative_row("N1")]
        result = self.build_with_fake_tools(rows)
        family = family_profile(result)
        self.assertEqual(family["model_layer"], DISCOVERY_LAYER)
        self.assertEqual(family["functional_calibration_status"], BLOCKED_FUNCTIONAL_STATUS)
        self.assertNotIn("N1", family["training_accessions"])
        audit = [
            item for item in read_tsv(result["training_audit_path"])
            if item["profile_id"] == family["profile_id"]
        ]
        self.assertEqual({item["accession"] for item in audit}, {"A1", "A2", "A3", "N1"})
        negative = next(item for item in audit if item["accession"] == "N1")
        self.assertEqual(negative["training_selected"], "false")
        self.assertEqual(negative["exclusion_reason"], "evidence_status_not_discovery_training_eligible")
        self.assertEqual(negative["blocks_functional_calibration"], "true")
        self.assertEqual(negative["block_reason"], "contradictory_experimental_negative")

    def test_two_eligible_sequences_with_a_negative_are_reference_only_and_blocked(self):
        rows = [annotation_row("A1"), annotation_row("A2"), negative_row("N1")]
        family = family_profile(self.build_with_fake_tools(rows))
        self.assertEqual(family["model_layer"], REFERENCE_LAYER)
        self.assertEqual(family["model_status"], "reference_only")
        self.assertEqual(family["functional_calibration_status"], BLOCKED_FUNCTIONAL_STATUS)
        self.assertEqual(family["hmm_path"], "")

    def test_manifest_declares_only_builder_owned_layers_and_never_bit_reproducibility(self):
        rows = [
            annotation_row("A1"), annotation_row("A2"), annotation_row("A3"),
            annotation_row("B1", family="family-2"), annotation_row("B2", family="family-2"),
        ]
        result = self.build_with_fake_tools(rows)
        manifest = read_tsv(result["profile_manifest"])
        self.assertGreaterEqual(len(manifest), 3)
        for item in manifest:
            self.assertIn(item["model_layer"], {REFERENCE_LAYER, DISCOVERY_LAYER})
            self.assertNotIn(item["model_layer"], VALIDATED_LAYERS)
            self.assertEqual(item["bit_reproducible"], "false")
            self.assertIn(
                item["functional_calibration_status"],
                {"not_function_calibrated", BLOCKED_FUNCTIONAL_STATUS},
            )
        trained = [item for item in manifest if item["model_layer"] == DISCOVERY_LAYER]
        untrained = [item for item in manifest if item["model_layer"] == REFERENCE_LAYER]
        self.assertTrue(trained)
        self.assertTrue(untrained)
        for item in trained:
            self.assertNotEqual(item["training_sequence_count"], "0")
            self.assertEqual(len(item["training_set_sha256"]), 64)
            self.assertTrue(item["alignment_sha256"])
            self.assertTrue(item["hmm_sha256"])
        for item in untrained:
            self.assertEqual(item["training_set_sha256"], "pending")
            self.assertEqual(item["hmm_path"], "")
        self.assertEqual(family_profile(result)["training_sequence_count"], 3)
        self.assertEqual(superfamily_profile(result)["training_sequence_count"], 5)

    def test_training_set_sha256_binds_the_written_training_fasta(self):
        rows = [annotation_row("A2"), annotation_row("A1"), annotation_row("A3")]
        result = self.build_with_fake_tools(rows)
        family = family_profile(result)
        manifest_row = next(
            item for item in read_tsv(result["profile_manifest"])
            if item["profile_id"] == family["profile_id"]
        )
        training_path = Path(result["profile_manifest"]).parent / "profiles" / f"{family['profile_id']}.training.faa"
        self.assertTrue(training_path.is_file())
        self.assertEqual(
            manifest_row["training_set_sha256"],
            hashlib.sha256(training_path.read_bytes()).hexdigest(),
        )
        self.assertEqual(manifest_row["training_accessions"], "A1;A2;A3")
        headers = [
            line[1:].strip()
            for line in training_path.read_text(encoding="ascii").splitlines()
            if line.startswith(">")
        ]
        self.assertEqual(headers, manifest_row["training_accessions"].split(";"))

    def test_with_lipase_family_is_built_but_never_claims_functional_discrimination(self):
        rows = [
            annotation_row("A1", family="DED_hfam_2", superfamily=WITH_LIPASE_SUPERFAMILY),
            annotation_row("A2", family="DED_hfam_2", superfamily=WITH_LIPASE_SUPERFAMILY),
            annotation_row("A3", family="DED_hfam_2", superfamily=WITH_LIPASE_SUPERFAMILY),
        ]
        result = self.build_with_fake_tools(rows)
        family = family_profile(result)
        self.assertEqual(family["phaded_superfamily"], WITH_LIPASE_SUPERFAMILY)
        self.assertEqual(family["model_layer"], DISCOVERY_LAYER)
        self.assertEqual(family["functional_calibration_status"], "not_function_calibrated")
        hmm_path = Path(result["profile_manifest"]).parent / "profiles" / f"{family['profile_id']}.hmm"
        self.assertEqual(family["hmm_sha256"], hashlib.sha256(hmm_path.read_bytes()).hexdigest())

    def test_manifest_columns_declare_the_mandated_discovery_provenance_fields(self):
        module = _load_module()
        for field in (
            "model_layer", "functional_calibration_status", "training_sequence_count",
            "experimental_anchor_count", "annotation_only_count", "training_accessions",
            "training_set_sha256", "alignment_sha256", "hmm_sha256", "mafft_version",
            "hmmer_version", "bit_reproducible",
        ):
            self.assertIn(field, module.MANIFEST_COLUMNS)
        for field in (
            "profile_id", "profile_kind", "phaded_superfamily", "phaded_family_id",
            "training_count", "model_status", "model_reason", "alignment_path", "hmm_path",
        ):
            self.assertIn(field, module.MANIFEST_COLUMNS)
        self.assertEqual(len(module.MANIFEST_COLUMNS), len(set(module.MANIFEST_COLUMNS)))

    def test_manifest_has_one_training_audit_row_per_ledger_row_and_profile(self):
        rows = [annotation_row("A1"), annotation_row("A2"), annotation_row("A3"), negative_row("N1")]
        result = self.build_with_fake_tools(rows)
        audit = read_tsv(result["training_audit_path"])
        manifest = read_tsv(result["profile_manifest"])
        for item in manifest:
            self.assertEqual(
                {entry["accession"] for entry in audit if entry["profile_id"] == item["profile_id"]},
                {"A1", "A2", "A3", "N1"},
            )
        selected = {entry["accession"] for entry in audit if entry["training_selected"] == "true"}
        self.assertEqual(selected, {"A1", "A2", "A3"})
        layers = {item["profile_id"]: item["model_layer"] for item in manifest}
        for entry in audit:
            self.assertEqual(entry["model_layer"], layers[entry["profile_id"]])
            if entry["model_layer"] == DISCOVERY_LAYER:
                self.assertEqual(
                    entry["functional_calibration_status"], BLOCKED_FUNCTIONAL_STATUS
                )
            if entry["training_selected"] == "true":
                self.assertEqual(entry["exclusion_reason"], "")
            else:
                self.assertNotEqual(entry["exclusion_reason"], "")


if __name__ == "__main__":
    unittest.main()
