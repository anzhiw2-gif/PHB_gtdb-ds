#!/usr/bin/env python3
"""Pre-register the falsification gates of the Cys discovery-layer HMM
(2026-09-17, Step 4).

This script must be run **before** any training or scoring happens.  It fails
closed (``RuntimeError``) if a training artefact (``results/cys_discovery.hmm``)
or a scoring output (``results/cys_discovery_scores.tsv``) already exists in the
run directory, so a gate can never be registered after the fact.

It writes, into ``inputs/``:

* ``falsification_preregistration.tsv`` - the machine-readable gate table;
* ``falsification_preregistration.md`` - the verbatim pre-registration text;
* ``falsification_preregistration.json`` - timestamp + bound constants.

Bound constants: the discovery E-value threshold (``1e-5``), the HMMER ``-Z``
comparison count (``109087``, so that reference-set, confounder-set and
candidate-set E-values are on one scale), the model layer label
(``discovery_hmm_uncalibrated``), the calibration gate (declared unchanged) and
the SHA-256 of ``pipeline/config/formal_scan_models.tsv`` *before* the run.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import os
from pathlib import Path

MODEL_LAYER = "discovery_hmm_uncalibrated"
DISCOVERY_EVALUE_THRESHOLD = "1e-05"
HMMSEARCH_DATABASE_SIZE_Z = 109087
GATE_FIELDS = [
    "gate_id",
    "gate_name",
    "statement",
    "metric",
    "threshold",
    "expected_observation",
    "fallback_action",
    "population",
    "explained_cross_talk_population",
    "explained_cross_talk_disclosure_required",
    "explained_cross_talk_counted_as_gate_failure",
    "registered_at_utc",
    "registered_at_local",
    "registered_before_training",
    "model_layer",
    "threshold_change_permitted_after_failure",
]

_GATES: list[dict] = [
    {
        "gate_id": "G1",
        "gate_name": "specificity",
        "statement": (
            "At the pre-registered discovery threshold the discovery HMM must not hit any "
            "sequence of the extracellular reference panel and must not hit the MCL "
            "depolymerase control Q84C08."
        ),
        "metric": "unexplained_hits",
        "threshold": "unexplained_hits == 0",
        "expected_observation": "0 hits among Q84C08 + the 365 extracellular references",
        "fallback_action": (
            "report falsified; do not lower the threshold; the discovery layer must not be "
            "used to claim specificity"
        ),
        "population": (
            "Q84C08 (1) + extracellular references (365 = 284 extracellular dPHASCL type 1 "
            "+ 71 extracellular dPHASCL type 2 + 6 extracellular dPHAMCL + 4 extracellular "
            "native-SCL/PhaZ7-like)"
        ),
        "explained_cross_talk_population": (
            "29 intracellular nPHASCL with lipase box references; any hit there is recorded as "
            "explained_cross_talk and disclosed per record, and is NOT counted as a G1 failure"
        ),
        "explained_cross_talk_disclosure_required": "true",
        "explained_cross_talk_counted_as_gate_failure": "false",
    },
    {
        "gate_id": "G2",
        "gate_name": "sensitivity",
        "statement": (
            "Recall of the discovery HMM over the integrity-filter-passing Cys references must "
            "be at least 0.90. This is an IN-SAMPLE measurement (the retained references are the "
            "training alignment members); it is explicitly not a held-out generalisation "
            "estimate, because the frozen ledger contains no held-out Cys positive."
        ),
        "metric": "recall",
        "threshold": "recall >= 0.90",
        "expected_observation": ">= 0.90 of integrity-filter-passing Cys references are hit",
        "fallback_action": (
            "report falsified; the discovery layer must not be presented as family-discriminative"
        ),
        "population": (
            "the 276 Cys superfamily references restricted to those passing the integrity "
            "filter (localization / lipase-box annotation conflict / minimum core length)"
        ),
        "explained_cross_talk_population": "not applicable",
        "explained_cross_talk_disclosure_required": "false",
        "explained_cross_talk_counted_as_gate_failure": "false",
    },
    {
        "gate_id": "G3",
        "gate_name": "stratification",
        "statement": (
            "Every output of this run is labelled discovery_hmm_uncalibrated and "
            "pipeline/config/formal_scan_models.tsv is byte-identical before and after the run."
        ),
        "metric": "formal_scan_models_sha256_unchanged AND unlabelled_outputs == 0",
        "threshold": "sha256_before == sha256_after AND unlabelled_outputs == 0",
        "expected_observation": "identical SHA-256 and zero unlabelled outputs",
        "fallback_action": (
            "report falsified; any registry/config change invalidates the whole run"
        ),
        "population": "every file produced by the run",
        "explained_cross_talk_population": "not applicable",
        "explained_cross_talk_disclosure_required": "false",
        "explained_cross_talk_counted_as_gate_failure": "false",
    },
    {
        "gate_id": "G4",
        "gate_name": "no_family_call",
        "statement": (
            "The run must not produce a family call or change any subtype_call: no output table "
            "may carry a family_call / subtype_call column, and the number of subtype_call rows "
            "changed by this run must be 0."
        ),
        "metric": "subtype_call_rows_changed",
        "threshold": "subtype_call_rows_changed == 0",
        "expected_observation": "0 rows changed and no family/subtype columns in any output",
        "fallback_action": (
            "report falsified; any family judgement derived from these hits must be withdrawn"
        ),
        "population": "all 109,087 candidates",
        "explained_cross_talk_population": "not applicable",
        "explained_cross_talk_disclosure_required": "false",
        "explained_cross_talk_counted_as_gate_failure": "false",
    },
]


def gate_definitions() -> list[dict]:
    """Return the four pre-registered gates (G1-G4)."""
    return [dict(gate) for gate in _GATES]


def _forbidden_present(paths: list[Path]) -> list[str]:
    return [str(p) for p in paths if p.exists()]


def write_preregistration(run_dir: os.PathLike[str] | str, *, formal_scan_models_sha256: str) -> dict:
    run_dir = Path(run_dir)
    inputs = run_dir / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    results = run_dir / "results"
    results.mkdir(parents=True, exist_ok=True)

    forbidden = _forbidden_present([
        results / "cys_discovery.hmm",
        results / "cys_discovery_scores.tsv",
        results / "gate_evaluation.json",
    ])
    if forbidden:
        raise RuntimeError(
            "refusing to pre-register gates after training/scoring artefacts exist: "
            + ", ".join(forbidden)
        )

    now_utc = dt.datetime.now(dt.timezone.utc)
    registered_at_utc = now_utc.isoformat(timespec="seconds").replace("+00:00", "Z")
    registered_at_local = now_utc.astimezone().isoformat(timespec="seconds")

    gates = []
    for gate in gate_definitions():
        row = dict(gate)
        row.update({
            "registered_at_utc": registered_at_utc,
            "registered_at_local": registered_at_local,
            "registered_before_training": "true",
            "model_layer": MODEL_LAYER,
            "threshold_change_permitted_after_failure": "false",
        })
        gates.append(row)

    tsv_path = inputs / "falsification_preregistration.tsv"
    with tsv_path.open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.DictWriter(handle, fieldnames=GATE_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for gate in gates:
            writer.writerow({name: gate.get(name, "") for name in GATE_FIELDS})

    payload = {
        "schema": "phaded-cys-discovery-preregistration-v1",
        "run_id": run_dir.name,
        "model_layer": MODEL_LAYER,
        "registered_at_utc": registered_at_utc,
        "registered_at_local": registered_at_local,
        "registered_before_training": "true",
        "discovery_evalue_threshold": DISCOVERY_EVALUE_THRESHOLD,
        "hmmsearch_database_size_Z": HMMSEARCH_DATABASE_SIZE_Z,
        "database_size_Z_rationale": (
            "hmmsearch E-values scale with the comparison count Z; fixing -Z 109087 for the "
            "reference, confounder and candidate searches keeps one E-value scale"
        ),
        "threshold_lowering_permitted": False,
        "threshold_change_permitted_after_failure": False,
        "calibration_gate_unchanged": True,
        "calibration_gate_text": (
            "per family: >=3 independent positives AND >=1 held-out AND >=1 family-resolved "
            "negative AND >=1 challenge AND zero unexplained hits; NOT touched by this run and "
            "still NOT met"
        ),
        "formal_scan_models_tsv_sha256_before": formal_scan_models_sha256,
        "family_call_made": False,
        "registry_eligible": False,
        "training_equals_calibration": False,
        "gates": gates,
    }
    json_path = inputs / "falsification_preregistration.json"
    with json_path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")

    md_lines = [
        "# 预注册证伪闸门 —— Cys 型发现层 HMM（`%s`）" % MODEL_LAYER,
        "",
        "**预注册时间（UTC）：** `%s`" % registered_at_utc,
        "**预注册时间（本地）：** `%s`" % registered_at_local,
        "**注册先于训练与打分：** `true`（若训练或打分产物已存在，本脚本 fail-closed 拒绝写入）",
        "**发现阈值（钉死）：** E-value < 1e-5；**HMMER 比较数：** `-Z %d`（三个搜索集同尺度）"
        % HMMSEARCH_DATABASE_SIZE_Z,
        "**阈值不得在失败后下调：** `true`",
        "**校准 gate：** 未变（≥3 独立阳性 + ≥1 held-out + ≥1 family-resolved 阴性 + ≥1 challenge + 零未解释命中），且**仍未达成**",
        "**`pipeline/config/formal_scan_models.tsv` 运行前 SHA-256：** `%s`" % formal_scan_models_sha256,
        "",
    ]
    for gate in gates:
        md_lines += [
            "## %s（%s）" % (gate["gate_id"], gate["gate_name"]),
            "",
            "- **陈述：** %s" % gate["statement"],
            "- **指标：** `%s`" % gate["metric"],
            "- **阈值：** `%s`" % gate["threshold"],
            "- **期望观测：** %s" % gate["expected_observation"],
            "- **总体：** %s" % gate["population"],
            "- **已解释交叉（不计入失败）：** %s" % gate["explained_cross_talk_population"],
            "- **违反时的回退：** %s" % gate["fallback_action"],
            "",
        ]
    md_lines += [
        "---",
        "",
        "*本预注册为 candidate-only 记录；发现层命中只表示候选同源或功能潜力，不产生 family 判定，"
        "不是已验证的 PHB/PHA 降解表型，且不进 registry。*",
        "",
    ]
    md_path = inputs / "falsification_preregistration.md"
    md_path.write_text("\n".join(md_lines), encoding="utf-8", newline="\n")

    payload["written_files"] = [str(tsv_path), str(md_path), str(json_path)]
    return payload


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--formal-scan-models", required=True)
    args = parser.parse_args(argv)

    from run_context import sha256_file  # local import so tests can load the module standalone

    payload = write_preregistration(
        args.run_dir, formal_scan_models_sha256=sha256_file(args.formal_scan_models)
    )
    print(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
